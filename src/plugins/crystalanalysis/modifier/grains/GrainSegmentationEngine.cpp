///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2018) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <core/utilities/concurrent/ParallelFor.h>
#include <plugins/particles/util/NearestNeighborFinder.h>
#include "GrainSegmentationEngine.h"
#include "GrainSegmentationModifier.h"
#include "DisjointSet.h"
#include "Graph.h"

#include <ptm/ptm_functions.h>
#include <ptm/ptm_quat.h>


namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

/******************************************************************************
* Constructor.
******************************************************************************/
GrainSegmentationEngine::GrainSegmentationEngine(
			ParticleOrderingFingerprint fingerprint, ConstPropertyPtr positions, const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, ConstPropertyPtr selection,
			FloatType rmsdCutoff, FloatType mergingThreshold,
			int minGrainAtomCount, bool orphanAdoption, bool outputBonds) :
	StructureIdentificationModifier::StructureIdentificationEngine(std::move(fingerprint), positions, simCell, std::move(typesToIdentify), std::move(selection)),
	_rmsdCutoff(rmsdCutoff),
	_mergingThreshold(mergingThreshold),
	_minGrainAtomCount(std::max(minGrainAtomCount, 1)),
	_rmsd(std::make_shared<PropertyStorage>(positions->size(), PropertyStorage::Float, 1, 0, QStringLiteral("RMSD"), false)),
	_orientations(ParticlesObject::OOClass().createStandardStorage(positions->size(), ParticlesObject::OrientationProperty, true)),
	_atomClusters(ParticlesObject::OOClass().createStandardStorage(positions->size(), ParticlesObject::ClusterProperty, true)),
	_orphanAdoption(orphanAdoption),
	_outputBonds(outputBonds)
{
	// Allocate memory for neighbor lists.
	_neighborLists = std::make_shared<PropertyStorage>(positions->size(), PropertyStorage::Int64, PTM_MAX_NBRS, 0, QStringLiteral("Neighbors"), false);
	std::fill(_neighborLists->dataInt64(), _neighborLists->dataInt64() + _neighborLists->size() * _neighborLists->componentCount(), -1);
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine::perform()
{
	if(positions()->size() == 0)
		return;	// No input particles, nothing to do.

	// Grain segmentation algorithm:
	if(!identifyAtomicStructures()) return;
	if(!buildNeighborGraph()) return;
	if(!mergeSuperclusters()) return;
	if(!regionMerging()) return;
return;

	//if(!randomizeClusterIDs()) return;
	if(!calculateAverageClusterOrientations()) return;

	// For final output, convert edge disorientation angles from radians to degrees.
	if(_outputBonds) {
		for(FloatType& angle : neighborDisorientationAngles()->floatRange())
			angle *= FloatType(180) / FLOATTYPE_PI;
	}

	if (_orphanAdoption) {
		if(!mergeOrphanAtoms()) return;
	}
}

/******************************************************************************
* Performs the PTM algorithm. Determines the local structure type and the
* local lattice orientation at each atomic site.
******************************************************************************/
bool GrainSegmentationEngine::identifyAtomicStructures()
{
	// Initialize the PTMAlgorithm object.
	boost::optional<PTMAlgorithm> _algorithm;
	_algorithm.emplace();
	_algorithm->setRmsdCutoff(0.0); // Note: We do our own RMSD threshold filtering in postProcessStructureTypes().

	// Specify the structure types the PTM should look for.
	for(int i = 0; i < typesToIdentify().size() && i < PTMAlgorithm::NUM_STRUCTURE_TYPES; i++) {
		_algorithm->setStructureTypeIdentification(static_cast<PTMAlgorithm::StructureType>(i), typesToIdentify()[i]);
	}

	if(!_algorithm->prepare(*positions(), cell(), selection().get(), task().get()))
		return false;

	// Initialize the neighbor finder for disordered atoms
	// Don't need more than 12 neighbors for defect atoms.
	#define MAX_DISORDERED_NEIGHBORS 12
	NearestNeighborFinder neighFinder(MAX_DISORDERED_NEIGHBORS);
	if(!neighFinder.prepare(*positions(), cell(), selection().get(), task().get()))
		return false;

	task()->setProgressValue(0);
	task()->setProgressMaximum(positions()->size());
	task()->setProgressText(GrainSegmentationModifier::tr("Pre-calculating neighbor ordering"));

	// Pre-order neighbors of each particle
	std::vector< uint64_t > cachedNeighbors(positions()->size());
	parallelForChunks(positions()->size(), *task(), [this, &cachedNeighbors, &_algorithm](size_t startIndex, size_t count, Task& task) {
		// Create a thread-local kernel for the PTM algorithm.
		PTMAlgorithm::Kernel kernel(*_algorithm);

		// Loop over input particles.
		size_t endIndex = startIndex + count;
		for(size_t index = startIndex; index < endIndex; index++) {

			// Update progress indicator.
			if((index % 256) == 0)
				task.incrementProgressValue(256);

			// Break out of loop when operation was canceled.
			if(task.isCanceled())
				break;

			// Skip particles that are not included in the analysis.
			if(selection() && !selection()->getInt(index)) {
				structures()->setInt(index, PTMAlgorithm::OTHER);
				rmsd()->setFloat(index, 0);
				continue;
			}

			// Calculate ordering of neighbors
			kernel.precacheNeighbors(index, &cachedNeighbors[index]);
		}
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	task()->setProgressValue(0);
	task()->setProgressText(GrainSegmentationModifier::tr("Performing polyhedral template matching"));


//FILE* out = fopen("graphene_positions.txt", "w");

	// Perform analysis on each particle.
	parallelForChunks(positions()->size(), *task(), [this, &cachedNeighbors, &_algorithm, &neighFinder](size_t startIndex, size_t count, Task& task) {
		// Create a thread-local kernel for the PTM algorithm.
		PTMAlgorithm::Kernel kernel(*_algorithm);

		// Create a neighbor query for disordered atoms
		NearestNeighborFinder::Query<MAX_DISORDERED_NEIGHBORS> neighQuery(neighFinder);

//size_t startIndex = 0, count = positions()->size();
		// Loop over input particles.
		size_t endIndex = startIndex + count;
		for(size_t index = startIndex; index < endIndex; index++) {

			// Update progress indicator.
			if((index % 256) == 0)
				task.incrementProgressValue(256);

			// Break out of loop when operation was canceled.
			if(task.isCanceled())
				break;

			// Skip particles that are not included in the analysis.
			if(selection() && !selection()->getInt(index)) {
				structures()->setInt(index, PTMAlgorithm::OTHER);
				rmsd()->setFloat(index, 0);
				continue;
			}

			// Perform the PTM analysis for the current particle.
			PTMAlgorithm::StructureType type = kernel.identifyStructure(index, cachedNeighbors, nullptr);

			// Store results in the output arrays.
			structures()->setInt(index, type);
			rmsd()->setFloat(index, kernel.rmsd());

//fprintf(out, "%f %f %f %f %f %f %f %f\n", positions()->getPoint3(index).x(), positions()->getPoint3(index).y(), positions()->getPoint3(index).z(),
//					kernel.rmsd(), kernel.orientation().w(), kernel.orientation().x(), kernel.orientation().y(), kernel.orientation().z());
//fflush(out);

			if(type != PTMAlgorithm::OTHER) {
				if(orientations()) orientations()->setQuaternion(index, kernel.orientation());

				// Store neighbor list.
				for(int j = 0; j < ptm_num_nbrs[type]; j++) {
					//OVITO_ASSERT(j < _neighborLists->componentCount());
					//OVITO_ASSERT(mapping[j + 1] >= 1);
					//OVITO_ASSERT(mapping[j + 1] <= numNeighbors);

					_neighborLists->setInt64Component(index, j, kernel._env.atom_indices[j + 1]);

					// Check if neighbor vector spans more than half of a periodic simulation cell.
					double* delta = kernel._env.points[j + 1];
					Vector3 neighborVector(delta[0], delta[1], delta[2]);
					for(size_t dim = 0; dim < 3; dim++) {
						if(cell().pbcFlags()[dim]) {
							if(std::abs(cell().inverseMatrix().prodrow(neighborVector, dim)) >= FloatType(0.5)+FLOATTYPE_EPSILON) {
								static const QString axes[3] = { QStringLiteral("X"), QStringLiteral("Y"), QStringLiteral("Z") };
								throw Exception(GrainSegmentationModifier::tr("Simulation box is too short along cell vector %1 (%2) to perform analysis. "
										"Please extend it first using the 'Replicate' modifier.").arg(dim+1).arg(axes[dim]));
							}
						}
					}
				}
			}
			else {
				rmsd()->setFloat(index, 0);

				neighQuery.findNeighbors(index);
				int numNeighbors = neighQuery.results().size();

				// Store neighbor list.
				for(int j = 0; j < numNeighbors; j++) {
					_neighborLists->setInt64Component(index, j, neighQuery.results()[j].index);
				}
			}
		}
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	// Determine histogram bin size based on maximum RMSD value.
	const size_t numHistogramBins = 100;
	_rmsdHistogram = std::make_shared<PropertyStorage>(numHistogramBins, PropertyStorage::Int64, 1, 0, GrainSegmentationModifier::tr("Count"), true, DataSeriesObject::YProperty);
	FloatType rmsdHistogramBinSize = FloatType(1.01) * *std::max_element(rmsd()->constDataFloat(), rmsd()->constDataFloat() + rmsd()->size()) / numHistogramBins;
	if(rmsdHistogramBinSize <= 0) rmsdHistogramBinSize = 1;
	_rmsdHistogramRange = rmsdHistogramBinSize * numHistogramBins;

	// Perform binning of RMSD values.
	for(size_t index = 0; index < structures()->size(); index++) {
		if(structures()->getInt(index) != PTMAlgorithm::OTHER) {
			OVITO_ASSERT(rmsd()->getFloat(index) >= 0);
			int binIndex = rmsd()->getFloat(index) / rmsdHistogramBinSize;
			if(binIndex < numHistogramBins)
				_rmsdHistogram->dataInt64()[binIndex]++;
		}
	}

	// Set the atomic structure to OTHER for atoms with an RMSD value above the threshold.
	if(_rmsdCutoff > 0) {
		for(size_t index = 0; index < structures()->size(); index++) {
			if(structures()->getInt(index) != PTMAlgorithm::OTHER) {
				if(rmsd()->getFloat(index) > _rmsdCutoff)
					structures()->setInt(index, PTMAlgorithm::OTHER);
			}
		}
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Builds the graph of neighbor atoms and calculates the misorientation angle
* for each graph edge (i.e. bond)
******************************************************************************/
bool GrainSegmentationEngine::buildNeighborGraph()
{
	// Generate bonds (edges) between neighboring lattice atoms.
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - edge generation"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(positions()->size());

	// Create the bonds connecting neighboring lattice atoms.
	for(size_t index = 0; index < positions()->size(); index++) {
		if(!task()->incrementProgressValue()) return false;

		for(size_t c = 0; c < _neighborLists->componentCount(); c++) {
			auto neighborIndex = _neighborLists->getInt64Component(index, c);
			if(neighborIndex == -1) break;

			// Skip every other atom pair so that we don't create the same bond twice.
			if(_neighborLists->getInt64Component(index, c) < index) continue;

			Bond bond;
			bond.index1 = index;
			bond.index2 = neighborIndex;
			// Determine PBC bond shift using minimum image convention.
			Vector3 delta = positions()->getPoint3(index) - positions()->getPoint3(neighborIndex);
			for(size_t dim = 0; dim < 3; dim++) {
				if(cell().pbcFlags()[dim])
					bond.pbcShift[dim] = (int)std::floor(cell().inverseMatrix().prodrow(delta, dim) + FloatType(0.5));
				else
					bond.pbcShift[dim] = 0;
			}

			// Skip every other bond to create only one bond per particle pair.
			if(!bond.isOdd())
				_latticeNeighborBonds.push_back(bond);
		}
	}

	// Compute disorientation angles associated with the neighbor graph edges.
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - misorientation calculation"));
	_neighborDisorientationAngles = std::make_shared<PropertyStorage>(_latticeNeighborBonds.size(), PropertyStorage::Float, 1, 0, QStringLiteral("Disorientation"), false);
	parallelFor(_latticeNeighborBonds.size(), *task(), [this](size_t bondIndex) {
		size_t index1 = latticeNeighborBonds()[bondIndex].index1;
		size_t index2 = latticeNeighborBonds()[bondIndex].index2;
		FloatType& disorientationAngle = *(neighborDisorientationAngles()->dataFloat() + bondIndex);
		disorientationAngle = std::numeric_limits<FloatType>::infinity();

		int structureTypeA = structures()->getInt(index1);
		int structureTypeB = structures()->getInt(index2);
		if(structureTypeA == structureTypeB) {

			int structureType = structureTypeA;
			const Quaternion& qA = orientations()->getQuaternion(index1);
			const Quaternion& qB = orientations()->getQuaternion(index2);

			double orientA[4] = { qA.w(), qA.x(), qA.y(), qA.z() };
			double orientB[4] = { qB.w(), qB.x(), qB.y(), qB.z() };

			if(structureType == PTMAlgorithm::SC || structureType == PTMAlgorithm::FCC || structureType == PTMAlgorithm::BCC || structureType == PTMAlgorithm::CUBIC_DIAMOND)
				disorientationAngle = (FloatType)ptm::quat_disorientation_cubic(orientA, orientB);
			else if(structureType == PTMAlgorithm::HCP || structureType == PTMAlgorithm::HEX_DIAMOND || structureType == PTMAlgorithm::GRAPHENE)
				disorientationAngle = (FloatType)ptm::quat_disorientation_hcp_conventional(orientA, orientB);
		}
	});

	return !task()->isCanceled();
}

/******************************************************************************
* Merges adjacent clusters with similar lattice orientations.
******************************************************************************/
bool GrainSegmentationEngine::mergeSuperclusters()
{
	size_t numAtoms = positions()->size();
	_superclusterSizes.resize(numAtoms, 1);
	_atomSuperclusters.resize(numAtoms, 0);
	DisjointSet uf(numAtoms);

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - supercluster merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(latticeNeighborBonds().size());

	FloatType threshold = _misorientationThreshold;

	// Merge superclusters.
	for(size_t bondIndex = 0; bondIndex < latticeNeighborBonds().size(); bondIndex++) {
		if(!task()->incrementProgressValue()) return false;

		const Bond& bond = latticeNeighborBonds()[bondIndex];
		// Skip high-angle edges.
		if(neighborDisorientationAngles()->getFloat(bondIndex) > _misorientationThreshold) continue;

		size_t parentA = uf.find(bond.index1);
		size_t parentB = uf.find(bond.index2);
		if(parentA != parentB) {
			uf.merge(parentA, parentB);
		}
	}
	if(task()->isCanceled())
		return false;

	for (size_t i=0;i<numAtoms;i++) {
		_atomSuperclusters[i] = uf.find(i);
		_superclusterSizes[i] = uf.sizes[i];
	}

	// Relabels the superclusters to obtain a contiguous sequence of cluster IDs.
	std::vector<size_t> superclusterRemapping(numAtoms);
	_numSuperclusters = 1;
	// Assign new consecutive IDs to root superclusters.
	for(size_t i = 0; i < numAtoms; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(_superclusterSizes[i] < _minGrainAtomCount) {
				superclusterRemapping[i] = 0;
			}
			else {
				superclusterRemapping[i] = _numSuperclusters;
				_superclusterSizes[_numSuperclusters] = _superclusterSizes[i];
				_numSuperclusters++;
			}
		}
	}

	// Determine new IDs for non-root superclusters.
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++)
		superclusterRemapping[particleIndex] = superclusterRemapping[uf.find(particleIndex)];

	// Relabel atoms after cluster IDs have changed.
	_superclusterSizes.resize(_numSuperclusters);
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++)
		_atomSuperclusters[particleIndex] = superclusterRemapping[particleIndex];

	// Supercluster 0 contains all atoms that are not part of a regular supercluster.
	_superclusterSizes[0] = 0;
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++)
		if(_atomSuperclusters[particleIndex] == 0)
			_superclusterSizes[0]++;

	return !task()->isCanceled();
}

class GraphEdge
{
public:
	GraphEdge(size_t _a, size_t _b, FloatType _w, size_t _superCluster)
		: a(_a), b(_b), w(_w), superCluster(_superCluster) {}

	size_t a;
	size_t b;
	FloatType w;
	size_t superCluster;
};

class DendrogramNode
{
public:
	DendrogramNode(size_t _a, size_t _b, FloatType _d, size_t _size)
		: a(_a), b(_b), d(_d), size(_size) {}

	size_t a;
	size_t b;
	FloatType d;
	size_t size;
};

/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine::regionMerging()
{
	size_t numAtoms = positions()->size();

	if (_numSuperclusters == 1) {
		_numClusters = 1;
		_clusterSizes.resize(_numClusters, 0);
		_clusterSizes[0] = _superclusterSizes[0];
		for (size_t particleIndex=0;particleIndex<numAtoms;particleIndex++) {
			atomClusters()->setInt64(particleIndex, _atomSuperclusters[particleIndex]);
		}
		return true;
	}

	_clusterSizes.resize(numAtoms, 1);
	DisjointSet uf(numAtoms);

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(latticeNeighborBonds().size());

	// Calculate a contiguous atom index mapping
	std::vector< std::tuple< size_t, size_t > > atomIntervals(_numSuperclusters);
	{
		size_t start = 0;
		for (size_t i=0;i<_numSuperclusters;i++) {
			std::get<0>(atomIntervals[i]) = start;
			std::get<1>(atomIntervals[i]) = start + _superclusterSizes[i];
			start += _superclusterSizes[i];
		}
	}

	if(task()->isCanceled())
		return false;

	// Calculate a contiguous bond index mapping
	std::vector< GraphEdge > initial_graph;

	// Build initial graph
	const FloatType* disorientation = neighborDisorientationAngles()->constDataFloat();
	for(const Bond& bond : latticeNeighborBonds()) {
		if(!task()->incrementProgressValue()) return false;

		FloatType dis = *disorientation++;

		// Skip high-angle edges.
		if(dis > _misorientationThreshold) continue;

		// Convert disorientations to graph weights
		FloatType deg = dis * FloatType(180) / FLOATTYPE_PI;
		FloatType weight = std::exp(-FloatType(1)/3 * deg * deg);		//this is fairly arbitrary but it works well

		size_t a = bond.index1;
		size_t b = bond.index2;
		size_t sc = _atomSuperclusters[a];
		initial_graph.emplace_back(a, b, weight, sc);
	}

	if(task()->isCanceled())
		return false;

	std::sort(initial_graph.begin(), initial_graph.end(), [](GraphEdge& e, GraphEdge& f) { return e.superCluster < f.superCluster; });
	std::map< size_t, std::tuple< size_t, size_t, int > > bondIntervals;

	size_t start = 0;
	for (size_t i=0;i<initial_graph.size();i++) {
		size_t sc = initial_graph[start].superCluster;
		if (initial_graph[i].superCluster != sc) {

			int type = structures()->getInt(initial_graph[i].a);
			bondIntervals[sc] = std::make_tuple(start, i, type);
			start = i;
		}
	}

	{
		size_t sc = initial_graph.back().superCluster;
		int type = structures()->getInt(initial_graph.back().a);
		bondIntervals[sc] = std::make_tuple(start, initial_graph.size(), type);
	}

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(initial_graph.size());

	std::vector< DendrogramNode > dendrogram;			// dendrogram as list of merges
	FloatType totalWeight = 0;
	for (auto edge : initial_graph) {
		totalWeight += edge.w;
	}

	// TODO: parallelize this. Use atomic push_back on dendrogram.
	//parallelFor(_numSuperclusters, *task(), [this, &numAtoms, &bondIntervals, &initial_graph, &uf](size_t sc) {
	//	if (sc == 0) return;
	for (size_t sc=1;sc<_numSuperclusters;sc++) {
		std::map< size_t, std::tuple< size_t, size_t, int > >::iterator it = bondIntervals.find(sc);
		auto value = it->second;
		size_t start = std::get<0>(value);
		size_t end = std::get<1>(value);
		int structureType = std::get<2>(value);

		Graph graph;
		for (size_t i=start;i<end;i++) {
			auto edge = initial_graph[i];
			graph.add_edge(edge.a, edge.b, edge.w, true);
		}
		graph.wtotal = totalWeight;

		//std::vector< std::tuple< size_t, size_t > > components;	// connected components

		size_t n = graph.num_nodes();
		while (graph.num_nodes()) {

			// nearest-neighbor chain
			size_t node = graph.next_node();
			if (node == (size_t)(-1)) {
				printf("node is -1\n");
				exit(3);
			}

			std::vector< size_t> chain{node};
			while (chain.size()) {

				size_t a = chain.back();
				chain.pop_back();
				if (a == (size_t)(-1)) {
					printf("a is -1\n");
					exit(3);
				}

				auto result = graph.nearest_neighbor(a);
				FloatType d = std::get<0>(result);
				size_t b = std::get<1>(result);
				if (b == (size_t)(-1)) {
					if (chain.size() != 0) {
						printf("non-zero chain size!\n"); fflush(stdout);
						//exit(3);
					}
				}

				if (b == (size_t)(-1)) {
					// remove the connected component
					size_t sa = graph.snode[a];
					//components.push_back(std::make_tuple(graph.rep[a], sa));
					graph.remove_node(a);
				}
				else if (chain.size()) {
					size_t c = chain.back();
					chain.pop_back();

					if (b == c) {
						size_t size = graph.snode[a] + graph.snode[b];
						//dendrogram.push_back(	DendrogramNode(	std::min(graph.rep[a], graph.rep[b]),
						//					std::max(graph.rep[a], graph.rep[b]),
						//					d,
						//					size)
						dendrogram.push_back(	DendrogramNode(	std::min(a, b),
											std::max(a, b),
											d,
											size)
									);
						if (size == 0) {
							printf("zero size\n");
							exit(3);
						}
						graph.contract_edge(a, b);
					}
					else {
						chain.push_back(c);
						chain.push_back(a);
						chain.push_back(b);
						if (a == (size_t)(-1)) {printf("!a is -1\n"); exit(3);}
						if (b == (size_t)(-1)) {printf("!b is -1\n"); exit(3);}
						if (c == (size_t)(-1)) {printf("!c is -1\n"); exit(3);}
					}
				}
				else if (b != (size_t)(-1)) {
					chain.push_back(a);
					chain.push_back(b);
					if (a == (size_t)(-1)) {printf("#a is -1\n"); exit(3);}
					if (b == (size_t)(-1)) {printf("#b is -1\n"); exit(3);}
				}
			}
		}

		// add connected components to the dendrogram
		//size_t u = graph.next;
		//a, s = components.pop()
		//for b, t in components:
		//	s += t
		//	D.append([min(a, b), max(a, b), float("inf"), s])
		//	a = u
		//	u += 1
	}

	std::sort(dendrogram.begin(), dendrogram.end(),
			[](DendrogramNode& a, DendrogramNode& b)
			{return a.d < b.d;});

	// Create PropertyStorage output objects
	_mergeDistance = std::make_shared<PropertyStorage>(dendrogram.size(), PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), true, DataSeriesObject::XProperty);
	_mergeSize = std::make_shared<PropertyStorage>(dendrogram.size(), PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Merge size"), true, DataSeriesObject::YProperty);

	// Scan through the entire merge list to determine merge sizes
	for (size_t i=0;i<dendrogram.size();i++) {
		auto node = dendrogram[i];
		size_t sa = uf.sizes[uf.find(node.a)];
		size_t sb = uf.sizes[uf.find(node.b)];
		uf.merge(node.a, node.b);
		size_t dsize = std::min(sa, sb);

		// output the data
		_mergeDistance->dataFloat()[i] = log(node.d);
		_mergeSize->dataFloat()[i] = dsize;
	}

	uf.clear();

	// Now iterate through merge list until distance cutoff is met
	for (auto node : dendrogram) {
		if (log(node.d) >= _mergingThreshold) {
			break;
		}

		uf.merge(node.a, node.b);
	}

	if(task()->isCanceled())
		return false;

	for(size_t i = 0; i < numAtoms; i++) {
		_clusterSizes[i] = uf.sizes[i];
	}

	// Relabels the clusters to obtain a contiguous sequence of cluster IDs.
	std::vector<size_t> clusterRemapping(numAtoms);

	// Assign new consecutive IDs to root clusters.
	_numClusters = 1;
	for(size_t i = 0; i < numAtoms; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(_clusterSizes[i] < _minGrainAtomCount) {
				clusterRemapping[i] = 0;
			}
			else {
				clusterRemapping[i] = _numClusters;
				_clusterSizes[_numClusters] = _clusterSizes[i];
				_numClusters++;
			}
		}
	}

	// Determine new IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++)
		clusterRemapping[particleIndex] = clusterRemapping[uf.find(particleIndex)];

	// Relabel atoms after cluster IDs have changed.
	_clusterSizes.resize(_numClusters);
	std::fill(_clusterSizes.begin(), _clusterSizes.end(), 0);
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {

		size_t gid = clusterRemapping[particleIndex];
		atomClusters()->setInt64(particleIndex, gid);
		_clusterSizes[gid]++;
	}

	// Relabel clusters by size (large to small)
	if (_numClusters > 1) {
		std::vector< std::tuple< size_t, size_t > > lut;

		for (size_t i=0;i<_numClusters;i++) {
			lut.emplace_back(std::make_tuple(i, _clusterSizes[i]));
		}

		// Sort by size, leaving zeroth cluster in place
		std::sort(lut.begin() + 1, lut.end(),
			[](std::tuple< size_t, size_t > a, std::tuple< size_t, size_t > b)
			{return std::get<1>(b) < std::get<1>(a);});

		std::vector< size_t > indices(_numClusters);
		for (size_t gid=0;gid<_numClusters;gid++) {
			indices[std::get<0>(lut[gid])] = gid;
		}

		for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {

			size_t gid = atomClusters()->getInt64(particleIndex);
			atomClusters()->setInt64(particleIndex, indices[gid]);
		}

		for (size_t i=0;i<_numClusters;i++) {
			_clusterSizes[i] = std::get<1>(lut[i]);
		}
	}

//TODO: output grain sizes as a property
printf("grain sizes:\n");
for (int i=0;i<_numClusters;i++) {
	printf("\t%d\t%lld\n", i, _clusterSizes[i]);
}

	_clusterOrientations.resize(_numClusters);
	return !task()->isCanceled();
}

/******************************************************************************
* Randomize cluster IDs for testing purposes (giving more color contrast).
******************************************************************************/
bool GrainSegmentationEngine::randomizeClusterIDs()
{
	size_t numAtoms = positions()->size();
	int numClusters = _numClusters;

	std::vector<size_t> clusterRandomMapping(numClusters+1);
	std::iota(clusterRandomMapping.begin(), clusterRandomMapping.end(), 0);
	std::mt19937 rng(1);
	std::shuffle(clusterRandomMapping.begin() + 1, clusterRandomMapping.end(), rng);

	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {

		int clusterId = atomClusters()->getInt64(particleIndex);
		int randomizedId = clusterRandomMapping[clusterId];
		atomClusters()->setInt64(particleIndex, randomizedId);
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Merges any orphan atoms into the closest cluster.
******************************************************************************/
bool GrainSegmentationEngine::mergeOrphanAtoms()
{
	// Build list of orphan atoms.
	std::vector<size_t> orphanAtoms;
	for(size_t i = 0; i < positions()->size(); i++) {
		if(atomClusters()->getInt64(i) == 0)
			orphanAtoms.push_back(i);
	}

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - merging orphan atoms"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(orphanAtoms.size());

	// Add orphan atoms to the grains.
	size_t oldOrphanCount = orphanAtoms.size();
	for(;;) {
		std::vector<size_t> newlyAssignedClusters(orphanAtoms.size(), 0);
		for(size_t i = 0; i < orphanAtoms.size(); i++) {
			if(task()->isCanceled()) return false;

			// Find the closest cluster atom in the neighborhood.
			FloatType minDistSq = FLOATTYPE_MAX;
			for(size_t c = 0; c < _neighborLists->componentCount(); c++) {
				auto neighborIndex = _neighborLists->getInt64Component(orphanAtoms[i], c);
				if(neighborIndex == -1) break;
				auto clusterId = atomClusters()->getInt64(neighborIndex);
				if(clusterId == 0) continue;

				// Determine interatomic vector using minimum image convention.
				Vector3 delta = cell().wrapVector(positions()->getPoint3(neighborIndex) - positions()->getPoint3(orphanAtoms[i]));
				FloatType distSq = delta.squaredLength();
				if(distSq < minDistSq) {
					minDistSq = distSq;
					newlyAssignedClusters[i] = clusterId;
				}
			}
		}

		// Assign atoms to closest cluster and compress orphan list.
		size_t newOrphanCount = 0;
		for(size_t i = 0; i < orphanAtoms.size(); i++) {
			atomClusters()->setInt64(orphanAtoms[i], newlyAssignedClusters[i]);
			if(newlyAssignedClusters[i] == 0) {
				orphanAtoms[newOrphanCount++] = orphanAtoms[i];
			}
			else {
				_clusterSizes[newlyAssignedClusters[i] - 1]++;
				if(!task()->incrementProgressValue()) return false;
			}
		}
		orphanAtoms.resize(newOrphanCount);
		if(newOrphanCount == oldOrphanCount)
			break;
		oldOrphanCount = newOrphanCount;
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Computes the average lattice orientation of each cluster.
******************************************************************************/
bool GrainSegmentationEngine::calculateAverageClusterOrientations()
{
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - average cluster orientation"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(positions()->size());

	// Allocate cluster orientation and size arrays.
	// We will create one cluster for each basin in the distance transform field.
	_clusterOrientations.resize(_numClusters, Quaternion(0,0,0,0));
	_clusterSizes.resize(_numClusters, 0);

	// Stores the seed atom index of each cluster.
	std::vector<qlonglong> firstClusterAtom(_numClusters, -1);

	for(size_t particleIndex = 0; particleIndex < positions()->size(); particleIndex++) {
		if(!task()->incrementProgressValue()) return false;

		qlonglong clusterId = atomClusters()->getInt64(particleIndex);
		if(clusterId == 0) continue;

		// Cluster IDs start at 1. Need to subtract 1 to get cluster index.
		qlonglong clusterIndex = clusterId - 1;

		_clusterSizes[clusterIndex]++;
		if(firstClusterAtom[clusterIndex] == -1)
			firstClusterAtom[clusterIndex] = particleIndex;

		const Quaternion& orient0 = orientations()->getQuaternion(firstClusterAtom[clusterIndex]);
		const Quaternion& orient = orientations()->getQuaternion(particleIndex);

		Quaternion qrot = orient0.inverse() * orient;
		double qrot_[4] = { qrot.w(), qrot.x(), qrot.y(), qrot.z() };

		int structureType = structures()->getInt(particleIndex);
		if(structureType == PTMAlgorithm::SC || structureType == PTMAlgorithm::FCC || structureType == PTMAlgorithm::BCC || structureType == PTMAlgorithm::CUBIC_DIAMOND)
			ptm::rotate_quaternion_into_cubic_fundamental_zone(qrot_);
		else if(structureType == PTMAlgorithm::HCP || structureType == PTMAlgorithm::HEX_DIAMOND || structureType == PTMAlgorithm::GRAPHENE)
			ptm::rotate_quaternion_into_hcp_conventional_fundamental_zone(qrot_);

		Quaternion qclosest = orient0 * Quaternion(qrot_[1], qrot_[2], qrot_[3], qrot_[0]);
		_clusterOrientations[clusterIndex] += qclosest;
	}
	for(auto& qavg : _clusterOrientations) {
		if(qavg.dot(qavg) > FLOATTYPE_EPSILON)
			qavg.normalize();
	}

	return !task()->isCanceled();
}

}	// End of namespace
}	// End of namespace
}	// End of namespace
