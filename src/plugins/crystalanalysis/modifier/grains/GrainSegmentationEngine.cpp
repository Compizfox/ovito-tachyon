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

			// Calculate ordering of neighbors
			kernel.precacheNeighbors(index, &cachedNeighbors[index]);
		}
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	task()->setProgressValue(0);
	task()->setProgressText(GrainSegmentationModifier::tr("Performing polyhedral template matching"));


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

class DisjointSet
{
public:
	DisjointSet(size_t n)
	{
		ranks.resize(n);

		parents.resize(n);
		std::iota(parents.begin(), parents.end(), (size_t)0);

		sizes.resize(n);
		std::fill(sizes.begin(), sizes.end(), 1);

		weights.resize(n);
		std::fill(weights.begin(), weights.end(), 0);
	}

	// "Find" part of Union-Find.
	size_t find(size_t index) {

		// Find root and make root as parent of i (path compression)
		size_t parent = parents[index];
		while(parent != parents[parent]) {
			parent = parents[parent];
		}

		parents[index] = parent;
		return parent;
	}

	// "Union" part of Union-Find.
	void merge(size_t index1, size_t index2) {
		size_t parentA = find(index1);
		size_t parentB = find(index2);
		if(parentA == parentB) return;

		// Attach smaller rank tree under root of high rank tree (Union by Rank)
		if(ranks[parentA] < ranks[parentB]) {
			parents[parentA] = parentB;
			sizes[parentB] += sizes[parentA];
			weights[parentB] += weights[parentA];
		}
		else {
			parents[parentB] = parentA;
			sizes[parentA] += sizes[parentB];
			weights[parentA] += weights[parentB];

			// If ranks are same, then make one as root and increment its rank by one
			if(ranks[parentA] == ranks[parentB])
				ranks[parentA]++;
		}
	}

	std::vector<size_t> sizes;
	std::vector<FloatType> weights;

private:
	std::vector<size_t> parents;
	std::vector<size_t> ranks;
};

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

printf("supercluster sizes:\n");
for(size_t i = 0; i < _numSuperclusters; i++)
	printf("\t%lu\t%lu\n", i, _superclusterSizes[i]);

	return !task()->isCanceled();
}


class GraphEdge
{
public:
	GraphEdge(size_t _a, size_t _b, FloatType _w, FloatType _mergeQuality, size_t _superCluster)
		: a(_a), b(_b), w(_w), mergeQuality(_mergeQuality), superCluster(_superCluster) {}

	size_t a;
	size_t b;
	FloatType w;
	FloatType mergeQuality;
	size_t superCluster;
};

static FloatType fastPower(FloatType x, int dim) {

	// Computes:	x^(1/2) if dim == 2
	//		x^(2/3) if dim == 3
	if (dim == 2) {
		return std::sqrt(x);
	}
	else if (dim == 3) {
		return std::cbrt(x * x);
	}
}

static FloatType mergeQuality(int type, FloatType curWeightA, FloatType curWeightB, FloatType w) {

	//TODO: use number of surface bonds in Wullf construction to calculate a good normalization factor, so threshold values are similar across structures.
	const int coordinations[PTMAlgorithm::NUM_STRUCTURE_TYPES] = {0, 12, 12, 14, 12, 6, 16, 16, 9};
	const int dimensions[PTMAlgorithm::NUM_STRUCTURE_TYPES] = {0, 3, 3, 3, 3, 3, 3, 3, 2};

	int dim = dimensions[type];
	FloatType hc = coordinations[type] / 2.0;
	FloatType k = 0.5;

	FloatType qualityA = (k * hc + w) / (1 + fastPower(curWeightA / hc, dim)) / coordinations[type];
	FloatType qualityB = (k * hc + w) / (1 + fastPower(curWeightB / hc, dim)) / coordinations[type];
	return std::max(qualityA, qualityB);
}


/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine::regionMerging()
{
	size_t numAtoms = positions()->size();
	_clusterSizes.resize(numAtoms, 1);
	DisjointSet uf(numAtoms);

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(latticeNeighborBonds().size());

	std::vector< GraphEdge > initial_graph;

FloatType mw = 10000;

	// Build initial graph
	const FloatType* disorientation = neighborDisorientationAngles()->constDataFloat();
	for(const Bond& bond : latticeNeighborBonds()) {
		if(!task()->incrementProgressValue()) return false;

		FloatType dis = *disorientation++;

		// Skip high-angle edges.
		if(dis > _misorientationThreshold) continue;

		// Convert disorientations to graph weights
		FloatType deg = dis * FloatType(180) / FLOATTYPE_PI;
		FloatType weight = std::exp(-FloatType(1)/3. * deg * deg);		//this is fairly arbitrary but it works well

		size_t sc = _atomSuperclusters[bond.index1];
		initial_graph.emplace_back(bond.index1, bond.index2, weight, -1, sc);
	}

	if(task()->isCanceled())
		return false;

	std::sort(initial_graph.begin(), initial_graph.end(), [](GraphEdge& e, GraphEdge& f) { return e.superCluster < f.superCluster; });

	std::vector< std::tuple< size_t, size_t > > intervals;
	size_t start = 0;
	for (size_t i=0;i<initial_graph.size();i++) {
		if (initial_graph[i].superCluster != initial_graph[start].superCluster) {
			intervals.emplace_back(std::make_tuple(start, i));
			start = i;
		}
	}
	intervals.emplace_back(std::make_tuple(start, initial_graph.size()));

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(initial_graph.size());

clock_t total_time[9] = {0};

	parallelFor(_numSuperclusters, *task(), [this, &numAtoms, &initial_graph, &intervals, &uf, &total_time](size_t sc) {
		if (sc == 0) return;

		std::vector< GraphEdge > graph;
		size_t start = std::get<0>(intervals[sc]);
		size_t end = std::get<1>(intervals[sc]);

		for (size_t i=start;i<end;i++) {
			graph.push_back(initial_graph[i]);
		}

		std::vector<bool> hit(numAtoms, 0);

		int iteration = 0;
		bool merged = true;
		while (merged) {
			merged = false;
			printf("iteration: %d %lu\n", iteration++, graph.size());

clock_t time[9];
time[0] = clock();
			std::vector<FloatType> curWeights(uf.weights.begin(), uf.weights.end());
time[1] = clock();

			for (size_t i=0;i<graph.size();i++) {	//can't use auto since we want to set object elements

				auto edge = graph[i];

				size_t a = edge.a;
				size_t b = edge.b;
				FloatType w = edge.w;
				int type = structures()->getInt(a);

				graph[i].mergeQuality = mergeQuality(type, curWeights[a], curWeights[b], w);
			}

time[2] = clock();

			// Sort edges by merge quality
			std::sort(graph.begin(), graph.end(), [](GraphEdge e, GraphEdge f) {return e.mergeQuality > f.mergeQuality;});
time[3] = clock();

			// Perform greedy matching of unmatched vertices
			std::fill(hit.begin(), hit.end(), 0);
time[4] = clock();

			for (auto edge: graph) {

				size_t a = edge.a;
				size_t b = edge.b;

				int num_hit = hit[a] + hit[b];
				if (edge.mergeQuality >= _mergingThreshold && num_hit == 0) {
					uf.merge(a, b);
					hit[a] = true;
					hit[b] = true;
					merged = true;

					//if(!task()->incrementProgressValue()) return false;
				}
			}

time[5] = clock();

			// Contract graph
			for (size_t i=0;i<graph.size();i++) {

				auto edge = graph[i];

				size_t a = edge.a;
				size_t b = edge.b;
				FloatType w = edge.w;

				size_t parentA = uf.find(a);
				size_t parentB = uf.find(b);
				if (parentA == parentB) {
					// Edge endpoints are now in same cluster
					// Add the edge weight to the cluster weight
					uf.weights[parentA] += w;

					// Remove the edge from the graph (not included in contracted graph)
					graph[i] = graph.back();
					graph.pop_back();
					i--;
				}
				else {
					// Edge endpoints are still in different clusters
					size_t keymin = parentA, keymax = parentB;
					if (keymin > keymax) {
						std::swap(keymin, keymax);
					}

					graph[i].a = keymin;
					graph[i].b = keymax;
				}
			}

time[6] = clock();

			// Sort graph such that edges connecting the same pair of clusters are adjacent
			std::sort(graph.begin(), graph.end(), [](GraphEdge& e, GraphEdge& f) { return e.a == f.a ? e.b < f.b : e.a < f.a; });

time[7] = clock();

			// Combine edges which connect the same pair of clusters
			std::vector< GraphEdge > temp;
			for (auto edge: graph) {

				bool same = temp.size() > 0 && temp.back().a == edge.a && temp.back().b == edge.b;
				if (same) {
					temp.back().w += edge.w;
				}
				else {
					temp.push_back(edge);
				}
			}
			graph.clear();
			for (auto edge: temp) {
				graph.push_back(edge);
			}

time[8] = clock();

for (int i=0;i<8;i++)
	total_time[i] += time[i + 1] - time[i];

#if 0
printf("!\t");
for(size_t i = 0; i < numAtoms; i++) {
	printf("%lu ", uf.find(i));
}
printf("\n");
#endif

		}
	});

#if 0
#endif
for (int i=0;i<8;i++) {
	long int clockTicksTaken = total_time[i];
	double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
	printf("\ttime taken: %f\n", timeInSeconds);
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
