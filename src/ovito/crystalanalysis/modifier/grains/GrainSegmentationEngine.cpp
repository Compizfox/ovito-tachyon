////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
//  Copyright 2019 Peter Mahler Larsen
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify it either under the
//  terms of the GNU General Public License version 3 as published by the Free Software
//  Foundation (the "GPL") or, at your option, under the terms of the MIT License.
//  If you do not alter this notice, a recipient may use your version of this
//  file under either the GPL or the MIT License.
//
//  You should have received a copy of the GPL along with this program in a
//  file LICENSE.GPL.txt.  You should have received a copy of the MIT License along
//  with this program in a file LICENSE.MIT.txt
//
//  This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND,
//  either express or implied. See the GPL or the MIT License for the specific language
//  governing rights and limitations.
//
////////////////////////////////////////////////////////////////////////////////////////

#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/stdobj/series/DataSeriesObject.h>
#include <ovito/particles/util/NearestNeighborFinder.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include "GrainSegmentationEngine.h"
#include "GrainSegmentationModifier.h"
#include "NodePairSampling.h"
#include "DisjointSet.h"

#include <ptm/ptm_functions.h>
#include <ptm/ptm_quat.h>

namespace Ovito { namespace CrystalAnalysis {

/******************************************************************************
* Constructor.
******************************************************************************/
GrainSegmentationEngine::GrainSegmentationEngine(
			ParticleOrderingFingerprint fingerprint, ConstPropertyPtr positions, const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, ConstPropertyPtr selection,
			FloatType rmsdCutoff, bool algorithmType,
			int minGrainAtomCount, bool orphanAdoption, bool outputBonds) :
	StructureIdentificationModifier::StructureIdentificationEngine(std::move(fingerprint), positions, simCell, std::move(typesToIdentify), std::move(selection)),
	_rmsdCutoff(rmsdCutoff),
	_algorithmType(algorithmType),
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
	if(!determineMergeSequence()) return;

#if 0
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
#endif
}

/******************************************************************************
* Performs the PTM algorithm. Determines the local structure type and the
* local lattice orientation at each atomic site.
******************************************************************************/
bool GrainSegmentationEngine::identifyAtomicStructures()
{
	// Initialize the PTMAlgorithm object.
	PTMAlgorithm algorithm;
	algorithm.setRmsdCutoff(0.0); // Note: We'll do our own RMSD threshold filtering below.

	// Specify the structure types the PTM should look for.
	for(int i = 0; i < typesToIdentify().size() && i < PTMAlgorithm::NUM_STRUCTURE_TYPES; i++) {
		algorithm.setStructureTypeIdentification(static_cast<PTMAlgorithm::StructureType>(i), typesToIdentify()[i]);
	}

	if(!algorithm.prepare(*positions(), cell(), selection().get(), task().get()))
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
	parallelForChunks(positions()->size(), *task(), [this, &cachedNeighbors, &algorithm](size_t startIndex, size_t count, Task& task) {
		// Create a thread-local kernel for the PTM algorithm.
		PTMAlgorithm::Kernel kernel(algorithm);

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
	parallelForChunks(positions()->size(), *task(), [this, &cachedNeighbors, &algorithm, &neighFinder](size_t startIndex, size_t count, Task& task) {
		// Create a thread-local kernel for the PTM algorithm.
		PTMAlgorithm::Kernel kernel(algorithm);

		// Create a neighbor query for disordered atoms
		NearestNeighborFinder::Query<MAX_DISORDERED_NEIGHBORS> neighQuery(neighFinder);

		// Loop over a range of input particles.
		for(size_t index = startIndex, endIndex = startIndex + count; index < endIndex; index++) {

			// Update progress indicator (only occasionally).
			if((index % 256) == 0)
				task.incrementProgressValue(256);

			// Break out of loop when computation was canceled.
			if(task.isCanceled())
				break;

			// Skip input particles that are excluded from the analysis.
			if(selection() && !selection()->getInt(index)) {
				structures()->setInt(index, PTMAlgorithm::OTHER);
				rmsd()->setFloat(index, 0);
				continue;
			}

			// Perform the PTM analysis for the current particle.
			PTMAlgorithm::StructureType type = kernel.identifyStructure(index, cachedNeighbors, nullptr);

			// Store results in the output property arrays.
			structures()->setInt(index, type);
			rmsd()->setFloat(index, kernel.rmsd());

//fprintf(out, "%f %f %f %f %f %f %f %f\n", positions()->getPoint3(index).x(), positions()->getPoint3(index).y(), positions()->getPoint3(index).z(),
//					kernel.rmsd(), kernel.orientation().w(), kernel.orientation().x(), kernel.orientation().y(), kernel.orientation().z());
//fflush(out);

			if(type != PTMAlgorithm::OTHER) {
				// Store computed local lattice orientation in the output property array.
				if(orientations())
					orientations()->setQuaternion(index, kernel.orientation());

				// Store neighbor list for later use.
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

FloatType GrainSegmentationEngine::calculate_disorientation(int structureType, std::vector< Quaternion >& qsum, size_t a, size_t b)
{
	Quaternion na = qsum[a].normalized();
	double qtarget[4] = {na.w(), na.x(), na.y(), na.z()};

	Quaternion nb = qsum[b].normalized();
	double q[4] = {nb.w(), nb.x(), nb.y(), nb.z()};

	// Convert structure type back to PTM representation
	int type = 0;
	if (structureType == PTMAlgorithm::FCC)			type = PTM_MATCH_FCC;
	else if (structureType == PTMAlgorithm::HCP)		type = PTM_MATCH_HCP;
	else if (structureType == PTMAlgorithm::BCC)		type = PTM_MATCH_BCC;
	else if (structureType == PTMAlgorithm::SC)		type = PTM_MATCH_SC;
	else if (structureType == PTMAlgorithm::CUBIC_DIAMOND)	type = PTM_MATCH_DCUB;
	else if (structureType == PTMAlgorithm::HEX_DIAMOND)	type = PTM_MATCH_DHEX;
	else if (structureType == PTMAlgorithm::GRAPHENE)	type = PTM_MATCH_GRAPHENE;

	double disorientation = 0;
	int8_t dummy_mapping[PTM_MAX_POINTS];
	if (ptm_remap_template(type, true, 0, qtarget, q, &disorientation, dummy_mapping, nullptr) < 0) {
		qWarning() << "remap failure";
		OVITO_ASSERT(false); // remap failure
	}

	FloatType norm = qsum[b].norm();
	qsum[a].w() += q[0] * norm;
	qsum[a].x() += q[1] * norm;
	qsum[a].y() += q[2] * norm;
	qsum[a].z() += q[3] * norm;
	return disorientation;
}

bool GrainSegmentationEngine::minimum_spanning_tree_clustering(std::vector<GraphEdge>& initial_graph, DisjointSet& uf, size_t start, size_t end,
								DendrogramNode* dendrogram, int structureType, std::vector<Quaternion>& qsum)
{
	Graph graph;
	for(size_t i = start; i < end; i++) {
		const auto& edge = initial_graph[i];
		auto a = edge.a;
		auto b = edge.b;
		auto w = edge.w;
		if(uf.find(a) != uf.find(b)) {
			size_t parent = uf.merge(a, b);
			FloatType disorientation;
			if(parent == a) {
				disorientation = calculate_disorientation(structureType, qsum, a, b);
			}
			else {
				disorientation = calculate_disorientation(structureType, qsum, b, a);
			}
			*dendrogram++ = DendrogramNode(std::min(a, b), std::max(a, b), w, disorientation, 1);
		}
	}

	return true;
}

/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine::determineMergeSequence()
{
	size_t numAtoms = positions()->size();

	if(_numSuperclusters == 1) {
		_numClusters = 1;
		for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {
			atomClusters()->setInt64(particleIndex, _atomSuperclusters[particleIndex]);
		}
		return true;
	}

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(latticeNeighborBonds().size());

	// Calculate a contiguous atom index mapping
	std::vector< std::tuple< size_t, size_t > > atomIntervals(_numSuperclusters);
	for(size_t i = 0, start = 0; i < _numSuperclusters; i++) {
		std::get<0>(atomIntervals[i]) = start;
		std::get<1>(atomIntervals[i]) = start + _superclusterSizes[i];
		start += _superclusterSizes[i];
	}

	if(task()->isCanceled())
		return false;

	// Build initial graph
	FloatType totalWeight = 0;
	std::vector< GraphEdge > initial_graph;
	std::vector< size_t > bondCount(_numSuperclusters, 0);	// number of bonds in each supercluster
	const FloatType* disorientation = neighborDisorientationAngles()->constDataFloat();

	for(size_t i = 0; i < latticeNeighborBonds().size(); i++) {
		if(!task()->incrementProgressValue()) return false;

		const Bond& bond = latticeNeighborBonds()[i];
		FloatType dis = disorientation[i];
		size_t a = bond.index1;
		size_t b = bond.index2;
		size_t sc = _atomSuperclusters[a];

		// Skip high-angle edges.
		if (sc == 0 || dis > _misorientationThreshold) {
			bondCount[0]++;
		}
		else {
			// Convert disorientations to graph weights
			FloatType deg = qRadiansToDegrees(dis);
			initial_graph.emplace_back(a, b, deg, sc);
			bondCount[sc]++;

			FloatType weight = std::exp(-FloatType(1)/3 * deg * deg);	//this is fairly arbitrary but it works well
			totalWeight += weight;
		}
	}

	size_t c = 0;
	std::vector< size_t > bondStart(_numSuperclusters, 0);
	for (size_t i=1;i<_numSuperclusters;i++) {
		bondStart[i] = c;
		c += bondCount[i];
	}

	std::sort(initial_graph.begin(), initial_graph.end(), [](GraphEdge& e, GraphEdge& f)
	{
		if (e.superCluster == f.superCluster) {
			return e.w < f.w;
		}
		else {
			return e.superCluster < f.superCluster;
		}
	});

	if(task()->isCanceled())
		return false;

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(initial_graph.size());

	c = 0;
	std::vector< size_t > dendrogramOffsets(_numSuperclusters, 0);
	for (size_t sc=1;sc<_numSuperclusters;sc++) {
		dendrogramOffsets[sc] = c;
		if (_superclusterSizes[sc] > 0) {		// this should always be non-zero
			c += _superclusterSizes[sc] - 1;	// number of edges in a tree = number of vertices - 1
		}
	}

	DisjointSet uf(numAtoms);
	_dendrogram.resize(c);

	std::vector<Quaternion> qsum(numAtoms);
	for(size_t particleIndex=0; particleIndex < numAtoms; particleIndex++) {
		const Quaternion& q = orientations()->getQuaternion(particleIndex);
		qsum[particleIndex].w() = q.w();
		qsum[particleIndex].x() = q.x();
		qsum[particleIndex].y() = q.y();
		qsum[particleIndex].z() = q.z();
	}

	parallelFor(_numSuperclusters, *task(), [this, &initial_graph, &uf, &bondStart, &bondCount, &dendrogramOffsets, &totalWeight, &qsum](size_t sc) {
		if(sc == 0) return;
		size_t start = bondStart[sc];
		size_t end = bondStart[sc] + bondCount[sc];
		size_t index = dendrogramOffsets[sc];
		int structureType = structures()->getInt(initial_graph[start].a);

		if (_algorithmType == 0) {
			minimum_spanning_tree_clustering(initial_graph, uf, start, end, &_dendrogram[index], structureType, qsum);
		}
		else {
			node_pair_sampling_clustering(initial_graph, start, end, totalWeight, &_dendrogram[index], structureType, qsum);
		}
	});

	std::sort(_dendrogram.begin(), _dendrogram.end(),
			[](DendrogramNode& a, DendrogramNode& b) { return a.distance < b.distance; });

	// Scan through the entire merge list to determine merge sizes.
	size_t numPlot = 0;
	uf.clear();
	for(DendrogramNode& node : _dendrogram) {
		size_t sa = uf.sizes[uf.find(node.a)];
		size_t sb = uf.sizes[uf.find(node.b)];
		size_t dsize = std::min(sa, sb);
		uf.merge(node.a, node.b);
		node.disorientation = qRadiansToDegrees(node.disorientation);

		// We don't want to plot very small merges - they extend the x-axis by a lot and don't provide much useful information
		if(dsize < _minGrainAtomCount / 2) {
			node.size = 0;
		}
		else {
			node.size = dsize;
			numPlot++;
		}
	}

	// Create PropertyStorage output objects
	_mergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), false, DataSeriesObject::XProperty);
	_mergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Merge size"), false, DataSeriesObject::YProperty);

	c = 0;
	for(const DendrogramNode& node : _dendrogram) {
		if (node.size > 0) {	// keep only those which we want to plot
			_mergeDistance->setFloat(c, log(node.distance));
			_mergeSize->setFloat(c, node.size);
			c++;
		}
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Executes precomputed merge steps up to the threshold value set by the user.
******************************************************************************/
void GrainSegmentationEngine::executeMergeSequence(FloatType mergingThreshold)
{
	size_t numAtoms = _atomSuperclusters.size();

	// Iterate through merge list until distance cutoff is met.
	DisjointSet uf(numAtoms);
	for(const DendrogramNode& node : _dendrogram) {
		if(log(node.distance) >= mergingThreshold)
			break;

		uf.merge(node.a, node.b);
	}

	// Relabels the clusters to obtain a contiguous sequence of cluster IDs.
	std::vector<size_t> clusterRemapping(numAtoms);

	// Assign new consecutive IDs to root clusters.
	_numClusters = 1;
	for(size_t i = 0; i < numAtoms; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(uf.sizes[i] < _minGrainAtomCount) {
				clusterRemapping[i] = 0;
			}
			else {
				clusterRemapping[i] = _numClusters;
				_numClusters++;
			}
		}
	}

	// Determine new IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++)
		clusterRemapping[particleIndex] = clusterRemapping[uf.find(particleIndex)];

	// Relabel atoms after cluster IDs have changed.
	std::vector<qlonglong> clusterSizes(_numClusters, 0);
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {
		size_t gid = clusterRemapping[particleIndex];
		atomClusters()->setInt64(particleIndex, gid);
		clusterSizes[gid]++;
	}

	// Relabel clusters by size (large to small)
	if(_numClusters > 1) {
		std::vector<std::tuple<size_t, size_t>> lut;

		for(size_t i = 0; i < _numClusters; i++) {
			lut.emplace_back(std::make_tuple(i, clusterSizes[i]));
		}

		// Sort by size, leaving zeroth cluster in place
		std::sort(lut.begin() + 1, lut.end(),
			[](const std::tuple<size_t, size_t>& a, const std::tuple<size_t, size_t>& b)
			{ return std::get<1>(b) < std::get<1>(a); });

		std::vector<size_t> indices(_numClusters);
		for(size_t gid = 0; gid < _numClusters; gid++) {
			indices[std::get<0>(lut[gid])] = gid;
		}

		for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {
			size_t gid = atomClusters()->getInt64(particleIndex);
			atomClusters()->setInt64(particleIndex, indices[gid]);
		}

		for(size_t i = 0; i < _numClusters; i++) {
			clusterSizes[i] = std::get<1>(lut[i]);
		}
	}

//TODO: output grain sizes as a property
//printf("grain sizes:\n");
//for (int i=0;i<_numClusters;i++) {
//	printf("\t%d\t%lld\n", i, clusterSizes[i]);
//}
}

#if 0
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
#endif

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
//				clusterSizes[newlyAssignedClusters[i] - 1]++;
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

#if 0
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
#endif

}	// End of namespace
}	// End of namespace
