////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
//  Copyright 2020 Peter Mahler Larsen
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
#include <ovito/stdobj/table/DataTable.h>
#include <ovito/particles/util/NearestNeighborFinder.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include <ovito/particles/util/PTMNeighborFinder.h>
#include "GrainSegmentationEngine.h"
#include "GrainSegmentationModifier.h"
#include "DisjointSet.h"
#include "ThresholdSelection.h"

#include <boost/heap/priority_queue.hpp>
#include <ptm/ptm_functions.h>

#define DEBUG_OUTPUT 0
#if DEBUG_OUTPUT
#include <sys/time.h>
#endif

namespace Ovito { namespace CrystalAnalysis {

#ifndef Q_CC_MSVC
constexpr int GrainSegmentationEngine1::MAX_DISORDERED_NEIGHBORS; // Definition is required by C++14 but not C++17 or later.
#endif

/******************************************************************************
* Constructor.
******************************************************************************/
GrainSegmentationEngine1::GrainSegmentationEngine1(
			ParticleOrderingFingerprint fingerprint, 
			ConstPropertyPtr positions,
			ConstPropertyPtr structureProperty,
			ConstPropertyPtr orientationProperty,
			ConstPropertyPtr correspondenceProperty,
			const SimulationCell& simCell,
			GrainSegmentationModifier::MergeAlgorithm algorithmType, 
			bool handleCoherentInterfaces,
			bool outputBonds) :
	_inputFingerprint(std::move(fingerprint)),
	_positions(std::move(positions)),
	_simCell(simCell),
	_algorithmType(algorithmType),
	_handleBoundaries(handleCoherentInterfaces),
	_structureTypes(structureProperty),
	_orientations(orientationProperty),
	_correspondences(correspondenceProperty),
	_outputBondsToPipeline(outputBonds)
{
	_numParticles = _positions->size();
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine1::perform()
{
	// First phase of grain segmentation algorithm:
	if(!identifyAtomicStructures()) return;
	if(!rotateHexagonalAtoms()) return;
	if(!computeDisorientationAngles()) return;
	if(!determineMergeSequence()) return;

	// Release data that is no longer needed.
	_positions.reset();

	//if(!_outputBondsToPipeline)
	//	decltype(_neighborBonds){}.swap(_neighborBonds);
}


/******************************************************************************
* Performs the PTM algorithm. Determines the local structure type and the
* local lattice orientation at each atomic site.
******************************************************************************/
bool GrainSegmentationEngine1::identifyAtomicStructures()
{
	NearestNeighborFinder neighFinder(PTMAlgorithm::MAX_INPUT_NEIGHBORS);
	if(!neighFinder.prepare(*positions(), cell(), nullptr, this))
		return false;

	setProgressValue(0);
	setProgressMaximum(_numParticles);
	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building neighbor lists"));

	// Mutex is needed to synchronize access to bonds list in parallelized loop.
	std::mutex bondsMutex;

	// Perform analysis on each particle.
	parallelForChunks(_numParticles, *this, [&](size_t startIndex, size_t count, Task& task) {

		// Construct local neighbor list builder.
		NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> neighQuery(neighFinder);
        auto ptmNeighQuery = PTMNeighborFinder();
		ptmNeighQuery.prepare(&neighQuery, correspondences(), structureTypes(), this);

		// Thread-local list of generated bonds connecting neighboring lattice atoms.
		std::vector<NeighborBond> threadlocalNeighborBonds;

		// Loop over a range of input particles.
		for(size_t index = startIndex, endIndex = startIndex + count; index < endIndex; index++) {

			// Update progress indicator (only occasionally).
			if((index % 256) == 0)
				task.incrementProgressValue(256);

			// Break out of loop when computation was canceled.
			if(task.isCanceled())
				break;

			// Get PTM information
			ptmNeighQuery.findNeighbors(index);
			auto structureType = ptmNeighQuery.structureType;
			int numNeighbors = ptmNeighQuery.results().size();
			if (structureType == PTMAlgorithm::OTHER) {
				numNeighbors = std::min(numNeighbors, (int)MAX_DISORDERED_NEIGHBORS);
			}

			for(int j = 0; j < numNeighbors; j++) {
				size_t neighborIndex = ptmNeighQuery.results()[j].index;
				FloatType length = sqrt(ptmNeighQuery.results()[j].distanceSq);

// TODO: apply canonical selection here rather than just using particle indices
				// Create a bond to the neighbor, but skip every other bond to create just one bond per particle pair.
				if(index < neighborIndex)
					threadlocalNeighborBonds.push_back({index,
                                                        neighborIndex,
                                                        std::numeric_limits<FloatType>::infinity(),
                                                        length});

				// Check if neighbor vector spans more than half of a periodic simulation cell.
				Vector3 neighborVector = ptmNeighQuery.results()[j].delta;
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

		// Append thread-local bonds to global bonds list.
		std::lock_guard<std::mutex> lock(bondsMutex);
		_neighborBonds.insert(_neighborBonds.end(), threadlocalNeighborBonds.cbegin(), threadlocalNeighborBonds.cend());
	});
	if(isCanceled())
		return false;

	return !isCanceled();
}

bool GrainSegmentationEngine1::interface_cubic_hex(NeighborBond& bond, bool parent_fcc, bool parent_dcub,
												   FloatType& disorientation, Quaternion& output, size_t& index)
{
	disorientation = std::numeric_limits<FloatType>::infinity();
	index = std::numeric_limits<size_t>::max();
	// Check for a coherent interface (or a crystalline bond, which we check below)
	if (!isCrystallineBond(bond)) {
		return false;
	}

	auto a = bond.a;
	auto b = bond.b;
	if (_adjustedStructureTypes[a] == _adjustedStructureTypes[b]) {
		return false;
	}

	// We want ordering of (a, b) to be (parent phase, defect phase)
	bool flipped = false;
	flipped |= parent_fcc   && _adjustedStructureTypes[a] == PTMAlgorithm::HCP;
	flipped |= !parent_fcc  && _adjustedStructureTypes[a] == PTMAlgorithm::FCC;
	flipped |= parent_dcub  && _adjustedStructureTypes[a] == PTMAlgorithm::HEX_DIAMOND;
	flipped |= !parent_dcub && _adjustedStructureTypes[a] == PTMAlgorithm::CUBIC_DIAMOND;
	if (flipped) {
		std::swap(a, b);
	}

	disorientation = PTMAlgorithm::calculate_interfacial_disorientation(_adjustedStructureTypes[a],
																		_adjustedStructureTypes[b],
																		_adjustedOrientations[a],
																		_adjustedOrientations[b],
																		output);
	index = b;
	return true;
}

/******************************************************************************
* Rotates hexagonal atoms (HCP and hex-diamond) to an equivalent cubic orientation.
******************************************************************************/
bool GrainSegmentationEngine1::rotateHexagonalAtoms()
{
	ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray(structureTypes());
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());

	// Make a copy of structure types and orientations.
	_adjustedStructureTypes = std::vector<PTMAlgorithm::StructureType>(structuresArray.cbegin(), structuresArray.cend());
	_adjustedOrientations = std::vector<Quaternion>(orientationsArray.cbegin(), orientationsArray.cend());

	// Only rotate hexagonal atoms if handling of coherent interfaces is enabled
	if (!_handleBoundaries)
		return true;

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - rotating minority atoms"));

	// Count structure types
	int structureCounts[PTMAlgorithm::NUM_STRUCTURE_TYPES] = {0};
	for (auto structureType: structuresArray) {
		structureCounts[(int)structureType]++;
	}

	bool parent_fcc = structureCounts[(size_t)PTMAlgorithm::FCC] >= structureCounts[(size_t)PTMAlgorithm::HCP];
	bool parent_dcub = structureCounts[(size_t)PTMAlgorithm::CUBIC_DIAMOND] >= structureCounts[(size_t)PTMAlgorithm::HEX_DIAMOND];

	// Set structure targets (i.e. which way a structure will flip)
	PTMAlgorithm::StructureType target[PTMAlgorithm::NUM_STRUCTURE_TYPES];
	if (parent_fcc) {
		target[(size_t)PTMAlgorithm::HCP] = PTMAlgorithm::FCC;
	}
	else {
		target[(size_t)PTMAlgorithm::FCC] = PTMAlgorithm::HCP;
	}

	if (parent_dcub) {
		target[(size_t)PTMAlgorithm::HEX_DIAMOND] = PTMAlgorithm::CUBIC_DIAMOND;
	}
	else {
		target[(size_t)PTMAlgorithm::CUBIC_DIAMOND] = PTMAlgorithm::HEX_DIAMOND;
	}


	NearestNeighborFinder neighFinder(PTMAlgorithm::MAX_INPUT_NEIGHBORS);
	if(!neighFinder.prepare(*positions(), cell(), nullptr, this))
		return false;

	// Construct local neighbor list builder.
	NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> neighQuery(neighFinder);
    auto ptmNeighQuery = PTMNeighborFinder();
	if(!ptmNeighQuery.prepare(&neighQuery, correspondences(), structureTypes(), this))
		return false;


	// TODO: replace comparator with a lambda function
	boost::heap::priority_queue<NeighborBond, boost::heap::compare<PriorityQueueCompare>> pq;

	size_t index;
	Quaternion rotated;
	FloatType disorientation;

	// Populate priority queue with bonds at an cubic-hexagonal interface
	for (auto bond : _neighborBonds) {
		if (interface_cubic_hex(bond, parent_fcc, parent_dcub, disorientation, rotated, index)
			&& disorientation < _misorientationThreshold) {
			pq.push({bond.a, bond.b, disorientation});
		}
	}

	while (pq.size()) {
		auto bond = *pq.begin();
		pq.pop();

		if (!interface_cubic_hex(bond, parent_fcc, parent_dcub, disorientation, rotated, index)) {
			continue;
		}

		// flip structure from 'defect' phase to parent phase and adjust orientation
		auto defectStructureType = _adjustedStructureTypes[index];
		auto targetStructureType = target[(size_t)defectStructureType];
		_adjustedStructureTypes[index] = targetStructureType;
		_adjustedOrientations[index] = rotated;

		// find neighbors to add to the queue
		ptmNeighQuery.findNeighbors(index);
		int numNeighbors = ptmNeighQuery.results().size();
		for(int j = 0; j < numNeighbors; j++) {
			size_t neighborIndex = ptmNeighQuery.results()[j].index;

			size_t dummy;
			if (interface_cubic_hex(bond, parent_fcc, parent_dcub, disorientation, rotated, dummy)
				&& disorientation < _misorientationThreshold) {
				pq.push({index, neighborIndex, disorientation});
			}
		}
	}

	return !isCanceled();
}

/******************************************************************************
* Calculates the disorientation angle for each graph edge (i.e. bond).
******************************************************************************/
bool GrainSegmentationEngine1::computeDisorientationAngles()
{
	// Compute disorientation angles associated with the neighbor graph edges.
	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - misorientation calculation"));

	parallelFor(_neighborBonds.size(), *this, [&](size_t bondIndex) {
		NeighborBond& bond = _neighborBonds[bondIndex];
		bond.disorientation = PTMAlgorithm::calculate_disorientation(_adjustedStructureTypes[bond.a],
																	 _adjustedStructureTypes[bond.b],
																	 _adjustedOrientations[bond.a],
																	 _adjustedOrientations[bond.b]);
	});
	if(isCanceled()) return false;

	// Sort graph edges by disorientation.
	boost::sort(_neighborBonds, [](const NeighborBond& a, const NeighborBond& b) {
		return a.disorientation < b.disorientation;
	});

	return !isCanceled();
}

/******************************************************************************
* Computes the disorientation angle between two crystal clusters of the 
* given lattice type. Furthermore, the function computes the weighted average
* of the two cluster orientations. The norm of the two input quaternions 
* and the output quaternion represents the size of the clusters.
******************************************************************************/
FloatType GrainSegmentationEngine1::calculate_disorientation(int structureType, Quaternion& qa, const Quaternion& qb)
{
	FloatType qa_norm = qa.norm();
	FloatType qb_norm = qb.norm();
	double qtarget[4] = { qa.w()/qa_norm, qa.x()/qa_norm, qa.y()/qa_norm, qa.z()/qa_norm };
	double q[4]	   = { qb.w()/qb_norm, qb.x()/qb_norm, qb.y()/qb_norm, qb.z()/qb_norm };

	// Convert structure type back to PTM representation
	int type = 0;
	if(structureType == PTMAlgorithm::OTHER) {
		qWarning() << "Grain segmentation: remap failure - disordered structure input";
		return std::numeric_limits<FloatType>::max();
	}
	else if(structureType == PTMAlgorithm::FCC) type = PTM_MATCH_FCC;
	else if(structureType == PTMAlgorithm::HCP) type = PTM_MATCH_HCP;
	else if(structureType == PTMAlgorithm::BCC) type = PTM_MATCH_BCC;
	else if(structureType == PTMAlgorithm::SC) type = PTM_MATCH_SC;
	else if(structureType == PTMAlgorithm::CUBIC_DIAMOND) type = PTM_MATCH_DCUB;
	else if(structureType == PTMAlgorithm::HEX_DIAMOND) type = PTM_MATCH_DHEX;
	else if(structureType == PTMAlgorithm::GRAPHENE) type = PTM_MATCH_GRAPHENE;

	double disorientation = 0;
	int8_t dummy_mapping[PTM_MAX_POINTS];
	if(ptm_remap_template(type, true, 0, qtarget, q, &disorientation, dummy_mapping, nullptr) < 0) {
		qWarning() << "Grain segmentation: remap failure";
		OVITO_ASSERT(false); // remap failure
	}

	qa.w() += q[0] * qb_norm;
	qa.x() += q[1] * qb_norm;
	qa.y() += q[2] * qb_norm;
	qa.z() += q[3] * qb_norm;
	return disorientation;
}

/******************************************************************************
* Clustering using minimum spanning tree algorithm.
******************************************************************************/
bool GrainSegmentationEngine1::minimum_spanning_tree_clustering(
		std::vector<Quaternion>& qsum, DisjointSet& uf)
{
	size_t progress = 0;
	for(const NeighborBond& edge : _neighborBonds) {

		if (edge.disorientation < _misorientationThreshold) {
			size_t pa = uf.find(edge.a);
			size_t pb = uf.find(edge.b);
			if(pa != pb && isCrystallineBond(edge)) {
				size_t parent = uf.merge(pa, pb);
				size_t child = (parent == pa) ? pb : pa;
				FloatType disorientation = calculate_disorientation(_adjustedStructureTypes[parent], qsum[parent], qsum[child]);
				OVITO_ASSERT(edge.a < edge.b);
				_dendrogram.emplace_back(parent, child, edge.disorientation, disorientation, 1, qsum[parent]);
			}
		}

		// Update progress indicator.
		if((progress++ % 1024) == 0) {
			if(!incrementProgressValue(1024)) 
				return false;
		}
	}

	return !isCanceled();
}

/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine1::determineMergeSequence()
{
	// Build graph.
	if(_algorithmType == GrainSegmentationModifier::GraphClusteringAutomatic || _algorithmType == GrainSegmentationModifier::GraphClusteringManual) {

		setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
		setProgressValue(0);
		setProgressMaximum(neighborBonds().size());

		size_t progress = 0;
		for (auto edge: neighborBonds()) {
			if (isCrystallineBond(edge) && edge.disorientation < _misorientationThreshold) {
				FloatType weight = calculateGraphWeight(edge.disorientation);
				graph.add_edge(edge.a, edge.b, weight);
			}

			if((progress++ % 1024) == 0) {
				if(!incrementProgressValue(1024)) 
					return false;
			}
		}
	}

	// Build dendrogram.
	std::vector<Quaternion> qsum(_adjustedOrientations.cbegin(), _adjustedOrientations.cend());
	DisjointSet uf(_numParticles);
	_dendrogram.resize(0);

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	setProgressValue(0);
	setProgressMaximum(_numParticles);  //TODO: make this num. crystalline particles

	if(_algorithmType == GrainSegmentationModifier::GraphClusteringAutomatic || _algorithmType == GrainSegmentationModifier::GraphClusteringManual) {
		node_pair_sampling_clustering(graph, qsum);
	}
	else {
		minimum_spanning_tree_clustering(qsum, uf);
	}
	if(isCanceled())
		return false;

	// Sort dendrogram entries by distance.
	boost::sort(_dendrogram, [](const DendrogramNode& a, const DendrogramNode& b) { return a.distance < b.distance; });

	if(isCanceled())
		return false;

#if DEBUG_OUTPUT
char filename[128];
struct timeval tp;
gettimeofday(&tp, NULL);
long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
sprintf(filename, "dump_%lu.txt", ms);
FILE* fout = fopen(filename, "w");
#endif

	// Scan through the entire merge list to determine merge sizes.
	size_t numPlot = 0;
	uf.clear();
	for(DendrogramNode& node : _dendrogram) {
		size_t sa = uf.nodesize(uf.find(node.a));
		size_t sb = uf.nodesize(uf.find(node.b));
		size_t dsize = std::min(sa, sb);
		node.gm_size = sqrt(sa * sb);
		uf.merge(node.a, node.b);

#if DEBUG_OUTPUT
if (fout)
	fprintf(fout, "%lu %lu %lu %lu %lu %e\n", node.a, node.b, sa, sb, dsize, node.distance);
#endif

		// We don't want to plot very small merges - they extend the x-axis by a lot and don't provide much useful information
		node.size = dsize;
		if(dsize >= _minPlotSize) {
			numPlot++;
		}
	}

#if DEBUG_OUTPUT
fclose(fout);
#endif

	if(_algorithmType == GrainSegmentationModifier::GraphClusteringAutomatic || _algorithmType == GrainSegmentationModifier::GraphClusteringManual) {

		// Create PropertyStorage objects for the output plot.
		PropertyAccess<FloatType> mergeDistanceArray = _mergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), false, DataTable::XProperty);
		PropertyAccess<FloatType> mergeSizeArray = _mergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Delta merge size"), false, DataTable::YProperty);

		// Generate output data plot points from dendrogram data.
		FloatType* mergeDistanceIter = mergeDistanceArray.begin();
		FloatType* mergeSizeIter = mergeSizeArray.begin();
		for(const DendrogramNode& node : _dendrogram) {
			if(node.size >= _minPlotSize) {
				*mergeDistanceIter++ = std::log(node.distance);
				*mergeSizeIter++ = node.size;
			}
		}

		auto regressor = ThresholdSelection::Regressor(_dendrogram);
		_suggestedMergingThreshold = regressor.calculate_threshold(_dendrogram, 1.5);

		// Create PropertyStorage objects for the output plot.
		numPlot = 0;
		for (auto y: regressor.ys) {
			numPlot += (y > 0) ? 1 : 0; // plot positive distances only, for clarity
		}

		PropertyAccess<FloatType> logMergeSizeArray = _logMergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log geometric merge size"), false, DataTable::XProperty);
		PropertyAccess<FloatType> logMergeDistanceArray = _logMergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), false, DataTable::YProperty);

		// Generate output data plot points from dendrogram data.
		FloatType* logMergeDistanceIter = logMergeDistanceArray.begin();
		FloatType* logMergeSizeIter = logMergeSizeArray.begin();
		for (size_t i=0;i<regressor.residuals.size();i++) {
			if (regressor.ys[i] > 0) {
				*logMergeSizeIter++ = regressor.xs[i];
				*logMergeDistanceIter++ = regressor.ys[i];
			}
		}

	}
	else {
		// Create PropertyStorage objects for the output plot.
		PropertyAccess<FloatType> mergeDistanceArray = _mergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Misorientation (degrees)"), false, DataTable::XProperty);
		PropertyAccess<FloatType> mergeSizeArray = _mergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Merge size"), false, DataTable::YProperty);

		// Generate output data plot points from dendrogram data.
		FloatType* mergeDistanceIter = mergeDistanceArray.begin();
		FloatType* mergeSizeIter = mergeSizeArray.begin();
		for(const DendrogramNode& node : _dendrogram) {
			if(node.size >= _minPlotSize) {
				*mergeDistanceIter++ = node.distance;
				*mergeSizeIter++ = node.size;
			}
		}
	}

	return !isCanceled();
}

/******************************************************************************
* Creates another engine that performs the next stage of the computation. 
******************************************************************************/
std::shared_ptr<AsynchronousModifier::Engine> GrainSegmentationEngine1::createContinuationEngine(ModifierApplication* modApp, const PipelineFlowState& input)
{
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(modApp->modifier());

	return std::make_shared<GrainSegmentationEngine2>(
		static_pointer_cast<GrainSegmentationEngine1>(shared_from_this()),
		modifier->mergingThreshold(),
		modifier->orphanAdoption(),
		modifier->minGrainAtomCount()
	);
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine2::perform()
{
	// Second phase: Execute merge steps up to the threshold set by the user or the adaptively determined threshold.
	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - merging clusters"));

	// Either use user-defined merge threshold or automatically computed threshold.
	FloatType mergingThreshold = _mergingThreshold;
	if(_engine1->_algorithmType == GrainSegmentationModifier::GraphClusteringAutomatic) {
		mergingThreshold = _engine1->suggestedMergingThreshold();
	}

	if(_engine1->_algorithmType == GrainSegmentationModifier::MinimumSpanningTree) {
		mergingThreshold = log(mergingThreshold);
	}

	const std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram = _engine1->_dendrogram;

#if 0
	// Refine the graph partitions
	if(_engine1->_algorithmType == GrainSegmentationModifier::GraphClusteringAutomatic) {
		auto graph = _engine1->graph;
		FloatType gamma = 1 / mergingThreshold;	 // resolution parameter

		for(auto node = dendrogram.crbegin(); node != dendrogram.crend(); ++node) {
			graph.reinstate_edge(node->a, node->b);
			//printf("%e %e\n", node->distance, graph.wnode[node->a] * graph.wnode[node->b] / graph.adj[node->a][node->b]);

			if (std::log(node->distance) < mergingThreshold) {
				// test all options for reassignment
			}
		}
	}
#endif

	ConstPropertyAccess<Quaternion> orientationsArray(_engine1->orientations());
	std::vector<Quaternion> meanOrientation(orientationsArray.cbegin(), orientationsArray.cend());

	// Iterate through merge list until distance cutoff is met.
	DisjointSet uf(_numParticles);
	for(auto node = dendrogram.cbegin(); node != dendrogram.cend(); ++node) {
		if(isCanceled()) 
			return;

		if(std::log(node->distance) > mergingThreshold)
			break;

		uf.merge(node->a, node->b);
		size_t parent = uf.find(node->a);
		OVITO_ASSERT(node->orientation.norm() > FLOATTYPE_EPSILON);
		meanOrientation[parent] = node->orientation;
	}

	// Relabels the clusters to obtain a contiguous sequence of cluster IDs.
	std::vector<size_t> clusterRemapping(_numParticles);

	// Assign new consecutive IDs to root clusters.
	_numClusters = 1;
	ConstPropertyAccess<int> structuresArray(_engine1->structureTypes());
	std::vector<int> clusterStructureTypes;
	std::vector<Quaternion> clusterOrientations;
	for(size_t i = 0; i < _numParticles; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(uf.nodesize(i) < _minGrainAtomCount || structuresArray[i] == PTMAlgorithm::OTHER) {
				clusterRemapping[i] = 0;
			}
			else {
				clusterRemapping[i] = _numClusters;
				_numClusters++;
				clusterStructureTypes.push_back(structuresArray[i]);
				clusterOrientations.push_back(meanOrientation[i].normalized());
			}
		}
	}
	if(isCanceled()) 
		return;

	// Allocate and fill output array storing the grain IDs (1-based identifiers). 
	_grainIds =  std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int64, 1, 0, QStringLiteral("Grain Identifier"), false, DataTable::XProperty);
	boost::algorithm::iota_n(PropertyAccess<qlonglong>(_grainIds).begin(), size_t(1), _grainIds->size());
	if(isCanceled()) 
		return;

	// Allocate output array storing the grain sizes.
	_grainSizes = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int64, 1, 0, QStringLiteral("Grain Size"), true, DataTable::YProperty);

	// Allocate output array storing the structure type of grains.
	_grainStructureTypes = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int, 1, 0, QStringLiteral("Structure Type"), false);
	boost::copy(clusterStructureTypes, PropertyAccess<int>(_grainStructureTypes).begin());
	if(isCanceled()) 
		return;

	// Allocate output array with each grain's unique color.
	// Fill it with random color values (using constant random seed to keep it reproducible).
	_grainColors = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 3, 0, QStringLiteral("Color"), false, 0, QStringList() << QStringLiteral("R") << QStringLiteral("G") << QStringLiteral("B"));
	std::default_random_engine rng(1);
	std::uniform_real_distribution<FloatType> uniform_dist(0, 1);
	boost::generate(PropertyAccess<Color>(_grainColors), [&]() { return Color::fromHSV(uniform_dist(rng), 1.0 - uniform_dist(rng) * 0.8, 1.0 - uniform_dist(rng) * 0.5); });
	if(isCanceled()) 
		return;

	// Allocate output array storing the mean lattice orientation of grains (represented by a quaternion).
	_grainOrientations = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 4, 0, QStringLiteral("Orientation"), true, 0, QStringList() << QStringLiteral("X") << QStringLiteral("Y") << QStringLiteral("Z") << QStringLiteral("W"));
	boost::copy(clusterOrientations, PropertyAccess<Quaternion>(_grainOrientations).begin());

	// Determine new IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		clusterRemapping[particleIndex] = clusterRemapping[uf.find(particleIndex)];
	if(isCanceled()) 
		return;

	// Relabel atoms after cluster IDs have changed.
	// Also count the number of atoms in each cluster.
	PropertyAccess<qlonglong> atomClustersArray(atomClusters());
	PropertyAccess<qlonglong> grainSizeArray(_grainSizes);
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++) {
		size_t gid = clusterRemapping[particleIndex];
		atomClustersArray[particleIndex] = gid;
		if(gid != 0) grainSizeArray[gid - 1]++;
	}
	if(isCanceled()) 
		return;

	// Reorder grains by size (large to small).
	if(_numClusters > 1) {

		// Determine the index remapping for reordering the grain list by size.
		std::vector<size_t> mapping(_numClusters - 1);
		std::iota(mapping.begin(), mapping.end(), size_t(0));
		std::sort(mapping.begin(), mapping.end(), [&](size_t a, size_t b) {
			return grainSizeArray[a] > grainSizeArray[b];
		});
		if(isCanceled()) 
			return;

		// Use index map to reorder grain data arrays.

		PropertyPtr originalGrainSizes = _grainSizes;
		PropertyStorage::makeMutable(_grainSizes);
		originalGrainSizes->mappedCopyTo(*_grainSizes, mapping);

		PropertyPtr originalGrainStructureTypes = _grainStructureTypes;
		PropertyStorage::makeMutable(_grainStructureTypes);
		originalGrainStructureTypes->mappedCopyTo(*_grainStructureTypes, mapping);

		PropertyPtr originalGrainOrientations = _grainOrientations;
		PropertyStorage::makeMutable(_grainOrientations);
		originalGrainOrientations->mappedCopyTo(*_grainOrientations, mapping);
		if(isCanceled()) 
			return;

		// Invert the grain index map. 

		std::vector<size_t> inverseMapping(_numClusters);
		inverseMapping[0] = 0; // Keep cluster ID 0 in place.
		for(size_t i = 1; i < _numClusters; i++)
			inverseMapping[mapping[i-1]+1] = i;

		// Remap per-particle grain IDs.

		for(auto& id : atomClustersArray)
			id = inverseMapping[id];
		if(isCanceled()) 
			return;

		// Adopt orphan atoms.
		if(_adoptOrphanAtoms)
			mergeOrphanAtoms();
	}
}

/******************************************************************************
* Merges any orphan atoms into the closest cluster.
******************************************************************************/
bool GrainSegmentationEngine2::mergeOrphanAtoms()
{
	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - merging orphan atoms"));
	setProgressValue(0);

	PropertyAccess<qlonglong> atomClustersArray(atomClusters());
	PropertyAccess<qlonglong> grainSizeArray(_grainSizes);

	/// The bonds connecting neighboring non-crystalline atoms.
	std::vector<GrainSegmentationEngine1::NeighborBond> noncrystallineBonds;
	for (auto nb: _engine1->neighborBonds()) {
		if (atomClustersArray[nb.a] == 0 || atomClustersArray[nb.b] == 0) {
			// Add bonds for both atoms
			noncrystallineBonds.push_back(nb);

            std::swap(nb.a, nb.b);
			noncrystallineBonds.push_back(nb);
		}
	}
	if(isCanceled())
		return false;

    boost::sort(noncrystallineBonds,
                 [](const GrainSegmentationEngine1::NeighborBond& a, const GrainSegmentationEngine1::NeighborBond& b)
                 {return a.a < b.a;});

	boost::heap::priority_queue<PQNode, boost::heap::compare<PQCompareLength>> pq;

	// Populate priority queue with bonds at a crystalline-noncrystalline interface
	for (auto bond: _engine1->neighborBonds()) {
        auto clusterA = atomClustersArray[bond.a];
        auto clusterB = atomClustersArray[bond.b];

		if (clusterA != 0 && clusterB == 0) {
			pq.push({clusterA, bond.b, bond.length});
		}
		else if (clusterA == 0 && clusterB != 0) {
			pq.push({clusterB, bond.a, bond.length});
		}
	}

	while (pq.size()) {
		auto node = *pq.begin();
		pq.pop();

        if (atomClustersArray[node.particleIndex] != 0)
            continue;

        atomClustersArray[node.particleIndex] = node.cluster;
        grainSizeArray[node.cluster - 1]++;

		// Get the range of bonds adjacent to the current atom.
		auto bondsRange = boost::range::equal_range(noncrystallineBonds, GrainSegmentationEngine1::NeighborBond{node.particleIndex, 0, 0, 0},
			[](const GrainSegmentationEngine1::NeighborBond& a, const GrainSegmentationEngine1::NeighborBond& b)
            { return a.a < b.a; });

		// Find the closest cluster atom in the neighborhood (using PTM ordering).
		for(const GrainSegmentationEngine1::NeighborBond& bond : boost::make_iterator_range(bondsRange.first, bondsRange.second)) {
			OVITO_ASSERT(bond.a == node.particleIndex);

			auto neighborIndex = bond.b;
			if(neighborIndex == std::numeric_limits<size_t>::max()) break;
			if(atomClustersArray[neighborIndex] != 0) continue;

			pq.push({node.cluster, neighborIndex, node.length + bond.length});
		}
	}

	return !isCanceled();
}

}	// End of namespace
}	// End of namespace
