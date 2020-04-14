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
#include "GrainSegmentationEngine.h"
#include "GrainSegmentationModifier.h"
#include "DisjointSet.h"

#include <ptm/ptm_functions.h>
#include <ptm/ptm_quat.h>

#define DEBUG_OUTPUT 0
#if DEBUG_OUTPUT
#include <sys/time.h>
#endif

namespace Ovito { namespace CrystalAnalysis {

/******************************************************************************
* Constructor.
******************************************************************************/
GrainSegmentationEngine1::GrainSegmentationEngine1(
			ParticleOrderingFingerprint fingerprint, 
			ConstPropertyPtr positions, 
			const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, 
			ConstPropertyPtr selection,
			FloatType rmsdCutoff, 
			GrainSegmentationModifier::MergeAlgorithm algorithmType, 
			bool outputBonds) :
	StructureIdentificationModifier::StructureIdentificationEngine(std::move(fingerprint), positions, simCell, std::move(typesToIdentify), std::move(selection)),
	_numParticles(positions->size()),
	_rmsdCutoff(rmsdCutoff),
	_algorithmType(algorithmType),
	_rmsd(std::make_shared<PropertyStorage>(_numParticles, PropertyStorage::Float, 1, 0, QStringLiteral("RMSD"), false)),
	_orientations(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::OrientationProperty, true)),
	_outputBondsToPipeline(outputBonds)
{
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine1::perform()
{
	// First phase of grain segmentation algorithm:
	if(!identifyAtomicStructures()) return;
	if(!computeDisorientationAngles()) return;
	if(!determineMergeSequence()) return;

	// Release data that is no longer needed.
	releaseWorkingData();

	//if(!_outputBondsToPipeline)
	//	decltype(_neighborBonds){}.swap(_neighborBonds);
}

/******************************************************************************
* Performs the PTM algorithm. Determines the local structure type and the
* local lattice orientation at each atomic site.
******************************************************************************/
bool GrainSegmentationEngine1::identifyAtomicStructures()
{
	// Initialize the PTMAlgorithm object.
	PTMAlgorithm algorithm;
	algorithm.setRmsdCutoff(0.0); // Note: We'll do our own RMSD threshold filtering below.

	// Specify the structure types the PTM should look for.
	for(int i = 0; i < typesToIdentify().size() && i < PTMAlgorithm::NUM_STRUCTURE_TYPES; i++) {
		algorithm.setStructureTypeIdentification(static_cast<PTMAlgorithm::StructureType>(i), typesToIdentify()[i]);
	}

	if(!algorithm.prepare(*positions(), cell(), selection(), this))
		return false;

	setProgressValue(0);
	setProgressMaximum(_numParticles);
	setProgressText(GrainSegmentationModifier::tr("Pre-calculating neighbor ordering"));

	// Pre-order neighbors of each particle.
	std::vector<uint64_t> cachedNeighbors(_numParticles);
	parallelForChunks(_numParticles, *this, [this, &cachedNeighbors, &algorithm](size_t startIndex, size_t count, Task& task) {
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

			// Calculate ordering of neighbors
			kernel.cacheNeighbors(index, &cachedNeighbors[index]);
		}
	});
	if(isCanceled() || positions()->size() == 0)
		return false;

	setProgressValue(0);
	setProgressText(GrainSegmentationModifier::tr("Performing polyhedral template matching"));

	// Prepare access to output memory arrays.
	PropertyAccess<int> structuresArray(structures());
	PropertyAccess<FloatType> rmsdArray(rmsd());
	PropertyAccess<Quaternion> orientationsArray(orientations());

	// Mutex is needed to synchronize access to bonds list in parallelized loop.
	std::mutex bondsMutex;

	// Perform analysis on each particle.
	parallelForChunks(_numParticles, *this, [&](size_t startIndex, size_t count, Task& task) {
		// Create a thread-local kernel for the PTM algorithm.
		PTMAlgorithm::Kernel kernel(algorithm);

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

			// Perform the PTM analysis for the current particle.
			PTMAlgorithm::StructureType type = kernel.identifyStructure(index, cachedNeighbors, nullptr);

			// Store identification result in the output property array.
			structuresArray[index] = type;

			int numNeighbors = 0;
			if(type == PTMAlgorithm::OTHER) {
				rmsdArray[index] = -1.0; // Store invalid RMSD value to exclude it from the RMSD histogram.

				kernel.resetNeighbors(index, cachedNeighbors);

				// Don't need more than 8 nearest neighbors to establish connectivity between non-crystalline atoms.
				numNeighbors = std::min(8, kernel.numGoodNeighbors());
			}
			else {
				numNeighbors = ptm_num_nbrs[type];
				rmsdArray[index] = kernel.rmsd();
				orientationsArray[index] = kernel.orientation().normalized();

				if(_rmsdCutoff != 0.0 && kernel.rmsd() > _rmsdCutoff) {
					// Mark atom as OTHER, but still store its RMSD value. It'll be needed to build the RMSD histogram below.
					structuresArray[index] = PTMAlgorithm::OTHER;
				}
			}

			for(int j = 0; j < numNeighbors; j++) {

				size_t neighborIndex = kernel._env.atom_indices[j + 1];

				// Create a bond to the neighbor, but skip every other bond to create just one bond per particle pair.
				if(index < neighborIndex)
					threadlocalNeighborBonds.push_back({index, neighborIndex});

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

		// Append thread-local bonds to global bonds list.
		std::lock_guard<std::mutex> lock(bondsMutex);
		_neighborBonds.insert(_neighborBonds.end(), threadlocalNeighborBonds.cbegin(), threadlocalNeighborBonds.cend());
	});
	if(isCanceled())
		return false;

	// Determine histogram bin size based on maximum RMSD value.
	const size_t numHistogramBins = 100;
	_rmsdHistogram = std::make_shared<PropertyStorage>(numHistogramBins, PropertyStorage::Int64, 1, 0, GrainSegmentationModifier::tr("Count"), true, DataTable::YProperty);
	FloatType rmsdHistogramBinSize = (_numParticles != 0) ? (FloatType(1.01) * *boost::max_element(rmsdArray) / numHistogramBins) : 0.01;
	if(rmsdHistogramBinSize <= 0) rmsdHistogramBinSize = 1;
	_rmsdHistogramRange = rmsdHistogramBinSize * numHistogramBins;

	// Bin RMSD values.
	PropertyAccess<qlonglong> histogramCounts(_rmsdHistogram);
	for(FloatType& rmsd : rmsdArray) {
		if(rmsd >= 0.0) {
			int binIndex = rmsd / rmsdHistogramBinSize;
			histogramCounts[binIndex]++;
		}
		else rmsd = 0.0;
	}
	if(isCanceled())
		return false;

	return !isCanceled();
}

/******************************************************************************
* Calculates the disorientation angle for each graph edge (i.e. bond).
******************************************************************************/
bool GrainSegmentationEngine1::computeDisorientationAngles()
{
	// Compute disorientation angles associated with the neighbor graph edges.
	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - misorientation calculation"));
	ConstPropertyAccess<int> structuresArray(structures());
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());

	parallelFor(_neighborBonds.size(), *this, [&](size_t bondIndex) {
		NeighborBond& bond = _neighborBonds[bondIndex];
		bond.disorientation = std::numeric_limits<FloatType>::max();

		int a = bond.a;
		int b = bond.b;
		if (structuresArray[b] < structuresArray[a]) {
			std::swap(a, b);
		}

		if(structuresArray[a] == structuresArray[b]) {

			int structureType = structuresArray[a];
			const Quaternion& qA = orientationsArray[a];
			const Quaternion& qB = orientationsArray[b];

			double orientA[4] = { qA.w(), qA.x(), qA.y(), qA.z() };
			double orientB[4] = { qB.w(), qB.x(), qB.y(), qB.z() };
			if(structureType == PTMAlgorithm::SC || structureType == PTMAlgorithm::FCC || structureType == PTMAlgorithm::BCC || structureType == PTMAlgorithm::CUBIC_DIAMOND)
				bond.disorientation = (FloatType)ptm::quat_disorientation_cubic(orientA, orientB);
			else if(structureType == PTMAlgorithm::HCP || structureType == PTMAlgorithm::HEX_DIAMOND || structureType == PTMAlgorithm::GRAPHENE)
				bond.disorientation = (FloatType)ptm::quat_disorientation_hcp_conventional(orientA, orientB);

            bond.disorientation = qRadiansToDegrees(bond.disorientation);
		}
#if 0
		else if(structuresArray[a] == PTMAlgorithm::FCC && structuresArray[b] == PTMAlgorithm::HCP) {

			const Quaternion& qA = orientationsArray[a];
			const Quaternion& qB = orientationsArray[b];

			double orientA[4] = { qA.w(), qA.x(), qA.y(), qA.z() };
			double orientB[4] = { qB.w(), qB.x(), qB.y(), qB.z() };

			double map_hcp_to_fcc[49][4] = {{0.11591690,  0.3647052, 0.27984814,  0.88047624},
											{0.45576804, -0.5406251, 0.70455634, -0.06000300}};
			for (int i=0;i<2;i++) {
				double testB[4];
				ptm::quat_rot(orientB, map_hcp_to_fcc[i], testB);
				FloatType disorientation = (FloatType)ptm::quat_disorientation_cubic(orientA, testB);
				bond.disorientation = std::min(bond.disorientation, disorientation);
			}
		}
#endif
	});
	if(isCanceled()) return false;

	// Sort graph edges by disorientation.
	boost::sort(_neighborBonds, [](NeighborBond& a, NeighborBond& b) {
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
        std::vector<NeighborBond>& neighborBonds, ConstPropertyAccess<int>& structuresArray,
        std::vector<Quaternion>& qsum, DisjointSet& uf)
{
	size_t progress = 0;
	for(const NeighborBond& edge : neighborBonds) {
        size_t pa = uf.find(edge.a);
        size_t pb = uf.find(edge.b);
		if(pa != pb && isCrystallineBond(structuresArray, edge)) {
			size_t parent = uf.merge(pa, pb);
			size_t child = (parent == pa) ? pb : pa;
			FloatType disorientation = calculate_disorientation(structuresArray[parent], qsum[parent], qsum[child]);
			OVITO_ASSERT(edge.a < edge.b);
			_dendrogram.emplace_back(parent, child, edge.disorientation, disorientation, 1, qsum[parent]);

			// Update progress indicator.
			if((progress++ % 1024) == 0) {
				if(!incrementProgressValue(1024)) 
					return false;
			}
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
	ConstPropertyAccess<int> structuresArray(structures());
	if(_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic || _algorithmType == GrainSegmentationModifier::NodePairSamplingManual) {

	    setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	    setProgressValue(0);
	    setProgressMaximum(neighborBonds().size());

    	size_t progress = 0;
		for (auto edge: neighborBonds()) {
			if (isCrystallineBond(structuresArray, edge) && edge.disorientation < _misorientationThreshold) {
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
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());
	std::vector<Quaternion> qsum(orientationsArray.cbegin(), orientationsArray.cend());
	DisjointSet uf(_numParticles);
	_dendrogram.resize(0);

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	setProgressValue(0);
	setProgressMaximum(_numParticles);  //TODO: make this num. crystalline particles

	if(_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic || _algorithmType == GrainSegmentationModifier::NodePairSamplingManual) {
		node_pair_sampling_clustering(graph, structuresArray, qsum);
	}
	else {
		minimum_spanning_tree_clustering(neighborBonds(), structuresArray, qsum, uf);
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

	// Create PropertyStorage objects for the output plot.
	PropertyAccess<FloatType> mergeDistanceArray = _mergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), false, DataTable::XProperty);
	PropertyAccess<FloatType> mergeSizeArray = _mergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Merge size"), false, DataTable::YProperty);

	// Generate output data plot points from dendrogram data.
	FloatType* mergeDistanceIter = mergeDistanceArray.begin();
	FloatType* mergeSizeIter = mergeSizeArray.begin();
	for(const DendrogramNode& node : _dendrogram) {
		if(node.size >= _minPlotSize) {
			*mergeDistanceIter++ = std::log(node.distance);
			*mergeSizeIter++ = node.size;
		}
	}

	if(_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic) {
		_suggestedMergingThreshold = calculate_threshold_suggestion();
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
	if(_engine1->_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic) {
		mergingThreshold = _engine1->suggestedMergingThreshold();
    }

	const std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram = _engine1->_dendrogram;

	if(_engine1->_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic) {
        auto graph = _engine1->graph;
        for (int i=dendrogram.size() - 1;i>=0;i--) {
            auto node = &dendrogram[i];
            if (1) {
                graph.reinstate_edge(node->a, node->b);
                printf("%e %e\n", node->distance, graph.wnode[node->a] * graph.wnode[node->b] / graph.adj[node->a][node->b]);
            }
        }
    }


	PropertyAccess<Quaternion> orientationsArray(_engine1->orientations());
	std::vector<Quaternion> meanOrientation(orientationsArray.cbegin(), orientationsArray.cend());

	// Iterate through merge list until distance cutoff is met.
	DisjointSet uf(_numParticles);
	auto node = dendrogram.cbegin();
	for(; node != dendrogram.cend(); ++node) {
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
	ConstPropertyAccess<int> structuresArray(_engine1->structures());
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

	// Build list of orphan atoms.
	std::vector<size_t> orphanAtoms;
	for(size_t i = 0; i < _numParticles; i++) {
		if(atomClustersArray[i] == 0)
			orphanAtoms.push_back(i);
	}

	/// The bonds connecting neighboring non-crystalline atoms.
	std::vector<ParticleIndexPair> noncrystallineBonds;
	for (auto nb: _engine1->neighborBonds()) {
        if (atomClustersArray[nb.a] == 0 || atomClustersArray[nb.b] == 0) {
            // Add bonds for both atoms
            noncrystallineBonds.push_back({(qlonglong)nb.a, (qlonglong)nb.b});
            noncrystallineBonds.push_back({(qlonglong)nb.b, (qlonglong)nb.a});
        }
    }
	if(isCanceled())
		return false;

    boost::stable_sort(noncrystallineBonds);

	// Add orphan atoms to the grains.
	setProgressMaximum(orphanAtoms.size());
	size_t oldOrphanCount = orphanAtoms.size();
	for(;;) {
		std::vector<size_t> newlyAssignedClusters(orphanAtoms.size(), 0);
		for(size_t i = 0; i < orphanAtoms.size(); i++) {
			if(isCanceled()) return false;

			size_t index = orphanAtoms[i];

			// Get the range of bonds adjacent to the current atom.
			auto bondsRange = boost::range::equal_range(noncrystallineBonds, ParticleIndexPair{{(qlonglong)index,0}},
				[](const ParticleIndexPair& a, const ParticleIndexPair& b) { return a[0] < b[0]; });

			// Find the closest cluster atom in the neighborhood (using PTM ordering).
			for(const ParticleIndexPair& bond : boost::make_iterator_range(bondsRange.first, bondsRange.second)) {
				OVITO_ASSERT(bond[0] == index);

				auto neighborIndex = bond[1];
				if(neighborIndex == std::numeric_limits<size_t>::max()) break;
				auto grain = atomClustersArray[neighborIndex];
				if(grain != 0) {
					newlyAssignedClusters[i] = grain;
					break;
				}
			}
		}

		// Assign atoms to closest cluster and compress orphan list.
		size_t newOrphanCount = 0;
		for(size_t i = 0; i < orphanAtoms.size(); i++) {
			atomClustersArray[orphanAtoms[i]] = newlyAssignedClusters[i];
			if(newlyAssignedClusters[i] == 0) {
				orphanAtoms[newOrphanCount++] = orphanAtoms[i];
			}
			else {
				grainSizeArray[newlyAssignedClusters[i] - 1]++;
				if(!incrementProgressValue()) return false;
			}
		}

		orphanAtoms.resize(newOrphanCount);
		if(newOrphanCount == oldOrphanCount)
			break;
		oldOrphanCount = newOrphanCount;
	}

	return !isCanceled();
}

}	// End of namespace
}	// End of namespace
