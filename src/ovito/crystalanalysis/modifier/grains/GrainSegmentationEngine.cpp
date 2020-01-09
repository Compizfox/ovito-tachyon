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
#include <ovito/stdobj/table/DataTable.h>
#include <ovito/particles/util/NearestNeighborFinder.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include "GrainSegmentationEngine.h"
#include "GrainSegmentationModifier.h"
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
	_numParticles(positions->size()),
	_rmsdCutoff(rmsdCutoff),
	_algorithmType(algorithmType),
	_minGrainAtomCount(std::max(minGrainAtomCount, 1)),
	_rmsd(std::make_shared<PropertyStorage>(_numParticles, PropertyStorage::Float, 1, 0, QStringLiteral("RMSD"), false)),
	_orientations(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::OrientationProperty, true)),
	_atomClusters(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::ClusterProperty, true)),
	_atomSuperclusters(_numParticles),
	_orphanAdoption(orphanAdoption),
	_outputBondsToPipeline(outputBonds)
{
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine::perform()
{
	if(_numParticles == 0)
		return;	// No input particles, nothing to do.

	// Grain segmentation algorithm:
	if(!identifyAtomicStructures()) return;
	if(!computeDisorientationAngles()) return;
	if(!formSuperclusters()) return;
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

	if(!algorithm.prepare(*positions(), cell(), selection(), task().get()))
		return false;

	task()->setProgressValue(0);
	task()->setProgressMaximum(_numParticles);
	task()->setProgressText(GrainSegmentationModifier::tr("Pre-calculating neighbor ordering"));

	// Pre-order neighbors of each particle
	std::vector<uint64_t> cachedNeighbors(_numParticles);
	parallelForChunks(_numParticles, *task(), [this, &cachedNeighbors, &algorithm](size_t startIndex, size_t count, Task& task) {
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
			kernel.precacheNeighbors(index, &cachedNeighbors[index]);
		}
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	task()->setProgressValue(0);
	task()->setProgressText(GrainSegmentationModifier::tr("Performing polyhedral template matching"));

	// Prepare access to output memory arrays.
	PropertyAccess<int> structuresArray(structures());
	PropertyAccess<FloatType> rmsdArray(rmsd());
	PropertyAccess<Quaternion> orientationsArray(orientations());

	// Mutex is needed to synchronize access to bonds list in parallelized loop.
	std::mutex bondsMutex;

	// Perform analysis on each particle.
	parallelForChunks(_numParticles, *task(), [&](size_t startIndex, size_t count, Task& task) {
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

			if(type != PTMAlgorithm::OTHER) {
				rmsdArray[index] = kernel.rmsd();

				// Store computed local lattice orientation in the output property array.
				if(orientationsArray)
					orientationsArray[index] = kernel.orientation().normalized();

				// Store neighbor list for later use.
				for(int j = 0; j < ptm_num_nbrs[type]; j++) {

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
			else {
				rmsdArray[index] = 0;
#if 0
				// Don't need more than 12 nearest neighbors to establish connectiveity between defect atoms.
				int numNeighbors = std::min(12, kernel.numNearestNeighbors());

				// Store neighbor list.
				for(int j = 0; j < numNeighbors; j++) {
					_neighborLists[index][j] = kernel.getNearestNeighbor(j).index;
				}
#endif
			}
		}

		// Append thread-local bonds to global bonds list.
		std::lock_guard<std::mutex> lock(bondsMutex);
		_neighborBonds.insert(_neighborBonds.end(), threadlocalNeighborBonds.cbegin(), threadlocalNeighborBonds.cend());
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	// Determine histogram bin size based on maximum RMSD value.
	const size_t numHistogramBins = 100;
	_rmsdHistogram = std::make_shared<PropertyStorage>(numHistogramBins, PropertyStorage::Int64, 1, 0, GrainSegmentationModifier::tr("Count"), true, DataTable::YProperty);
	FloatType rmsdHistogramBinSize = FloatType(1.01) * *boost::max_element(rmsdArray) / numHistogramBins;
	if(rmsdHistogramBinSize <= 0) rmsdHistogramBinSize = 1;
	_rmsdHistogramRange = rmsdHistogramBinSize * numHistogramBins;

	// Bin RMSD values.
	PropertyAccess<qlonglong> histogramCounts(_rmsdHistogram);
	for(size_t index = 0; index < _numParticles; index++) {
		if(structuresArray[index] != PTMAlgorithm::OTHER) {
			int binIndex = rmsdArray[index] / rmsdHistogramBinSize;
			if(binIndex < numHistogramBins)
				histogramCounts[binIndex]++;
		}
	}

	// Reset the identified structure type to OTHER for particles having an RMSD value above the threshold.
	if(_rmsdCutoff > 0) {
		for(size_t index = 0; index < _numParticles; index++) {
			if(structuresArray[index] != PTMAlgorithm::OTHER) {
				if(rmsdArray[index] > _rmsdCutoff)
					structuresArray[index] = PTMAlgorithm::OTHER;
			}
		}
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Calculates the disorientation angle for each graph edge (i.e. bond).
******************************************************************************/
bool GrainSegmentationEngine::computeDisorientationAngles()
{
	// Compute disorientation angles associated with the neighbor graph edges.
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - misorientation calculation"));
	ConstPropertyAccess<int> structuresArray(structures());
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());

	parallelFor(_neighborBonds.size(), *task(), [&](size_t bondIndex) {
		NeighborBond& bond = _neighborBonds[bondIndex];
		bond.disorientation = std::numeric_limits<FloatType>::infinity();

		int structureTypeA = structuresArray[bond.a];
		int structureTypeB = structuresArray[bond.b];
		if(structureTypeA == structureTypeB) {

			int structureType = structureTypeA;
			const Quaternion& qA = orientationsArray[bond.a];
			const Quaternion& qB = orientationsArray[bond.b];

			double orientA[4] = { qA.w(), qA.x(), qA.y(), qA.z() };
			double orientB[4] = { qB.w(), qB.x(), qB.y(), qB.z() };

			if(structureType == PTMAlgorithm::SC || structureType == PTMAlgorithm::FCC || structureType == PTMAlgorithm::BCC || structureType == PTMAlgorithm::CUBIC_DIAMOND)
				bond.disorientation = (FloatType)ptm::quat_disorientation_cubic(orientA, orientB);
			else if(structureType == PTMAlgorithm::HCP || structureType == PTMAlgorithm::HEX_DIAMOND || structureType == PTMAlgorithm::GRAPHENE)
				bond.disorientation = (FloatType)ptm::quat_disorientation_hcp_conventional(orientA, orientB);
		}
	});

	return !task()->isCanceled();
}

/******************************************************************************
* Groups lattice atoms with similar orientations into superclusters
******************************************************************************/
bool GrainSegmentationEngine::formSuperclusters()
{
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - supercluster merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(neighborBonds().size());

	// Group connected particles having similar lattice orientations into superclusters.
	DisjointSet uf(_numParticles);
	size_t progress = 0;
	for(const NeighborBond& bond : neighborBonds()) {
		if(!task()->setProgressValueIntermittent(progress++)) return false;

		// Skip high-angle edges.
		if(bond.disorientation > _misorientationThreshold) continue;

		uf.merge(bond.a, bond.b);
	}

	// Relabel the superclusters to obtain a contiguous sequence of cluster IDs.
	_superclusterSizes.resize(1);
	std::vector<size_t> superclusterRemapping(_numParticles);
	// Assign new consecutive IDs to root superclusters.
	for(size_t i = 0; i < _numParticles; i++) {
		if(uf.find(i) == i) {
			superclusterRemapping[i] = _superclusterSizes.size();
			_superclusterSizes.push_back(uf.nodesize(i));
		}
	}
	_numSuperclusters = _superclusterSizes.size();

	if(task()->isCanceled()) return false;

	// Supercluster 0 contains all atoms that are not part of a regular supercluster.
	_superclusterSizes[0] = _numParticles - std::accumulate(_superclusterSizes.begin() + 1, _superclusterSizes.end(), (size_t)0);

	// Determine supercluster IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		superclusterRemapping[particleIndex] = superclusterRemapping[uf.find(particleIndex)];

	if(task()->isCanceled()) return false;

	// Relabel atoms after cluster IDs have changed.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		_atomSuperclusters[particleIndex] = superclusterRemapping[particleIndex];

	return !task()->isCanceled();
}

/******************************************************************************
* Computes the disorientation angle between two crystal clusters of the 
* given lattice type. Furthermore, the function computes the weighted average
* of the two cluster orientations. The norm of the two input quaternions 
* and the output quaternion represents the size of the clusters.
******************************************************************************/
FloatType GrainSegmentationEngine::calculate_disorientation(int structureType, Quaternion& qa, const Quaternion& qb)
{
	FloatType qa_norm = qa.norm();
	FloatType qb_norm = qb.norm();
	double qtarget[4] = { qa.w()/qa_norm, qa.x()/qa_norm, qa.y()/qa_norm, qa.z()/qa_norm };
	double q[4]       = { qb.w()/qb_norm, qb.x()/qb_norm, qb.y()/qb_norm, qb.z()/qb_norm };

	// Convert structure type back to PTM representation
	int type = 0;
	if(structureType == PTMAlgorithm::FCC) type = PTM_MATCH_FCC;
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
bool GrainSegmentationEngine::minimum_spanning_tree_clustering(
		boost::iterator_range<std::vector<NeighborBond>::iterator> edgeRange,
		DendrogramNode* dendrogram, int structureType, std::vector<Quaternion>& qsum, DisjointSet& uf)
{
	// Sort graph edges by weight.
	boost::sort(edgeRange, [](NeighborBond& a, NeighborBond& b) {
		return a.weight < b.weight;
	});
	if(task()->isCanceled()) return false;

	size_t progress = 0;
	for(const NeighborBond& edge : edgeRange) {
		if(uf.find(edge.a) != uf.find(edge.b)) {
			size_t parent = uf.merge(edge.a, edge.b);
			size_t child = (parent == edge.a) ? edge.b : edge.a;
			FloatType disorientation  = calculate_disorientation(structureType, qsum[parent], qsum[child]);
			OVITO_ASSERT(edge.a < edge.b);
			*dendrogram++ = DendrogramNode(edge.a, edge.b, edge.weight, disorientation, 1);

			// Update progress indicator.
			if((progress++ % 1024) == 0) {
				if(!task()->incrementProgressValue(1024)) 
					return false;
			}
		}
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine::determineMergeSequence()
{
	// There is not much to do if there are no crystalline atoms at all (just supercluster 0).
	if(_numSuperclusters == 1) {
		_numClusters = 1;
		boost::copy(_atomSuperclusters, PropertyAccess<qlonglong>(atomClusters()).begin());
		return true;
	}

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(neighborBonds().size());

	// Build initial graph.
	FloatType totalWeight = 0;
	std::vector<size_t> bondCount(_numSuperclusters, 0); // Number of bonds in each supercluster.

	size_t progress = 0;
	for(NeighborBond& bond : neighborBonds()) {
		if(!task()->setProgressValueIntermittent(progress++)) return false;

		bond.superCluster = _atomSuperclusters[bond.a];

		// Skip high-angle edges.
		if(bond.superCluster != 0 && bond.disorientation <= _misorientationThreshold) {
			// Convert disorientations to graph weights.
			FloatType deg = qRadiansToDegrees(bond.disorientation);
			if(_algorithmType == 0)
				bond.weight = std::exp(-FloatType(1)/3 * deg * deg);	// This is fairly arbitrary but it works well.
			else
				bond.weight = deg;

			totalWeight += bond.weight;
		}
		else {
			bond.weight = 0;
			bond.superCluster = 0;
		}
		bondCount[bond.superCluster]++;
	}

	// setting the total weight to 1 is an effective multi-frame normalization
	totalWeight = 1;

	// Group graph edges by supercluster.
	boost::sort(neighborBonds(), [](NeighborBond& a, NeighborBond& b) {
		return a.superCluster < b.superCluster;
	});

	// Compute the start index in the global edge list for each supercluster.
	std::vector<size_t> bondStart(_numSuperclusters);
	bondStart[0] = 0;
	std::partial_sum(bondCount.cbegin(), bondCount.cend() - 1, bondStart.begin() + 1);

	if(task()->isCanceled())
		return false;

	// Allocate memory for the dendrograms.
	size_t dendrogramSize = 0;
	std::vector<size_t> dendrogramOffsets(_numSuperclusters);
	dendrogramOffsets[0] = 0;
	for(size_t sc = 1; sc < _numSuperclusters; sc++) {
		dendrogramOffsets[sc] = dendrogramSize;
		dendrogramSize += _superclusterSizes[sc] - 1;	// Number of edges in a tree = number of vertices - 1.
	}
	_dendrogram.resize(dendrogramSize);

	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(dendrogramSize);

	// Build dendrograms.
	ConstPropertyAccess<int> structuresArray(structures());
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());
	std::vector<Quaternion> qsum(orientationsArray.cbegin(), orientationsArray.cend());
	DisjointSet uf(_numParticles);

	// Parallelize dendrogram computation over superclusters.
	parallelFor(_numSuperclusters - 1, [&](size_t sc) {
		sc++;
		if(_superclusterSizes[sc] >= _minGrainAtomCount) {
			size_t start = bondStart[sc];
			size_t count = bondCount[sc];
			size_t index = dendrogramOffsets[sc];
			int structureType = structuresArray[neighborBonds()[start].a];

			if(_algorithmType == 0) {
				node_pair_sampling_clustering(
					boost::make_iterator_range_n(neighborBonds().cbegin() + start, count), 
					&_dendrogram[index], structureType, qsum, totalWeight);
			}
			else {
				minimum_spanning_tree_clustering(
					boost::make_iterator_range_n(neighborBonds().begin() + start, count), 
					&_dendrogram[index], structureType, qsum, uf);
			}
		}
	});
	if(task()->isCanceled())
		return false;

	// Sort dendrogram entries by distance.
	boost::sort(_dendrogram, [](const DendrogramNode& a, const DendrogramNode& b) { return a.distance < b.distance; });

	if(task()->isCanceled())
		return false;

	// Scan through the entire merge list to determine merge sizes.
	size_t numPlot = 0;
	uf.clear();
	for(DendrogramNode& node : _dendrogram) {
		size_t sa = uf.nodesize(uf.find(node.a));
		size_t sb = uf.nodesize(uf.find(node.b));
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

	// Create PropertyStorage objects for the output plot.
	PropertyAccess<FloatType> mergeDistanceArray = _mergeDistance = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Log merge distance"), false, DataTable::XProperty);
	PropertyAccess<FloatType> mergeSizeArray = _mergeSize = std::make_shared<PropertyStorage>(numPlot, PropertyStorage::Float, 1, 0, GrainSegmentationModifier::tr("Merge size"), false, DataTable::YProperty);

	// Generate output data plot points from dendrogram data.
	FloatType* mergeDistanceIter = mergeDistanceArray.begin();
	FloatType* mergeSizeIter = mergeSizeArray.begin();
	for(const DendrogramNode& node : _dendrogram) {
		if(node.size > 0) {
			*mergeDistanceIter++ = std::log(node.distance);
			*mergeSizeIter++ = node.size;
		}
	}

	return !task()->isCanceled();
}

/******************************************************************************
* Executes precomputed merge steps up to the threshold value set by the user.
******************************************************************************/
void GrainSegmentationEngine::executeMergeSequence(int minGrainAtomCount, FloatType mergingThreshold)
{
	// Iterate through merge list until distance cutoff is met.
	DisjointSet uf(_numParticles);
	for(const DendrogramNode& node : _dendrogram) {
		if(std::log(node.distance) >= mergingThreshold)
			break;

		uf.merge(node.a, node.b);
	}

	// Relabels the clusters to obtain a contiguous sequence of cluster IDs.
	std::vector<size_t> clusterRemapping(_numParticles);

	// Assign new consecutive IDs to root clusters.
	_numClusters = 1;
	ConstPropertyAccess<int> structuresArray(structures());
	std::vector<int> clusterStructureTypes;
	for(size_t i = 0; i < _numParticles; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(uf.nodesize(i) < minGrainAtomCount) {
				clusterRemapping[i] = 0;
			}
			else {
				clusterRemapping[i] = _numClusters;
				_numClusters++;
				clusterStructureTypes.push_back(structuresArray[i]);
			}
		}
	}

	// Allocate and fill output array storing the grain IDs (1-based identifiers). 
	_grainIds =  std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int64, 1, 0, QStringLiteral("Grain Identifier"), false, DataTable::XProperty);
	boost::algorithm::iota_n(PropertyAccess<qlonglong>(_grainIds).begin(), size_t(1), _grainIds->size());

	// Allocate output array storing the grain sizes.
	_grainSizes = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int64, 1, 0, QStringLiteral("Grain Size"), true, DataTable::YProperty);

	// Allocate output array storing the structure type of grains.
	_grainStructureTypes = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Int, 1, 0, QStringLiteral("Structure Type"), false);
	boost::copy(clusterStructureTypes, PropertyAccess<int>(_grainStructureTypes).begin());

	// Allocate output array with each grain's unique color.
	// Fill it with random color values (using constant random seed to keep it reproducible).
	_grainColors = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 3, 0, QStringLiteral("Color"), false);
	std::default_random_engine rng(1);
	std::uniform_real_distribution<FloatType> uniform_dist(0, 1);
	boost::generate(PropertyAccess<Color>(_grainColors), [&]() { return Color::fromHSV(uniform_dist(rng), 1.0 - uniform_dist(rng) * 0.5, 1.0 - uniform_dist(rng) * 0.3); });

	// Allocate output array storing the mean lattice orientation of grains (represented by a quaternion).
	_grainOrientations = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 4, 0, QStringLiteral("Orientation"), true);

	// Determine new IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		clusterRemapping[particleIndex] = clusterRemapping[uf.find(particleIndex)];

	// Relabel atoms after cluster IDs have changed.
	// Also count the number of atoms in each cluster.
	PropertyAccess<qlonglong> atomClustersArray(atomClusters());
	PropertyAccess<qlonglong> grainSizeArray(_grainSizes);
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++) {
		size_t gid = clusterRemapping[particleIndex];
		atomClustersArray[particleIndex] = gid;
		if(gid != 0) grainSizeArray[gid - 1]++;
	}

	// Reorder grains by size (large to small).
	if(_numClusters > 1) {

		// Determine the index remapping for reordering the grain list by size.
		std::vector<size_t> mapping(_numClusters - 1);
		std::iota(mapping.begin(), mapping.end(), size_t(0));
		std::sort(mapping.begin(), mapping.end(), [&](size_t a, size_t b) {
			return grainSizeArray[a] > grainSizeArray[b];
		});

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

		// Invert the grain index map. 

		std::vector<size_t> inverseMapping(_numClusters);
		inverseMapping[0] = 0; // Keep cluster ID 0 in place.
		for(size_t i = 1; i < _numClusters; i++)
			inverseMapping[mapping[i-1]+1] = i;

		// Remap per-particle grain IDs.

		for(auto& id : atomClustersArray)
			id = inverseMapping[id];
	}
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

#if 0
/******************************************************************************
* Merges any orphan atoms into the closest cluster.
******************************************************************************/
bool GrainSegmentationEngine::mergeOrphanAtoms()
{
	PropertyAccess<qlonglong> atomClustersArray(atomClusters());
	ConstPropertyAccess<Point3> positionsArray(positions());

	// Build list of orphan atoms.
	std::vector<size_t> orphanAtoms;
	for(size_t i = 0; i < positions()->size(); i++) {
		if(atomClustersArray[i] == 0)
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
			for(size_t c = 0; c < _neighborLists[orphanAtoms[i]].size(); c++) {
				auto neighborIndex = _neighborLists[orphanAtoms[i]][c];
				if(neighborIndex == std::numeric_limits<size_t>::max()) break;
				auto clusterId = atomClustersArray[neighborIndex];
				if(clusterId == 0) continue;

				// Determine interatomic vector using minimum image convention.
				Vector3 delta = cell().wrapVector(positionsArray[neighborIndex] - positionsArray[orphanAtoms[i]]);
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
			atomClustersArray[orphanAtoms[i]] = newlyAssignedClusters[i];
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
#endif

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
