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
GrainSegmentationEngine::GrainSegmentationEngine(
			ParticleOrderingFingerprint fingerprint, ConstPropertyPtr positions, const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, ConstPropertyPtr selection,
			FloatType rmsdCutoff, GrainSegmentationModifier::MergeAlgorithm algorithmType, bool outputBonds) :
	StructureIdentificationModifier::StructureIdentificationEngine(std::move(fingerprint), positions, simCell, std::move(typesToIdentify), std::move(selection)),
	_numParticles(positions->size()),
	_rmsdCutoff(rmsdCutoff),
	_algorithmType(algorithmType),
	_rmsd(std::make_shared<PropertyStorage>(_numParticles, PropertyStorage::Float, 1, 0, QStringLiteral("RMSD"), false)),
	_orientations(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::OrientationProperty, true)),
	_atomClusters(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::ClusterProperty, true)),
	_atomSuperclusters(_numParticles),
	_outputBondsToPipeline(outputBonds)
{
}

/******************************************************************************
* The grain segmentation algorithm.
******************************************************************************/
void GrainSegmentationEngine::perform()
{
	// Grain segmentation algorithm:
	if(!identifyAtomicStructures()) return;
	if(!computeDisorientationAngles()) return;
	if(!formSuperclusters()) return;
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
bool GrainSegmentationEngine::identifyAtomicStructures()
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
bool GrainSegmentationEngine::computeDisorientationAngles()
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

	return !isCanceled();
}

/******************************************************************************
* Groups lattice atoms with similar orientations into superclusters
******************************************************************************/
bool GrainSegmentationEngine::formSuperclusters()
{
	ConstPropertyAccess<int> structuresArray(structures());

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - supercluster merging"));
	setProgressValue(0);
	setProgressMaximum(neighborBonds().size());

	// Group connected particles having similar lattice orientations into superclusters.
	DisjointSet uf(_numParticles);
	size_t progress = 0;
	for(const NeighborBond& bond : neighborBonds()) {
		if(!setProgressValueIntermittent(progress++)) return false;

		// Skip high-angle edges.
		if(bond.disorientation > _misorientationThreshold) continue;

		OVITO_ASSERT(structuresArray[bond.a] != PTMAlgorithm::OTHER);
		//OVITO_ASSERT(structuresArray[bond.a] == structuresArray[bond.b]);	 // TODO: fix this for stacking faults

		uf.merge(bond.a, bond.b);
	}

	// Relabel the superclusters to obtain a contiguous sequence of cluster IDs.
	_superclusterSizes.resize(1);
	std::vector<size_t> superclusterRemapping(_numParticles);
	// Assign new consecutive IDs to root superclusters.
	for(size_t i = 0; i < _numParticles; i++) {
		if(uf.find(i) == i && structuresArray[i] != PTMAlgorithm::OTHER) {
			superclusterRemapping[i] = _superclusterSizes.size();
			_superclusterSizes.push_back(uf.nodesize(i));
		}
	}
	_numSuperclusters = _superclusterSizes.size();

	if(isCanceled()) return false;

	// Supercluster 0 contains all atoms that are not part of a regular supercluster.
	_superclusterSizes[0] = _numParticles - std::accumulate(_superclusterSizes.begin() + 1, _superclusterSizes.end(), (size_t)0);

	// Determine supercluster IDs for non-root clusters.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		superclusterRemapping[particleIndex] = superclusterRemapping[uf.find(particleIndex)];

	if(isCanceled()) return false;

	// Relabel atoms after cluster IDs have changed.
	for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
		_atomSuperclusters[particleIndex] = superclusterRemapping[particleIndex];

	return !isCanceled();
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
	double q[4]	   = { qb.w()/qb_norm, qb.x()/qb_norm, qb.y()/qb_norm, qb.z()/qb_norm };

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
	if(isCanceled()) return false;

	size_t progress = 0;
	for(const NeighborBond& edge : edgeRange) {
		if(uf.find(edge.a) != uf.find(edge.b)) {
			size_t pa = uf.find(edge.a);
			size_t pb = uf.find(edge.b);
			size_t parent = uf.merge(pa, pb);
			size_t child = (parent == pa) ? pb : pa;
			FloatType disorientation = calculate_disorientation(structureType, qsum[parent], qsum[child]);
			OVITO_ASSERT(edge.a < edge.b);
			*dendrogram++ = DendrogramNode(edge.a, edge.b, edge.weight, disorientation, 1, qsum[parent]);

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
bool GrainSegmentationEngine::determineMergeSequence()
{
	// There is not much to do if there are no crystalline atoms at all (just supercluster 0).
	if(_numSuperclusters == 1) {
		_numClusters = 1;
		boost::copy(_atomSuperclusters, PropertyAccess<qlonglong>(atomClusters()).begin());
		return true;
	}

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - building graph"));
	setProgressValue(0);
	setProgressMaximum(neighborBonds().size());

	// Build initial graph.
	FloatType totalWeight = 0;
	std::vector<size_t> bondCount(_numSuperclusters, 0); // Number of bonds in each supercluster.
	ConstPropertyAccess<int> structuresArray(structures());

	size_t progress = 0;
	for(NeighborBond& bond : neighborBonds()) {
		if(!setProgressValueIntermittent(progress++)) return false;

		bond.superCluster = _atomSuperclusters[bond.a];

		// Skip high-angle edges.
		if(bond.superCluster != 0 && bond.disorientation <= _misorientationThreshold) {
			// Convert disorientations to graph weights.
			FloatType deg = bond.disorientation;
			if(_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic || _algorithmType == GrainSegmentationModifier::NodePairSamplingManual) {
				bond.weight = std::exp(-FloatType(1)/3 * deg * deg);	// This is fairly arbitrary but it works well.
				if (structuresArray[bond.a] != structuresArray[bond.b]) {
					bond.weight /= 2;
				}
			}
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

	// Group graph edges by supercluster.
	boost::sort(neighborBonds(), [](const NeighborBond& a, const NeighborBond& b) {
		return a.superCluster < b.superCluster;
	});

	// Compute the start index in the global edge list for each supercluster.
	std::vector<size_t> bondStart(_numSuperclusters);
	bondStart[0] = 0;
	std::partial_sum(bondCount.cbegin(), bondCount.cend() - 1, bondStart.begin() + 1);

	if(isCanceled())
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

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	setProgressValue(0);
	setProgressMaximum(dendrogramSize);

	// Build dendrograms.
	ConstPropertyAccess<Quaternion> orientationsArray(orientations());
	std::vector<Quaternion> qsum(orientationsArray.cbegin(), orientationsArray.cend());
	DisjointSet uf(_numParticles);

	// Parallelize dendrogram computation over superclusters.
	parallelFor(_numSuperclusters - 1, [&](size_t sc) {
		sc++;
		size_t start = bondStart[sc];
		size_t count = bondCount[sc];
		if(count == 0) return;
		size_t index = dendrogramOffsets[sc];
		int structureType = structuresArray[neighborBonds()[start].a];

		if(_algorithmType == GrainSegmentationModifier::NodePairSamplingAutomatic || _algorithmType == GrainSegmentationModifier::NodePairSamplingManual) {
			// setting the total weight to 1 is an effective multi-frame normalization
			node_pair_sampling_clustering(
				boost::make_iterator_range_n(neighborBonds().cbegin() + start, count), 
				&_dendrogram[index], structureType, qsum, 1);
		}
		else {
			minimum_spanning_tree_clustering(
				boost::make_iterator_range_n(neighborBonds().begin() + start, count), 
				&_dendrogram[index], structureType, qsum, uf);
		}
	});
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

size_t count = 0;
for(size_t particleIndex = 0; particleIndex < _numParticles; particleIndex++)
	count += structuresArray[particleIndex] == PTMAlgorithm::OTHER ? 0 : 1;

if (fout)
	fprintf(fout, "%lu %e\n", count, totalWeight);
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
* Executes precomputed merge steps up to the threshold value set by the user.
******************************************************************************/
void GrainSegmentationEngine::executeMergeSequence(int minGrainAtomCount, FloatType mergingThreshold, bool adoptOrphanAtoms)
{
	PropertyAccess<Quaternion> orientationsArray(orientations());
	std::vector<Quaternion> meanOrientation(orientationsArray.cbegin(), orientationsArray.cend());

	// Iterate through merge list until distance cutoff is met.
	DisjointSet uf(_numParticles);
	auto node = _dendrogram.cbegin();
	for(; node != _dendrogram.cend(); ++node) {
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
	ConstPropertyAccess<int> structuresArray(structures());
	std::vector<int> clusterStructureTypes;
	std::vector<Quaternion> clusterOrientations;
	for(size_t i = 0; i < _numParticles; i++) {
		if(uf.find(i) == i) {
			// If the cluster's size is below the threshold, dissolve the cluster.
			if(uf.nodesize(i) < minGrainAtomCount || structuresArray[i] == PTMAlgorithm::OTHER) {
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

	// Continue iterating through merge list to merge dissolved grains into adjacent grains.
	// Makes sure that sub-critical clusters that are isolated (i.e. not connected to any super-critical cluster)
	// remain dissolved even if they grow in size by eating other sub-critical clusters.
	for(; node != _dendrogram.cend(); ++node) {
		size_t clusterA = uf.find(node->a);
		size_t clusterB = uf.find(node->b);
		
		// Don't merge two super-critical clusters.
		if(clusterRemapping[clusterA] != 0 && clusterRemapping[clusterB] != 0)
			continue;

		// Merge the two clusters.
		uf.merge(node->a, node->b);
		size_t parent = uf.find(node->a);
		meanOrientation[parent] = node->orientation;

		// When merging a super-critical and a sub-critical cluster:
		if(clusterRemapping[clusterB] == 0) {
			clusterRemapping[parent] = clusterRemapping[clusterA];
		}
		else if(clusterRemapping[clusterA] == 0) {
			clusterRemapping[parent] = clusterRemapping[clusterB];
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
	_grainColors = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 3, 0, QStringLiteral("Color"), false, 0, QStringList() << QStringLiteral("R") << QStringLiteral("G") << QStringLiteral("B"));
	std::default_random_engine rng(1);
	std::uniform_real_distribution<FloatType> uniform_dist(0, 1);
	boost::generate(PropertyAccess<Color>(_grainColors), [&]() { return Color::fromHSV(uniform_dist(rng), 1.0 - uniform_dist(rng) * 0.8, 1.0 - uniform_dist(rng) * 0.5); });

	// Allocate output array storing the mean lattice orientation of grains (represented by a quaternion).
	_grainOrientations = std::make_shared<PropertyStorage>(_numClusters - 1, PropertyStorage::Float, 4, 0, QStringLiteral("Orientation"), true, 0, QStringList() << QStringLiteral("X") << QStringLiteral("Y") << QStringLiteral("Z") << QStringLiteral("W"));
	boost::copy(clusterOrientations, PropertyAccess<Quaternion>(_grainOrientations).begin());

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

		// Adopt orphan atoms.
		if(adoptOrphanAtoms)
			mergeOrphanAtoms();
	}
}

/******************************************************************************
* Merges any orphan atoms into the closest cluster.
******************************************************************************/
bool GrainSegmentationEngine::mergeOrphanAtoms()
{
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
	for (auto nb: neighborBonds()) {
        if (atomClustersArray[nb.a] == 0 || atomClustersArray[nb.b] == 0) {
            // Add bonds for both atoms
            noncrystallineBonds.push_back({(qlonglong)nb.a, (qlonglong)nb.b});
            noncrystallineBonds.push_back({(qlonglong)nb.b, (qlonglong)nb.a});
        }
    }

    boost::stable_sort(noncrystallineBonds);

	setProgressText(GrainSegmentationModifier::tr("Grain segmentation - merging orphan atoms"));
	setProgressValue(0);
	setProgressMaximum(orphanAtoms.size());

	// Add orphan atoms to the grains.
	size_t oldOrphanCount = orphanAtoms.size();
	for(;;) {
		std::vector<size_t> newlyAssignedClusters(orphanAtoms.size(), 0);
		for(size_t i = 0; i < orphanAtoms.size(); i++) {
			//if(task()->isCanceled()) return false;

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
				//if(!task()->incrementProgressValue()) return false;
			}
		}

		orphanAtoms.resize(newOrphanCount);
		if(newOrphanCount == oldOrphanCount)
			break;
		oldOrphanCount = newOrphanCount;
	}

	return true;//!task()->isCanceled();
}

}	// End of namespace
}	// End of namespace
