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
			int minGrainAtomCount) :
	StructureIdentificationModifier::StructureIdentificationEngine(std::move(fingerprint), positions, simCell, std::move(typesToIdentify), std::move(selection)),
	_rmsdCutoff(rmsdCutoff),
	_mergingThreshold(mergingThreshold),
	_minGrainAtomCount(std::max(minGrainAtomCount, 1)),
	_rmsd(std::make_shared<PropertyStorage>(positions->size(), PropertyStorage::Float, 1, 0, QStringLiteral("RMSD"), false)),
	_orientations(ParticlesObject::OOClass().createStandardStorage(positions->size(), ParticlesObject::OrientationProperty, true)),
	_atomClusters(ParticlesObject::OOClass().createStandardStorage(positions->size(), ParticlesObject::ClusterProperty, true))
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
	if(!regionMerging()) return;

	//if(!randomizeClusterIDs()) return;
	if(!calculateAverageClusterOrientations()) return;

	// For final output, convert edge disorientation angles from radians to degrees.
	for(FloatType& angle : neighborDisorientationAngles()->floatRange())
		angle *= FloatType(180) / FLOATTYPE_PI;

	//if(!mergeOrphanAtoms()) return;
}

/******************************************************************************
* Performs the PTM algorithm. Determines the local structure type and the 
* local lattice orientation at each atomic site.
******************************************************************************/
bool GrainSegmentationEngine::identifyAtomicStructures()
{
	boost::optional<PTMAlgorithm> _algorithm;

	_algorithm.emplace();
	_algorithm->setRmsdCutoff(0.0); // Note: We do our own RMSD threshold filtering in postProcessStructureTypes().


	// Specify the structure types the PTM should look for.
	for(int i = 0; i < typesToIdentify().size() && i < PTMAlgorithm::NUM_STRUCTURE_TYPES; i++) {
		_algorithm->setStructureTypeIdentification(static_cast<PTMAlgorithm::StructureType>(i), typesToIdentify()[i]);
	}

	// Initialize the algorithm object.
	if(!_algorithm->prepare(*positions(), cell(), selection().get(), task().get()))
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

					_neighborLists->setInt64Component(index, j, kernel._atom_indices[j + 1]);

#if 0
//TODO: put back in
					const Vector3& neighborVector = neighQuery.results()[mapping[j + 1] - 1].delta;
					// Check if neighbor vector spans more than half of a periodic simulation cell.
					for(size_t dim = 0; dim < 3; dim++) {
						if(cell().pbcFlags()[dim]) {
							if(std::abs(cell().inverseMatrix().prodrow(neighborVector, dim)) >= FloatType(0.5)+FLOATTYPE_EPSILON) {
								static const QString axes[3] = { QStringLiteral("X"), QStringLiteral("Y"), QStringLiteral("Z") };
								throw Exception(GrainSegmentationModifier::tr("Simulation box is too short along cell vector %1 (%2) to perform analysis. "
										"Please extend it first using the 'Replicate' modifier.").arg(dim+1).arg(axes[dim]));
							}
						}
					}
#endif
				}
			}
			else {
				rmsd()->setFloat(index, 0);

#if 0
//TODO: put back in.  use topological ordering if atom is disordered.
				// Store neighbor list.  Don't need more than 12 neighbors for defect atoms.
				numNeighbors = std::min(numNeighbors, PTM_MAX_NBRS);
				numNeighbors = std::min(numNeighbors, 12);
				OVITO_ASSERT(numNeighbors <= _neighborLists->componentCount());
				for(int j = 0; j < numNeighbors; j++) {
					_neighborLists->setInt64Component(index, j, neighQuery.results()[j].index);
				}
#endif
			}
		}
	});
	if(task()->isCanceled() || positions()->size() == 0)
		return false;

	// Determine histogram bin size based on maximum RMSD value.
	QVector<int> rmsdHistogramData(100, 0);
	FloatType rmsdHistogramBinSize = FloatType(1.01) * *std::max_element(rmsd()->constDataFloat(), rmsd()->constDataFloat() + rmsd()->size());
	rmsdHistogramBinSize /= rmsdHistogramData.size();
	if(rmsdHistogramBinSize <= 0) rmsdHistogramBinSize = 1;

	// Build RMSD histogram.
	for(size_t index = 0; index < structures()->size(); index++) {
		if(structures()->getInt(index) != PTMAlgorithm::OTHER) {
			OVITO_ASSERT(rmsd()->getFloat(index) >= 0);
			int binIndex = rmsd()->getFloat(index) / rmsdHistogramBinSize;
			if(binIndex < rmsdHistogramData.size())
				rmsdHistogramData[binIndex]++;
		}
	}
	setRmsdHistogram(std::move(rmsdHistogramData), rmsdHistogramBinSize);

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

	// Count number of bonds to create.
	size_t bondCount = 0;
	for(size_t index = 0; index < positions()->size(); index++) {
		for(size_t c = 0; c < _neighborLists->componentCount(); c++) {
			if(_neighborLists->getInt64Component(index, c) == -1) break;
			if(_neighborLists->getInt64Component(index, c) < index) continue;
			bondCount++;
		}
	}
	allocateBonds(bondCount);

	qlonglong* bondTopology = latticeNeighborBonds()->dataInt64();
	Vector3I* pbcShift = bondPBCShiftVectors()->dataVector3I();
	for(size_t index = 0; index < positions()->size(); index++) {
		if(!task()->incrementProgressValue()) return false;

		for(size_t c = 0; c < _neighborLists->componentCount(); c++) {
			auto neighborIndex = _neighborLists->getInt64Component(index, c);
			if(neighborIndex == -1) break;

			// Skip every other atom pair so that we don't create the same bond twice.
			if(_neighborLists->getInt64Component(index, c) < index) continue;

			// Determine PBC bond shift using minimum image convention.
			Vector3 delta = positions()->getPoint3(index) - positions()->getPoint3(neighborIndex);
			for(size_t dim = 0; dim < 3; dim++) {
				if(cell().pbcFlags()[dim])
					(*pbcShift)[dim] = (int)floor(cell().inverseMatrix().prodrow(delta, dim) + FloatType(0.5));
			}
			++pbcShift;
			*bondTopology++ = index;
			*bondTopology++ = neighborIndex;
		}
	}

	// Compute disorientation angles associated with the neighbor graph edges.
	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - misorientation calculation"));
	parallelFor(bondCount, *task(), [this](size_t bondIndex) {
		size_t index1 = latticeNeighborBonds()->getInt64Component(bondIndex, 0);
		size_t index2 = latticeNeighborBonds()->getInt64Component(bondIndex, 1);
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

/******************************************************************************
* Builds grains by iterative region merging
******************************************************************************/
bool GrainSegmentationEngine::regionMerging()
{
	size_t numAtoms = positions()->size();
	_clusterSizes.resize(numAtoms, 1);

	// Disjoint sets data structures.
	std::vector<size_t> ranks(numAtoms, 0);
	std::vector<size_t> parents(numAtoms);
	std::iota(parents.begin(), parents.end(), (size_t)0);

	std::vector<size_t> sizes(numAtoms, 1);
	std::vector<FloatType> weights(numAtoms, 0);

	// Disjoint-sets helper function. Find part of Union-Find
	auto findParent = [&parents](size_t index) {
		// Find root and make root as parent of i (path compression)
		size_t parent = parents[index];
	    while(parent != parents[parent]) {
	    	parent = parents[parent];
	    }
		parents[index] = parent;
	    return parent;
	};

	auto merge = [&findParent, &parents, &ranks, &sizes, &weights](size_t index1, size_t index2) {
		size_t parentA = findParent(index1);
		size_t parentB = findParent(index2);
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
	};


	task()->setProgressText(GrainSegmentationModifier::tr("Grain segmentation - region merging"));
	task()->setProgressValue(0);
	task()->setProgressMaximum(latticeNeighborBonds()->size());

	std::vector< std::tuple< size_t, size_t, FloatType, FloatType > > graph;
	FloatType threshold = 4 / FloatType(180) * FLOATTYPE_PI;

	// Build graph edges
	for(size_t bondIndex = 0; bondIndex < latticeNeighborBonds()->size(); bondIndex++) {
		if(!task()->incrementProgressValue()) return false;
		size_t index1 = latticeNeighborBonds()->getInt64Component(bondIndex, 0);
		size_t index2 = latticeNeighborBonds()->getInt64Component(bondIndex, 1);

		// Skip high-angle edges.
		FloatType disorientation = neighborDisorientationAngles()->getFloat(bondIndex);
		if(disorientation > threshold) continue;

		FloatType deg = disorientation * FloatType(180) / FLOATTYPE_PI;
		FloatType weight = std::exp(-FloatType(1)/3 * deg * deg);		//this is fairly arbitrary but it works well

		graph.emplace_back(std::make_tuple(index1, index2, weight, -1));
	}

	if(task()->isCanceled())
		return false;

	std::vector<bool> hit(numAtoms, 0);
	std::vector< std::tuple< size_t, size_t, FloatType > > contracted;

	bool merged = true;
	while (merged) {
		merged = false;
		std::vector<size_t> components(parents.begin(), parents.end());
		std::vector<FloatType> curWeights(weights.begin(), weights.end());

		double wsum = 0;
		for (auto edge: graph)
			wsum += std::get<2>(edge);

		double cwsum = 0;
		for(size_t i = 0; i < numAtoms; i++)
			if (i == findParent(i))
				cwsum += weights[i];

		FloatType threshold = _mergingThreshold;
		for (size_t i=0;i<graph.size();i++) {	//can't use auto since we want to set tuple elements

			std::tuple< size_t, size_t, FloatType, FloatType > edge = graph[i];

			size_t a = std::get<0>(edge);
			size_t b = std::get<1>(edge);
			FloatType w = std::get<2>(edge);

			int coordinations[PTMAlgorithm::NUM_STRUCTURE_TYPES] = {0, 12, 12, 14, 12, 6, 16, 16, 9};
			int dimensions[PTMAlgorithm::NUM_STRUCTURE_TYPES] = {0, 3, 3, 3, 3, 3, 3, 3, 2};

			int type = structures()->getInt(a);
			int hc = coordinations[type] / 2;
			int dim = dimensions[type];

			FloatType power = (-1. + dim) / dim;
			FloatType za = (hc + w) / (1 + pow(curWeights[a] / hc, power));
			FloatType zb = (hc + w) / (1 + pow(curWeights[b] / hc, power));
			std::get<3>(graph[i]) = std::max(za, zb);
		}

		std::sort(graph.begin(), graph.end(),
			[](std::tuple< size_t, size_t, FloatType, FloatType > a, std::tuple< size_t, size_t, FloatType, FloatType > b)
			{return std::get<3>(b) < std::get<3>(a);});

		std::fill(hit.begin(), hit.end(), 0);

		for (auto edge: graph) {

			size_t a = std::get<0>(edge);
			size_t b = std::get<1>(edge);
			FloatType w = std::get<2>(edge);
			FloatType z = std::get<3>(edge);

			int num_hit = hit[a] + hit[b];
			if (z >= threshold && num_hit < 2) {
				merge(a, b);
				hit[a] = true;
				hit[b] = true;
				merged = true;
			}
		}

		std::map< std::tuple< size_t, size_t >, FloatType > contracted;

		for (auto edge: graph) {

			size_t a = std::get<0>(edge);
			size_t b = std::get<1>(edge);
			FloatType w = std::get<2>(edge);

			size_t parentA = findParent(a);
			size_t parentB = findParent(b);
			if (parentA == parentB) {
				weights[parentA] += w;
			}
			else {
				size_t keymin = parentA, keymax = parentB;
				if (keymin > keymax) {
					std::swap(keymin, keymax);
				}

				std::tuple< size_t, size_t > key = std::make_tuple(keymin, keymax);

				if (contracted.find(key) == contracted.end()) {
					contracted[key] = 0;
				}
				contracted[key] += w;
			}
		}

		graph.clear();
		for (auto it=contracted.begin(); it!=contracted.end(); it++) {

			auto key = it->first;
			size_t a = std::get<0>(key);
			size_t b = std::get<1>(key);
			auto weight = it->second;
			graph.emplace_back(std::make_tuple(a, b, weight, -1));
		}
	}

	for(size_t i = 0; i < numAtoms; i++) {
		_clusterSizes[i] = sizes[i];
	}

	// Relabels the clusters to obtain a contiguous sequence of cluster IDs.	
	std::vector<size_t> clusterRemapping(numAtoms);

	// Assign new consecutive IDs to root clusters.
	_numClusters = 1;
	for(size_t i = 0; i < numAtoms; i++) {
		if(findParent(i) == i) {
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
		clusterRemapping[particleIndex] = clusterRemapping[findParent(particleIndex)];

	// Relabel atoms after cluster IDs have changed.
	_clusterSizes.resize(_numClusters);
	std::fill(_clusterSizes.begin(), _clusterSizes.end(), 0);
	for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {

		size_t grainID = clusterRemapping[particleIndex];
		atomClusters()->setInt64(particleIndex, grainID);
		_clusterSizes[grainID]++;
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

		for(size_t particleIndex = 0; particleIndex < numAtoms; particleIndex++) {

			size_t grainID = atomClusters()->getInt64(particleIndex);
			atomClusters()->setInt64(particleIndex, std::get<0>(lut[grainID]));
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

}	// End of namespace
}	// End of namespace
}	// End of namespace
