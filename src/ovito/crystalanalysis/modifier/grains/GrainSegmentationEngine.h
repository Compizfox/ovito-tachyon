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

#pragma once


#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/particles/modifier/analysis/StructureIdentificationModifier.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/particles/modifier/analysis/ptm/PTMAlgorithm.h>
#include "DisjointSet.h"

namespace Ovito { namespace CrystalAnalysis {

/*
 * Computation engine of the GrainSegmentationModifier, which decomposes a polycrystalline microstructure into individual grains.
 */
class GrainSegmentationEngine : public StructureIdentificationModifier::StructureIdentificationEngine
{
public:

	/// Represents a single bond connecting two neighboring lattice atoms.
	struct NeighborBond {
		size_t a;
		size_t b;
		FloatType disorientation;
		FloatType weight;
		size_t superCluster;
	};

	struct DendrogramNode {

		DendrogramNode() = default;
		DendrogramNode(size_t _a, size_t _b, FloatType _distance, FloatType _disorientation, size_t _size, Quaternion _orientation)
			: a(_a), b(_b), distance(_distance), disorientation(_disorientation), size(_size), orientation(_orientation) {}

		size_t a = 0;
		size_t b = 0;
		FloatType distance = std::numeric_limits<FloatType>::lowest();
		FloatType disorientation = std::numeric_limits<FloatType>::lowest();
		size_t size = 0;
		Quaternion orientation;
	};

	/// Constructor.
	GrainSegmentationEngine(
			ParticleOrderingFingerprint fingerprint, ConstPropertyPtr positions, const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, ConstPropertyPtr selection,
			FloatType rmsdCutoff, bool algorithmType, bool outputBonds);

	/// Performs the computation.
	virtual void perform() override;

	/// Injects the computed results into the data pipeline.
	virtual void emitResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

	/// Returns the per-atom RMSD values computed by the PTM algorithm.
	const PropertyPtr& rmsd() const { return _rmsd; }

	/// Returns the RMSD value range of the histogram.
	FloatType rmsdHistogramRange() const { return _rmsdHistogramRange; }

	/// Returns the histogram of computed RMSD values.
	const PropertyPtr& rmsdHistogram() const { return _rmsdHistogram; }

	// Returns the merge distances for the scatter plot
	const PropertyPtr& mergeDistance() const { return _mergeDistance; }

	// Returns the merge sizes for the scatter plot
	const PropertyPtr& mergeSize() const { return _mergeSize; }

	/// Returns the computed per-particle lattice orientations.
	const PropertyPtr& orientations() const { return _orientations; }

	/// Returns the array storing the cluster ID of each particle.
	const PropertyPtr& atomClusters() const { return _atomClusters; }

private:

	/// Returns the list of bonds connecting neighboring lattice atoms.
	std::vector<NeighborBond>& neighborBonds() { return _neighborBonds; }

	/// Performs the PTM algorithm. Determines the local structure type and the local lattice orientation.
	bool identifyAtomicStructures();

	/// Calculates the disorientation angle for each graph edge (i.e. bond).
	bool computeDisorientationAngles();

	/// Groups lattice atoms with similar orientations into superclusters.
	bool formSuperclusters();

	/// Builds grains by iterative region merging.
	bool determineMergeSequence();

	/// Executes precomputed merge steps up to the threshold value set by the user.
	void executeMergeSequence(int minGrainAtomCount, FloatType mergingThreshold, bool adoptOrphanAtoms);

	/// Merges any orphan atoms into the closest cluster.
	bool mergeOrphanAtoms();

	/// Computes the disorientation angle between two crystal clusters of the given lattice type. 
	/// Furthermore, the function computes the weighted average of the two cluster orientations. 
	/// The norm of the two input quaternions and the output quaternion represents the size of the clusters.
	static FloatType calculate_disorientation(int structureType, Quaternion& qa, const Quaternion& qb);

	// Algorithm types:
	bool minimum_spanning_tree_clustering(boost::iterator_range<std::vector<NeighborBond>::iterator> edgeRange, DendrogramNode* dendrogram, int structureType, std::vector<Quaternion>& qsum, DisjointSet& uf);
	bool node_pair_sampling_clustering(boost::iterator_range<std::vector<NeighborBond>::const_iterator> edgeRange, DendrogramNode* dendrogram, int structureType, std::vector<Quaternion>& qsum, FloatType totalWeight);

	// Statistical functions
	FloatType calculate_median(std::vector< FloatType >& data);
	std::vector< FloatType > theil_sen_estimator(size_t num_samples, std::vector< std::tuple< FloatType, FloatType> >& data,
												 FloatType& gradient, FloatType& intercept);

    FloatType calculate_threshold_suggestion();

private:

	const size_t _minPlotSize = 20;

	/// The number of input particles.
	size_t _numParticles;

	/// Supercluster ID of each particle.
	std::vector<size_t> _atomSuperclusters;

	/// Counts the number of superclusters.
	size_t _numSuperclusters = 0;

	/// Stores the number of particle in each supercluster.
	std::vector<size_t> _superclusterSizes;

	/// The cutoff parameter used by the PTM algorithm.
	FloatType _rmsdCutoff;

	/// Counts the number of clusters
	size_t _numClusters = 0;

	/// The per-atom RMSD values computed by the PTM algorithm.
	const PropertyPtr _rmsd;

	/// Histogram of the RMSD values computed by the PTM algorithm.
	PropertyPtr _rmsdHistogram;

	/// The value range of the RMSD histogram.
	FloatType _rmsdHistogramRange;

	// The merge distances
	PropertyPtr _mergeDistance;

	// The merge sizes
	PropertyPtr _mergeSize;

	// The linkage criterion used in the merge algorithm
	bool _algorithmType;

	/// The computed per-particle lattice orientations.
	PropertyPtr _orientations;

	/// The particle to cluster assignment.
	PropertyPtr _atomClusters;

	/// The bonds connecting neighboring lattice atoms.
	std::vector<NeighborBond> _neighborBonds;

	/// The bonds connecting neighboring non-crystalline atoms.
	std::vector<ParticleIndexPair> _noncrystallineBonds;

	/// Controls the output of neighbor bonds to the data pipeline for visualization purposes.
	bool _outputBondsToPipeline;

	// A hardcoded cutoff used for defining superclusters
	const FloatType _misorientationThreshold = qDegreesToRadians(4.0);

	// Dendrogram as list of cluster merges.
	std::vector<DendrogramNode> _dendrogram;

	/// Tells for each atom (being an orphan) what its parent atom is. 
	std::vector<std::atomic<size_t>> _orphanParentAtoms;

	// The output list of grain IDs.
	PropertyPtr _grainIds;

	// The output list of grain sizes.
	PropertyPtr _grainSizes;

	/// The output list of per-grain structure types.
	PropertyPtr _grainStructureTypes;

	/// The output list of colors assigned to grains.
	PropertyPtr _grainColors;

	/// The output list of mean grain orientations.
	PropertyPtr _grainOrientations;
};

}	// End of namespace
}	// End of namespace
