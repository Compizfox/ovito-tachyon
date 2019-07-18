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

#pragma once


#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/particles/modifier/analysis/StructureIdentificationModifier.h>

#include <plugins/particles/Particles.h>
#include <plugins/particles/objects/ParticlesObject.h>
#include <plugins/particles/objects/BondsObject.h>

#include <plugins/particles/modifier/analysis/ptm/PTMAlgorithm.h>
#include <boost/optional/optional.hpp>


namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

/*
 * Computation engine of the GrainSegmentationModifier, which decomposes a polycrystalline microstructure into individual grains.
 */
class GrainSegmentationEngine : public StructureIdentificationModifier::StructureIdentificationEngine
{
public:
	/// Constructor.
	GrainSegmentationEngine(
			ParticleOrderingFingerprint fingerprint, ConstPropertyPtr positions, const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, ConstPropertyPtr selection,
			FloatType rmsdCutoff, FloatType misorientationThreshold,
			int minGrainAtomCount);

	/// Performs the computation.
	virtual void perform() override;

	/// This method is called by the system to free memory and release any working data after the 
	/// computation has been successfully completed.
	virtual void cleanup() override {
		_neighborLists.reset();
		decltype(_distanceSortedAtoms){}.swap(_distanceSortedAtoms);
		decltype(_clusterOrientations){}.swap(_clusterOrientations);
		decltype(_clusterSizes){}.swap(_clusterSizes);
		StructureIdentificationEngine::cleanup();
	}

	/// Injects the computed results into the data pipeline.
	virtual void emitResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;


	/// Returns the per-atom RMSD values computed by the PTM algorithm.
	const PropertyPtr& rmsd() const { return _rmsd; }

	/// Returns the computed RMSD histogram data.
	const QVector<int>& rmsdHistogramData() const { return _rmsdHistogramData; }

	/// Returns the bin size of the RMSD histogram.
	FloatType rmsdHistogramBinSize() const { return _rmsdHistogramBinSize; }

	/// Replaces the stored RMSD histogram data.
	void setRmsdHistogram(QVector<int> counts, FloatType rmsdHistogramBinSize) {
		_rmsdHistogramData = std::move(counts);
		_rmsdHistogramBinSize = rmsdHistogramBinSize;
	}

	/// Returns the computed per-particle lattice orientations.
	const PropertyPtr& orientations() const { return _orientations; }

	/// Returns the particle to cluster assignment.
	const PropertyPtr& atomClusters() const { return _atomClusters; }

	/// Returns the bonds generated between neighboring lattice atoms.
	const PropertyPtr& latticeNeighborBonds() const { return _latticeNeighborBonds; }

	/// Returns the PBC shift vectors of bonds.
	const PropertyPtr& bondPBCShiftVectors() const { return _bondPBCShiftVectors; }

	/// Allocates the bonds arrays.
	void allocateBonds(size_t count) { 
		_latticeNeighborBonds = BondsObject::OOClass().createStandardStorage(count, BondsObject::TopologyProperty, false);
		_bondPBCShiftVectors = BondsObject::OOClass().createStandardStorage(count, BondsObject::PeriodicImageProperty, true);
		_neighborDisorientationAngles = std::make_shared<PropertyStorage>(count, PropertyStorage::Float, 1, 0, QStringLiteral("Disorientation"), false);
	}

	/// Returns the computed disorientation angles between neighboring lattice atoms.
	const PropertyPtr& neighborDisorientationAngles() const { return _neighborDisorientationAngles; }

private:

	/// Performs the PTM algorithm. Determines the local structure type and the local lattice orientation.
	bool identifyAtomicStructures();

	/// Builds the graph of neighbor atoms and calculates the misorientation angle for each graph edge (i.e. bond).
	bool buildNeighborGraph();

	/// Computes the average lattice orientation of each cluster.
	bool calculateAverageClusterOrientations();

	/// Builds grains by iterative region merging
	bool regionMerging();

	/// Randomizes cluster IDs for testing purposes (giving more color contrast).
	bool randomizeClusterIDs();

	/// Merges any orphan atoms into the closest cluster.
	bool mergeOrphanAtoms();

private:

	/// The cutoff parameter used by the PTM algorithm.
	FloatType _rmsdCutoff;

	/// The merging criterion threshold.
	FloatType _mergingThreshold;

	/// The minimum number of crystalline atoms per grain.	
	int _minGrainAtomCount;

	/// Stores the list of neighbors of each lattice atom.
	PropertyPtr _neighborLists;

	/// Stores indices of crystalline atoms sorted by value of distance transform.
	std::vector<size_t> _distanceSortedAtoms;

	/// Counts the number of clusters
	size_t _numClusters = 0;

	/// Stores the average lattice orientation of each cluster.
	std::vector<Quaternion> _clusterOrientations;
	
	/// Stores the number of atoms in each cluster.
	std::vector<qlonglong> _clusterSizes;

	/// The per-atom RMSD values computed by the PTM algorithm.
	const PropertyPtr _rmsd;

	/// Histogram of the RMSD values computed by the PTM algorithm.
	QVector<int> _rmsdHistogramData;

	/// Bin size of the RMSD histogram.
	FloatType _rmsdHistogramBinSize;

	/// The computed per-particle lattice orientations.
	PropertyPtr _orientations;

	/// The particle to cluster assignment.
	PropertyPtr _atomClusters;

	/// The bonds generated between neighboring lattice atoms.
	PropertyPtr _latticeNeighborBonds;

	/// The PBC shift vectors of bonds.
	PropertyPtr _bondPBCShiftVectors;

	/// The computed disorientation angles between neighboring lattice atoms.
	PropertyPtr _neighborDisorientationAngles;
};

}	// End of namespace
}	// End of namespace
}	// End of namespace
