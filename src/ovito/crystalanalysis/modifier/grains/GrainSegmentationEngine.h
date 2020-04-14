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

#pragma once


#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/particles/modifier/analysis/StructureIdentificationModifier.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/particles/modifier/analysis/ptm/PTMAlgorithm.h>
#include "DisjointSet.h"
#include "GrainSegmentationModifier.h"

namespace Ovito { namespace CrystalAnalysis {

/*
 * Computation engine of the GrainSegmentationModifier, which decomposes a polycrystalline microstructure into individual grains.
 */
class GrainSegmentationEngine1 : public StructureIdentificationModifier::StructureIdentificationEngine
{
public:

    class Graph
    {
    public:
	    size_t next = 0;
	    std::map<size_t, FloatType> wnode;
	    std::map<size_t, std::map<size_t, FloatType>> adj;
	    std::map<size_t, std::map<size_t, FloatType>> deleted_adj;

	    size_t num_nodes() const {
		    return wnode.size();
	    }

	    size_t next_node() const {
		    return wnode.begin()->first;
	    }

	    std::tuple<FloatType, size_t> nearest_neighbor(size_t a) const {
		    FloatType dmin = std::numeric_limits<FloatType>::max();
		    size_t vmin = std::numeric_limits<size_t>::max();

		    OVITO_ASSERT(adj.find(a) != adj.end());
		    for (const auto& x : adj.find(a)->second) {
			    size_t v = x.first;
			    FloatType weight = x.second;

			    OVITO_ASSERT(v != a); // Graph has self loops.
			    if(v == a) {
				    qWarning() << "Graph has self loops";
				    exit(3);
			    }

			    OVITO_ASSERT(wnode.find(v) != wnode.end());
			    FloatType d = wnode.find(v)->second / weight;
			    OVITO_ASSERT(!std::isnan(d));

			    if (d < dmin) {
				    dmin = d;
				    vmin = v;
			    }
			    else if (d == dmin) {
				    vmin = std::min(vmin, v);
			    }
		    }

		    OVITO_ASSERT(wnode.find(a) != wnode.end());
		    FloatType check = dmin * wnode.find(a)->second;
		    OVITO_ASSERT(!std::isnan(check));

		    return std::make_tuple(dmin * wnode.find(a)->second, vmin);
	    }

	    void add_node(size_t u) {
		    next = u + 1;
		    wnode[u] = 0;
	    }

	    void add_edge(size_t u, size_t v, FloatType w) {

		    auto it_u = adj.find(u);
		    if (it_u == adj.end()) {
			    add_node(u);
			    it_u = adj.emplace(u, std::map<size_t, FloatType>{{{v,w}}}).first;
		    }
		    else it_u->second[v] = w;

		    auto it_v = adj.find(v);
		    if (it_v == adj.end()) {
			    add_node(v);
			    it_v = adj.emplace(v, std::map<size_t, FloatType>{{{u,w}}}).first;
		    }
		    else it_v->second[u] = w;

		    wnode[u] += w;
		    wnode[v] += w;
	    }

	    void remove_node(size_t u) {

		    for (auto const& x: adj[u]) {
			    size_t v = x.first;
			    adj[v].erase(u);
		    }

		    adj.erase(u);
		    wnode.erase(u);
            //deleted_adj[u] = adj[u];
            //adj.erase(u);
	    }

	    void reinstate_node(size_t u) {

            adj[u] = deleted_adj[u];
            deleted_adj.erase(u);

		    for (auto const& x: adj[u]) {
			    size_t v = x.first;
			    FloatType w = x.second;
			    (adj[v])[u] += w;
		    }
	    }

	    size_t contract_edge(size_t a, size_t b) {

		    if (adj[b].size() > adj[a].size()) {
			    std::swap(a, b);
		    }

		    for (auto const& x: adj[b]) {
			    size_t v = x.first;
			    FloatType w = x.second;
                if (v == a) continue;

			    (adj[a])[v] += w;
			    (adj[v])[a] += w;
		    }

		    adj[a].erase(b);
		    wnode[a] += wnode[b];
		    remove_node(b);
		    return a;
	    }

	    size_t reinstate_edge(size_t a, size_t b) {
            // TODO: check component sizes

		    for (auto const& x: deleted_adj[b]) {
			    size_t v = x.first;
			    FloatType w = x.second;
                if (v == a) continue;

			    (adj[a])[v] -= w;
			    (adj[v])[a] -= w;

                // TODO: if edge weight is now ~=zero, remove edge
		    }

		    (adj[a])[b] = (deleted_adj[b])[a];
		    wnode[a] -= wnode[b];
		    reinstate_node(b);
		    return a;
	    }
    };

	/// Represents a single bond connecting two neighboring lattice atoms.
	struct NeighborBond {
		size_t a;
		size_t b;
		FloatType disorientation;
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
	GrainSegmentationEngine1(
			ParticleOrderingFingerprint fingerprint, 
			ConstPropertyPtr positions, 
			const SimulationCell& simCell,
			const QVector<bool>& typesToIdentify, 
			ConstPropertyPtr selection,
			FloatType rmsdCutoff, 
			GrainSegmentationModifier::MergeAlgorithm algorithmType,
			bool outputBonds);

	/// Performs the computation.
	virtual void perform() override;

	/// Injects the computed results into the data pipeline.
	virtual void applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

	/// This method is called by the system whenever a parameter of the modifier changes.
	/// The method can be overriden by subclasses to indicate to the caller whether the engine object should be 
	/// discarded (false) or may be kept in the cache, because the computation results are not affected by the changing parameter (true). 
	virtual bool modifierChanged(const PropertyFieldEvent& event) override {

		// Avoid a recomputation if a parameters changes that does not affect this algorithm stage.
		if(event.field() == &PROPERTY_FIELD(GrainSegmentationModifier::colorParticlesByGrain)
				|| event.field() == &PROPERTY_FIELD(GrainSegmentationModifier::mergingThreshold) 
				|| event.field() == &PROPERTY_FIELD(GrainSegmentationModifier::minGrainAtomCount)
				|| event.field() == &PROPERTY_FIELD(GrainSegmentationModifier::orphanAdoption))
			return true;

		return StructureIdentificationModifier::StructureIdentificationEngine::modifierChanged(event);
	}

	/// Creates another engine that performs the next stage of the computation. 
	virtual std::shared_ptr<Engine> createContinuationEngine(ModifierApplication* modApp, const PipelineFlowState& input) override;

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

	/// Returns the adaptively determined merge threshold.
	FloatType suggestedMergingThreshold() const { return _suggestedMergingThreshold; }

private:

	/// Returns the list of bonds connecting neighboring lattice atoms.
	std::vector<NeighborBond>& neighborBonds() { return _neighborBonds; }

	/// Performs the PTM algorithm. Determines the local structure type and the local lattice orientation.
	bool identifyAtomicStructures();

	/// Calculates the disorientation angle for each graph edge (i.e. bond).
	bool computeDisorientationAngles();

	/// Builds grains by iterative region merging.
	bool determineMergeSequence();

	/// Computes the disorientation angle between two crystal clusters of the given lattice type. 
	/// Furthermore, the function computes the weighted average of the two cluster orientations. 
	/// The norm of the two input quaternions and the output quaternion represents the size of the clusters.
	static FloatType calculate_disorientation(int structureType, Quaternion& qa, const Quaternion& qb);

	// Algorithm types:
	bool minimum_spanning_tree_clustering(std::vector<NeighborBond>& neighborBonds, ConstPropertyAccess<int>& structuresArray, std::vector<Quaternion>& qsum, DisjointSet& uf);
	bool node_pair_sampling_clustering(Graph& graph, ConstPropertyAccess<int>& structuresArray, std::vector<Quaternion>& qsum);

	// Selects a threshold for Node Pair Sampling algorithm
    FloatType calculate_threshold_suggestion();

    // Determines if a bond is crystalline
    bool isCrystallineBond(ConstPropertyAccess<int>& structuresArray, const NeighborBond& bond)
    {
        return structuresArray[bond.a] != PTMAlgorithm::OTHER
               && structuresArray[bond.a] == structuresArray[bond.b];
    }

private:

	const size_t _minPlotSize = 20;

	/// The number of input particles.
	size_t _numParticles;

	/// The cutoff parameter used by the PTM algorithm.
	FloatType _rmsdCutoff;

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
	GrainSegmentationModifier::MergeAlgorithm _algorithmType;

	/// The computed per-particle lattice orientations.
	PropertyPtr _orientations;

	/// The bonds connecting neighboring lattice atoms.
	std::vector<NeighborBond> _neighborBonds;

	/// Controls the output of neighbor bonds to the data pipeline for visualization purposes.
	bool _outputBondsToPipeline;

	// Dendrogram as list of cluster merges.
	std::vector<DendrogramNode> _dendrogram;

	/// The adaptively computed merge threshold.
	FloatType _suggestedMergingThreshold = 0;

	friend class GrainSegmentationEngine2;
};

/*
 * Computation engine of the GrainSegmentationModifier, which decomposes a polycrystalline microstructure into individual grains.
 */
class GrainSegmentationEngine2 : public AsynchronousModifier::Engine
{
public:

	/// Constructor.
	GrainSegmentationEngine2(
			std::shared_ptr<GrainSegmentationEngine1> engine1,
			FloatType mergingThreshold, 
			bool adoptOrphanAtoms, 
			size_t minGrainAtomCount) :
		_engine1(std::move(engine1)),
		_numParticles(_engine1->_numParticles),
		_mergingThreshold(mergingThreshold),
		_adoptOrphanAtoms(adoptOrphanAtoms),
		_minGrainAtomCount(minGrainAtomCount),
		_atomClusters(ParticlesObject::OOClass().createStandardStorage(_numParticles, ParticlesObject::ClusterProperty, true)) {}
	
	/// Performs the computation.
	virtual void perform() override;

	/// Injects the computed results into the data pipeline.
	virtual void applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

	/// This method is called by the system whenever a parameter of the modifier changes.
	/// The method can be overriden by subclasses to indicate to the caller whether the engine object should be 
	/// discarded (false) or may be kept in the cache, because the computation results are not affected by the changing parameter (true). 
	virtual bool modifierChanged(const PropertyFieldEvent& event) override {

		// Avoid a recomputation if a parameters changes that does not affect the algorithm's results.
		if(event.field() == &PROPERTY_FIELD(GrainSegmentationModifier::colorParticlesByGrain))
			return true; // Indicate that the stored results are not affected by the parameter change.

		return Engine::modifierChanged(event);
	}

	/// Returns the array storing the cluster ID of each particle.
	const PropertyPtr& atomClusters() const { return _atomClusters; }

private:

	/// Merges any orphan atoms into the closest cluster.
	bool mergeOrphanAtoms();

private:

	/// Pointer to the first algorithm stage.
	std::shared_ptr<GrainSegmentationEngine1> _engine1;

	/// The number of input particles.
	size_t _numParticles;

	/// The particle to cluster assignment.
	PropertyPtr _atomClusters;

	/// Counts the number of clusters
	size_t _numClusters = 1;

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

	/// The user-defined merge threshold.
	FloatType _mergingThreshold;

	/// The minimum number of atoms a grain must have.
	size_t _minGrainAtomCount;

	/// Contrals the adoption of orphan atoms after the grains have been formed.
	bool _adoptOrphanAtoms;
};

}	// End of namespace
}	// End of namespace
