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
#include "GrainSegmentationEngine.h"

namespace Ovito { namespace CrystalAnalysis {

namespace {

class Graph
{
public:
	size_t next = 0;
	std::map<size_t, std::map<size_t, FloatType>> adj;
	std::map<size_t, FloatType> wnode;
	std::map<size_t, size_t> snode;

	size_t num_nodes() const {
		return adj.size();
	}

	size_t next_node() const {
		return adj.begin()->first;
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
		snode[u] = 1;
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
		snode.erase(u);
	}

	size_t contract_edge(size_t a, size_t b) {

		if (adj[b].size() > adj[a].size()) {
			std::swap(a, b);
		}

		adj[a].erase(b);
		adj[b].erase(a);

		for (auto const& x: adj[b]) {
			size_t v = x.first;
			FloatType w = x.second;

			(adj[a])[v] += w;
			(adj[v])[a] += w;
		}

		wnode[a] += wnode[b];
		snode[a] += snode[b];
		remove_node(b);
		return a;
	}
};

} // End of anonymous namespace

/******************************************************************************
* Clustering using pair sampling algorithm.
******************************************************************************/
bool GrainSegmentationEngine1::node_pair_sampling_clustering(
	boost::iterator_range<std::vector<NeighborBond>::const_iterator> edgeRange, 
	DendrogramNode* dendrogram, int structureType, std::vector<Quaternion>& qsum, FloatType totalWeight)
{
	Graph graph;
	for(const NeighborBond& edge : edgeRange) {
        // Calculate edge weight based on disorientation. This is fairly arbitrary but it works well.
		FloatType weight = std::exp(-FloatType(1)/3 * edge.disorientation * edge.disorientation);
		//if (structuresArray[bond.a] != structuresArray[bond.b]) {
		//	bond.weight /= 2;
		//}

		graph.add_edge(edge.a, edge.b, weight);
	}

	size_t progress = 0;
	size_t n = graph.num_nodes();
	while(graph.num_nodes()) {

		// nearest-neighbor chain
		size_t node = graph.next_node();
		OVITO_ASSERT(node != std::numeric_limits<size_t>::max());

		std::vector<size_t> chain{node};
		while(chain.size()) {

			size_t a = chain.back();
			chain.pop_back();
			OVITO_ASSERT(a != std::numeric_limits<size_t>::max());

			FloatType d;
			size_t b;
			std::tie(d, b) = graph.nearest_neighbor(a);
			if(b == std::numeric_limits<size_t>::max()) {
				OVITO_ASSERT(chain.size() == 0);
				// Remove the connected component
				size_t sa = graph.snode[a];
				graph.remove_node(a);
			}
			else if(chain.size()) {
				size_t c = chain.back();
				chain.pop_back();

				if(b == c) {
					size_t size = graph.snode[a] + graph.snode[b];
					OVITO_ASSERT(size != 0);
					size_t parent = graph.contract_edge(a, b);
					size_t child = (parent == a) ? b : a;

					FloatType disorientation = calculate_disorientation(structureType, qsum[parent], qsum[child]);

					*dendrogram++ = DendrogramNode(std::min(a, b), std::max(a, b), d / totalWeight, disorientation, size, qsum[parent]);

					// Update progress indicator.
					if((progress++ % 1024) == 0) {
						if(!incrementProgressValue(1024)) 
							return false;
					}
				}
				else {
					OVITO_ASSERT(a != std::numeric_limits<size_t>::max());
					OVITO_ASSERT(b != std::numeric_limits<size_t>::max());
					OVITO_ASSERT(c != std::numeric_limits<size_t>::max());
					chain.push_back(c);
					chain.push_back(a);
					chain.push_back(b);
				}
			}
			else {
				OVITO_ASSERT(a != std::numeric_limits<size_t>::max());
				OVITO_ASSERT(b != std::numeric_limits<size_t>::max());
				chain.push_back(a);
				chain.push_back(b);
			}
		}
	}

	return true;
}

}	// End of namespace
}	// End of namespace
