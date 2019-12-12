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

#include <plugins/crystalanalysis/CrystalAnalysis.h>

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

class Graph
{
public:
	size_t next = 0;
	std::map<size_t, std::map<size_t, FloatType>> adj;
	std::map<size_t, FloatType> wnode;
	std::map<size_t, size_t> snode;
	std::map<size_t, size_t> rep;

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
		rep[u] = u;
		next = u + 1;
		snode[u] = 1;
		wnode[u] = 0;
	}

	void add_edge(size_t u, size_t v, FloatType w) {

		auto it = adj.find(u);
		if (it == adj.end()) {
			add_node(u);
		}

		it = adj.find(v);
		if (it == adj.end()) {
			add_node(v);
		}

		(adj[u])[v] = w;
		(adj[v])[u] = w;

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

		rep[a] = next++;
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


}	// End of namespace
}	// End of namespace
}	// End of namespace
