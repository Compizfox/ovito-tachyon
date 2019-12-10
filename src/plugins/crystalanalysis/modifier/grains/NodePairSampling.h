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
	std::map<size_t, std::map<size_t, double> > adj;
	std::map<size_t, double> wnode;
	std::map<size_t, size_t> snode;
	std::map<size_t, size_t> rep;

	Graph() {}

	size_t num_nodes() {
		return adj.size();
	}

	size_t next_node() {
		return adj.begin()->first;
	}

	std::tuple< double, size_t> nearest_neighbor(size_t a) {
		double dmin = INFINITY;
		size_t vmin = (size_t)(-1);

		for (auto const& x: adj[a]) {
			size_t v = x.first;
			double weight = x.second;

			if (v == a) {
				printf("graph has self loops\n");
				exit(3);
			}

			double d = wnode[v] / weight;
			if (d != d) {
				printf("bad number: %lu %lu %e %e\n", a, v, d, weight);
				exit(3);
			}

			if (d < dmin) {
				dmin = d;
				vmin = v;
			}
			else if (d == dmin) {
				vmin = std::min(vmin, v);
			}
		}

		double check = dmin * wnode[a];
		if (check != check) {
			printf("e. bad number: %lu %lu %e %e\n", a, vmin, dmin, wnode[a]);
			exit(3);
		}

		return std::make_tuple(dmin * wnode[a], vmin);
	}

	void add_node(size_t u) {
		rep[u] = u;//next++;
		next = u + 1;
		snode[u] = 1;
		wnode[u] = 0;
	}

	void add_edge(size_t u, size_t v, double w) {

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
			double w = x.second;

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
