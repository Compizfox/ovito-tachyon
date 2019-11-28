#ifndef GRAPH_H
#define GRAPH_H

#include <core/utilities/FloatType.h>


namespace Ovito { namespace Plugins { namespace CrystalAnalysis {


class GraphEdge
{
public:
	GraphEdge(size_t _a, size_t _b, FloatType _w, size_t _superCluster)
		: a(_a), b(_b), w(_w), superCluster(_superCluster) {}

	size_t a;
	size_t b;
	FloatType w;
	size_t superCluster;
};

class DendrogramNode
{
public:
	DendrogramNode(size_t _a, size_t _b, FloatType _d, size_t _size)
		: a(_a), b(_b), d(_d), size(_size) {}

	DendrogramNode()
		: a(0), b(0), d(-INFINITY), size(0) {}

	size_t a;
	size_t b;
	FloatType d;
	size_t size;
};

class Graph
{
public:
	size_t next = 0;
	double wtotal = 0;
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

		double check = dmin * wnode[a] / wtotal;
		if (check != check) {
			printf("e. bad number: %lu %lu %e %e %e\n", a, vmin, dmin, wnode[a], wtotal);
			exit(3);
		}

		return std::make_tuple(dmin * wnode[a] / wtotal, vmin);
	}

	void add_node(size_t u) {
		rep[u] = u;//next++;
		next = u + 1;
		snode[u] = 1;
		wnode[u] = 0;
	}

	void add_edge(size_t u, size_t v, double w, bool update) {

		auto it = adj.find(u);
		if (it == adj.end()) {
			if (!update) {printf("update error\n"); exit(3);}
			add_node(u);
		}

		it = adj.find(v);
		if (it == adj.end()) {
			if (!update) {printf("update error\n"); exit(3);}
			add_node(v);
		}

		(adj[u])[v] = w;
		(adj[v])[u] = w;

		if (update) {
			wnode[u] += w;
			wnode[v] += w;
			if (u != v) {
				wtotal += w;
			}
		}
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

	void contract_edge(size_t a, size_t b) {

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
	}
};


}	// End of namespace
}	// End of namespace
}	// End of namespace

#endif

