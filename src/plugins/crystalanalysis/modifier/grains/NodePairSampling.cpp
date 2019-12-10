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

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include "GrainSegmentationEngine.h"
#include "NodePairSampling.h"

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {


bool GrainSegmentationEngine::node_pair_sampling_clustering(std::vector< GraphEdge >& initial_graph, size_t start, size_t end,
								FloatType totalWeight, DendrogramNode* dendrogram, int structureType, std::vector< Quaternion >& qsum)
{
	Graph graph;
	for (size_t i=start;i<end;i++) {
		auto edge = initial_graph[i];
		FloatType deg = edge.w;
		FloatType weight = std::exp(-FloatType(1)/3 * deg * deg);	//this is fairly arbitrary but it works well
		graph.add_edge(edge.a, edge.b, weight);
	}

	//std::vector< std::tuple< size_t, size_t > > components;	// connected components

	size_t n = graph.num_nodes();
	while (graph.num_nodes()) {

		// nearest-neighbor chain
		size_t node = graph.next_node();
		if (node == (size_t)(-1)) {
			printf("node is -1\n");
			//exit(3);
		}

		std::vector< size_t> chain{node};
		while (chain.size()) {

			size_t a = chain.back();
			chain.pop_back();
			if (a == (size_t)(-1)) {
				printf("a is -1\n");
				//exit(3);
			}

			auto result = graph.nearest_neighbor(a);
			FloatType d = std::get<0>(result) / totalWeight / totalWeight;
			size_t b = std::get<1>(result);
			if (b == (size_t)(-1)) {
				if (chain.size() != 0) {
					printf("non-zero chain size!\n"); fflush(stdout);
					//exit(3);
				}
			}

			if (b == (size_t)(-1)) {
				// remove the connected component
				size_t sa = graph.snode[a];
				//components.push_back(std::make_tuple(graph.rep[a], sa));
				graph.remove_node(a);
			}
			else if (chain.size()) {
				size_t c = chain.back();
				chain.pop_back();

				if (b == c) {
					size_t size = graph.snode[a] + graph.snode[b];
					if (size == 0) {
						printf("zero size\n");
						//exit(3);
					}
					size_t parent = graph.contract_edge(a, b);

					double disorientation = INFINITY;
					if (parent == a) {
						disorientation = calculate_disorientation(structureType, qsum, a, b);
					}
					else {
						disorientation = calculate_disorientation(structureType, qsum, b, a);
					}

					//dendrogram.push_back(	DendrogramNode(	std::min(graph.rep[a], graph.rep[b]),
					//					std::max(graph.rep[a], graph.rep[b]),
					//					d, disorientation,
					//					size)

					auto node = DendrogramNode(	std::min(a, b),
									std::max(a, b),
									d, disorientation,
									size);
					*dendrogram++ = node;
				}
				else {
					chain.push_back(c);
					chain.push_back(a);
					chain.push_back(b);
					if (a == (size_t)(-1)) {printf("!a is -1\n"); /*exit(3);*/}
					if (b == (size_t)(-1)) {printf("!b is -1\n"); /*exit(3);*/}
					if (c == (size_t)(-1)) {printf("!c is -1\n"); /*exit(3);*/}
				}
			}
			else if (b != (size_t)(-1)) {
				chain.push_back(a);
				chain.push_back(b);
				if (a == (size_t)(-1)) {printf("#a is -1\n"); /*exit(3);*/}
				if (b == (size_t)(-1)) {printf("#b is -1\n"); /*exit(3);*/}
			}
		}
	}

	// add connected components to the dendrogram
	//size_t u = graph.next;
	//a, s = components.pop()
	//for b, t in components:
	//	s += t
	//	D.append([min(a, b), max(a, b), float("inf"), s])
	//	a = u
	//	u += 1

	return true;
}

}	// End of namespace
}	// End of namespace
}	// End of namespace
