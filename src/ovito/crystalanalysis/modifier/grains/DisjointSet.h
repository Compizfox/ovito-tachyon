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

namespace Ovito { namespace CrystalAnalysis {

class DisjointSet
{
public:

	DisjointSet(size_t n) {
		ranks.resize(n);
		parents.resize(n);
		sizes.resize(n);
		weights.resize(n);
		clear();
	}

	void clear() {
		std::iota(parents.begin(), parents.end(), (size_t)0);
		std::fill(sizes.begin(), sizes.end(), 1);
		std::fill(weights.begin(), weights.end(), 0);
	}

	// "Find" part of Union-Find.
	size_t find(size_t index) {

		// Find root and make root as parent of i (path compression)
		size_t parent = parents[index];
		while(parent != parents[parent]) {
			parent = parents[parent];
		}

		parents[index] = parent;
		return parent;
	}

	// "Union" part of Union-Find.
	size_t merge(size_t index1, size_t index2) {
		size_t parentA = find(index1);
		size_t parentB = find(index2);
		if(parentA == parentB) return index1;

		// Attach smaller rank tree under root of high rank tree (Union by Rank)
		if(ranks[parentA] < ranks[parentB]) {
			parents[parentA] = parentB;
			sizes[parentB] += sizes[parentA];
			weights[parentB] += weights[parentA];
			return index2;
		}
		else {
			parents[parentB] = parentA;
			sizes[parentA] += sizes[parentB];
			weights[parentA] += weights[parentB];

			// If ranks are same, then make one as root and increment its rank by one
			if(ranks[parentA] == ranks[parentB])
				ranks[parentA]++;

			return index1;
		}
	}

	std::vector<size_t> sizes;
	std::vector<FloatType> weights;

private:
	std::vector<size_t> parents;
	std::vector<size_t> ranks;
};


}	// End of namespace
}	// End of namespace
