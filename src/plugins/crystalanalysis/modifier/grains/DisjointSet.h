#ifndef DISJOINTSET_H
#define DISJOINTSET_H

#include <core/utilities/FloatType.h>


namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

class DisjointSet
{
public:
	void clear() {
		std::iota(parents.begin(), parents.end(), (size_t)0);
		std::fill(sizes.begin(), sizes.end(), 1);
		std::fill(weights.begin(), weights.end(), 0);
	}

	DisjointSet(size_t n)
	{
		ranks.resize(n);
		parents.resize(n);
		sizes.resize(n);
		weights.resize(n);

		clear();
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
}	// End of namespace

#endif

