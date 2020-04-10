/*Copyright (c) 2020 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#include <cmath>
#include "ptm_constants.h"
#include "ptm_correspondences.h"


namespace ptm {

// taken from http://antoinecomeau.blogspot.com/2014/07/mapping-between-permutations-and.html
void index_to_permutation(int n, uint64_t k, int8_t* permuted)
{
	int elems[PTM_MAX_INPUT_POINTS];
	for(int i=0;i<n;i++)
		elems[i] = i;

	uint64_t m = k;
	for(int i=0;i<n;i++)
	{
		uint64_t ind = m % (n - i);
		m = m / (n - i);
		permuted[i] = elems[ind];
		elems[ind] = elems[n - i - 1];
	}
}

// taken from http://antoinecomeau.blogspot.com/2014/07/mapping-between-permutations-and.html
uint64_t permutation_to_index(int n, int8_t* permutation)
{
	int pos[PTM_MAX_INPUT_POINTS];
	int elems[PTM_MAX_INPUT_POINTS];
	for(int i=0;i<n;i++)
	{
		pos[i] = i;
		elems[i] = i;
	}

	uint64_t m = 1;
	uint64_t k = 0;
	for(int i=0;i<n-1;i++)
	{
		k += m * pos[permutation[i]];
		m = m * (n - i);
		pos[elems[n - i - 1]] = pos[permutation[i]];
		elems[pos[permutation[i]]] = elems[n - i - 1];
	}

	return k;
}

}

#ifdef __cplusplus
extern "C" {
#endif

void ptm_index_to_permutation(int n, uint64_t k, int8_t* permuted)
{
	return ptm::index_to_permutation(n, k, permuted);
}

uint64_t ptm_permutation_to_index(int n, int8_t* permutation)
{
	return ptm::permutation_to_index(n, permutation);
}

#ifdef __cplusplus
}
#endif

