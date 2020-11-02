/*Copyright (c) 2020 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#include <algorithm>
#include <cmath>
#include <cmath>
#include <cassert>
#include "ptm_constants.h"
#include "ptm_correspondences.h"


namespace ptm {

void index_to_permutation(int n, uint64_t encoded, int8_t* permutation)
{
    int base = n;
    uint64_t code[PTM_MAX_INPUT_POINTS] = {0};

    for (int i=0;i<base;i++) {
        code[i] = encoded % (base - i);
        encoded /= base - i;
    }

    for (int i=0;i<n;i++) {
        permutation[i] = i;
    }

    for (int i=0;i<base - 1;i++) {
        std::swap(permutation[i], permutation[i + code[i]]);
    }
}

uint64_t permutation_to_index(int n, int8_t* permutation)
{
	int p[PTM_MAX_INPUT_POINTS];
	int q[PTM_MAX_INPUT_POINTS];
    for (int i=0;i<n;i++) {
        p[i] = i;
        q[i] = i;
    }

    int8_t code[PTM_MAX_INPUT_POINTS] = {0};
    for (int i=0;i<n;i++) {
        int e = permutation[i];
        int d = q[e] - i;
        code[i] = d;
        if (d > 0) {
            int j = q[e];
            std::swap(q[p[i]], q[p[j]]);
            std::swap(p[i], p[j]);
        }
    }

    uint64_t encoded = 0;
    for (int i=0;i<n;i++) {
        uint64_t v = code[n - i - 1];
        encoded *= i + 1;
        encoded += v;
    }

    return encoded;
}

static bool is_single_shell(int type)
{
    if (type == PTM_MATCH_FCC
        || type == PTM_MATCH_HCP
        || type == PTM_MATCH_BCC
        || type == PTM_MATCH_ICO
        || type == PTM_MATCH_SC) {
        return true;
    }
    else if (type == PTM_MATCH_DCUB
             || type == PTM_MATCH_DHEX
             || type == PTM_MATCH_GRAPHENE) {
        return false;
    }
    else {
        assert(0);
        return false;
    }
}

static void vector_add(int n, int8_t* input, int8_t* transformed, int c)
{
    for (int i=0;i<n;i++) {
        transformed[i] = input[i] + c;
    }
}

void complete_correspondences(int n, int8_t* correspondences)
{
    bool hit[PTM_MAX_INPUT_POINTS] = {false};

    for (int i=0;i<n;i++) {
        int c = correspondences[i];
        hit[c] = true;
    }

    int c = n;
    for (int i=0;i<PTM_MAX_INPUT_POINTS;i++) {
        if (!hit[i]) {
            correspondences[c++] = i;
        }
    }
}

uint64_t encode_correspondences(int type, int8_t* correspondences)
{
    int num_nbrs = ptm_num_nbrs[type];
    int8_t transformed[PTM_MAX_INPUT_POINTS];

    if (is_single_shell(type)) {
        complete_correspondences(num_nbrs + 1, correspondences);
        vector_add(PTM_MAX_INPUT_POINTS - 1, &correspondences[1], transformed, -1);
        return permutation_to_index(PTM_MAX_INPUT_POINTS - 1, transformed);
    }
    else {
        int num_inner = 4, num_outer = 3;   //diamond types
        if (type == PTM_MATCH_GRAPHENE) {
            num_inner = 3;
            num_outer = 2;
        }

        for (int i=0;i<num_nbrs + 1;i++) {
            assert(correspondences[i] <= 17);
        }

        vector_add(num_inner, &correspondences[1], &transformed[0], -1);
        for (int i=0;i<num_inner;i++) {
            vector_add(num_outer,
                       &correspondences[1 + num_inner + i * num_outer],
                       &transformed[num_inner + i * num_outer], -1);
        }

        uint64_t encoded = 0;
        for (int i=0;i<num_nbrs;i++) {
            uint64_t t = transformed[i];
            encoded |= t << (4 * i);
        }

        return encoded;
    }
}

void decode_correspondences(int type, uint64_t encoded, int8_t* correspondences)
{
    int8_t decoded[PTM_MAX_INPUT_POINTS];

    if (is_single_shell(type)) {
        index_to_permutation(PTM_MAX_INPUT_POINTS - 1, encoded, decoded);

        correspondences[0] = 0;
        for (int i=0;i<PTM_MAX_INPUT_POINTS - 1;i++)
            correspondences[i + 1] = decoded[i] + 1;
    }
    else {
        int num_inner = 4, num_outer = 3;   //diamond types
        if (type == PTM_MATCH_GRAPHENE) {
            num_inner = 3;
            num_outer = 2;
        }

        int num_nbrs = ptm_num_nbrs[type];
        for (int i=0;i<num_nbrs;i++) {
            decoded[i] = encoded & 0xF;
            encoded >>= 4;
        }

        correspondences[0] = 0;
        vector_add(num_inner, &decoded[0], &correspondences[1], +1);
        for (int i=0;i<num_inner;i++) {
            vector_add(num_outer,
                       &decoded[num_inner + i * num_outer],
                       &correspondences[1 + num_inner + i * num_outer], +1);
        }

        complete_correspondences(num_nbrs + 1, correspondences);
    }
}

}

#ifdef __cplusplus
extern "C" {
#endif

uint64_t ptm_encode_correspondences(int type, int8_t* correspondences)
{
    return ptm::encode_correspondences(type, correspondences);
}

void ptm_decode_correspondences(int type, uint64_t encoded, int8_t* correspondences)
{
    return ptm::decode_correspondences(type, encoded, correspondences);
}

#ifdef __cplusplus
}
#endif

