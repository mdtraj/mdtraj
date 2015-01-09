#include "stdio.h"
#include <vector>
#include <pmmintrin.h>
#include "ssetools.h"
#include "msvccompat.h"
#include "geometryutils.h"
#include "neighbors.hpp"


/**
 * Compute the distance from atom `i` to atom `j`, whose coordinates
 * are the `i`th and `j`th row of the 2D array of cartesian coordinates
 * `frame_xyz` (optionally using periodic boundary conditions).
 */
template<bool periodic> float get_dist(float* frame_xyz, int i, int j,
                                       float* box_matrix)
{
    float result;
    __m128 x1, x2, r12, r12_2, s;

    x1 = load_float3(frame_xyz + 3*i);
    x2 = load_float3(frame_xyz + 3*j);
    /* r12 = x2 - x1 */
    r12 = _mm_sub_ps(x2, x1);

    if (periodic) {
        __m128 hinv[3];
        __m128 h[3];
        loadBoxMatrix(box_matrix, &h, &hinv);
        r12 = minimum_image(r12, &h, &hinv);
    }

    /* r12_2 = r12*r12 */
    r12_2 = _mm_mul_ps(r12, r12);
    /* horizontal add the components of d2 (last one is zero)*/
    s = _mm_hsum_ps(r12_2);
    /* sqrt our final answer */
    s = _mm_sqrt_ps(s);

    _mm_store_ss(&result, s);
    return result;
}



/**
 *
 *
 */
std::vector<int> _compute_neighbors(
    float* frame_xyz, int n_atoms, float cutoff,
    const std::vector<int>& query_indices,
    const std::vector<int>& haystack_indices,
    float* box_matrix)
{
    std::vector<int> result;

    std::vector<int>::const_iterator hit;
    for (hit = haystack_indices.begin(); hit != haystack_indices.end(); ++hit) {
        // is this haystack atom within cutoff of _any_ query atom?
        bool match = false;

        std::vector<int>::const_iterator qit;
        for (qit = query_indices.begin(); qit != query_indices.end(); ++qit) {
            // compute distance from haystack atom *hit to query atom *qit
            if (*hit == *qit) {
                continue;
            }
            float dist = 0;
            if (box_matrix == NULL) {
                dist = get_dist<false>(frame_xyz, *hit, *qit, box_matrix);
            } else {
                dist = get_dist<true>(frame_xyz, *hit, *qit, box_matrix);
            }

            if (dist < cutoff) {
                 match = true;
                 break;
             }
         }

         // this haystack atom is within cutoff of at least 1 query atom
         if (match) {
             result.push_back(*hit);
         }
    }

    return result;
}
