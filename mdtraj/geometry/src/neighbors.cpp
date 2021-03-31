#include <cmath>
#include <vector>
#include "msvccompat.h"
#include "vectorize.h"
#include "neighbors.hpp"
#include "math_patch.h"

std::vector<int> _compute_neighbors(
    float* frame_xyz, int n_atoms, float cutoff,
    const std::vector<int>& query_indices,
    const std::vector<int>& haystack_indices,
    float* box_matrix)
{
    float cutoff2 = cutoff*cutoff;
    bool periodic = (box_matrix != NULL);
    bool triclinic = periodic && (box_matrix[1] != 0 || box_matrix[2] != 0 ||
            box_matrix[3] != 0 || box_matrix[5] != 0 || box_matrix[6] != 0 || box_matrix[7] != 0);
    float recip_box_size[3] = {0.0f, 0.0f, 0.0f};
    fvec4 box_size, inv_box_size, box_vec1, box_vec2, box_vec3;
    if (periodic) {
        recip_box_size[0] = 1.0f/box_matrix[0];
        recip_box_size[1] = 1.0f/box_matrix[4];
        recip_box_size[2] = 1.0f/box_matrix[8];
        box_size = fvec4(box_matrix[0], box_matrix[4], box_matrix[8], 0);
        inv_box_size = fvec4(recip_box_size[0], recip_box_size[1], recip_box_size[2], 0);
        box_vec1 = fvec4(box_matrix[0], box_matrix[1], box_matrix[2], 0);
        box_vec2 = fvec4(box_matrix[3], box_matrix[4], box_matrix[5], 0);
        box_vec3 = fvec4(box_matrix[6], box_matrix[7], box_matrix[8], 0);
        box_vec3 -= box_vec2*roundf(box_vec3[1]/box_vec2[1]);
        box_vec3 -= box_vec1*roundf(box_vec3[0]/box_vec1[0]);
        box_vec2 -= box_vec1*roundf(box_vec2[0]/box_vec1[0]);
    }
    std::vector<int> result;
    std::vector<int>::const_iterator hit;
    for (hit = haystack_indices.begin(); hit != haystack_indices.end(); ++hit) {
        // Is this haystack atom within cutoff of _any_ query atom?
        
        int i = *hit;
        fvec4 pos1(frame_xyz[3*i], frame_xyz[3*i+1], frame_xyz[3*i+2], 0);
        std::vector<int>::const_iterator qit;
        for (qit = query_indices.begin(); qit != query_indices.end(); ++qit) {
            // Compute distance from haystack atom *hit to query atom *qit
            
            int j = *qit;
            if (i == j)
                continue;
            fvec4 pos2(frame_xyz[3*j], frame_xyz[3*j+1], frame_xyz[3*j+2], 0);
            fvec4 delta = pos1-pos2;
            if (triclinic) {
                delta -= box_vec3*roundf(delta[2]*recip_box_size[2]);
                delta -= box_vec2*roundf(delta[1]*recip_box_size[1]);
                delta -= box_vec1*roundf(delta[0]*recip_box_size[0]);
            }
            else if (periodic)
                delta -= round(delta*inv_box_size)*box_size;
            float dist2 = dot3(delta, delta);
            if (dist2 < cutoff2) {
                // The haystack atom is within cutoff of this query atom.
                
                result.push_back(i);
                break;
             }
        }
    }
    return result;
}
