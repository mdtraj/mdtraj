#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <cmath>
#include "vectorize_sse.h"

#ifdef __cplusplus
extern "C" {
#endif

int compute_distances(
    const float* xyz, const ssize_t* pairs, const float* box_matrix,
    float* distance_out, float* displacement_out,
    const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_pairs);

int compute_angles(
    const float* xyz, const ssize_t* triplets,
    const float* box_matrix, float* out,
    const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_angles);

int compute_dihedrals(
    const float* xyz, const ssize_t* quartets,
    const float* box_matrix, float* out,
    const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_quartets);


#ifdef __cplusplus
}
#endif
#endif
