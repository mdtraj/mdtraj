#ifndef MDTRAJ_NEIGHBORS_H
#define MDTRAJ_NEIGHBORS_H

std::vector<int> _compute_neighbors(
    float* frame_xyz, int n_atoms, float cutoff,
    const std::vector<int>& query_indices,
    const std::vector<int>& haystack_indices,
    float* box_matrix);

#endif