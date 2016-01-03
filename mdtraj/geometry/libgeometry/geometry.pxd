cdef extern from "geometry.h":
    int compute_distances(
        const float* xyz, const ssize_t* pairs, const float* box_matrix,
        float* distance_out, float* displacement_out,
        const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_pairs) nogil
    int compute_angles(
        const float* xyz, const ssize_t* triplets,
        const float* box_matrix, float* out,
        const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_angles) nogil
    int compute_dihedrals(
        const float* xyz, const ssize_t* quartets,
        const float* box_matrix, float* out,
        const ssize_t n_frames, const ssize_t n_atoms, const ssize_t n_quartets) nogil
