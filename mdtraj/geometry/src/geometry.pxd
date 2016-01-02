cdef extern from "geometry.h":
    int compute_distances(
        const float* xyz, const int* pairs, float* distance_out,
        float* displacement_out, const int n_frames, const int n_atoms,
        const int n_pairs) nogil
    int compute_distances_orthorhombic(
        const float* xyz, const int* pairs, const float* box_matrix,
        float* distance_out, float* displacement_out,
        const int n_frames, const int n_atoms, const int n_pairs) nogil
    int compute_distances_triclinic(
        const float* xyz, const int* pairs, const float* box_matrix,
        float* distance_out, float* displacement_out,
        const int n_frames, const int n_atoms, const int n_pairs) nogil

    int compute_angles(
        const float* xyz, const int* triplets, float* out,
        const int n_frames, const int n_atoms, const int n_angles) nogil
    int compute_angles_orthorhombic(
        const float* xyz, const int* triplets,
        const float* box_matrix, float* out,
        const int n_frames, const int n_atoms, const int n_angles) nogil

    int compute_dihedrals(
        const float* xyz, const int* quartets, float* out,
        const int n_frames, const int n_atoms, const int n_quartets) nogil
    int compute_dihedrals_orthorhombic(
        const float* xyz, const int* quartets,
        const float* box_matrix, float* out,
        const int n_frames, const int n_atoms, const int n_quartets) nogil
