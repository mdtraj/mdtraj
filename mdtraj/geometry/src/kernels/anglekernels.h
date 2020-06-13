/**
 *  Compute the angle between triples of atoms in every frame of
 *  xyz.
 *
 * Two versions of this function can be compiled, one of which takes an extra
 * `box_matrix` argument that uses the minimum image convention and periodic
 * boundary conditions.
 *
 *  Parameters
 *  ----------
 *  xyz : array, shape=(n_frames, n_atoms, 3)
 *      Cartesian coordinates of the atoms in every frame, in contiguous C order.
 *  triplets : array, shape=(n_angles, 3)
 *      The specific tripple of atoms whose angle you want to compute. The
 *      angle computed will be centered around the middle element (i.e aABC).
 *      A 2d array of indices, in C order.
 *  box_matrix : array, shape=(n_frames, 3, 3)
 *      The box matrix for a each frame in the trajectory, in contigous C order.
 *  out : array, shape=(n_frames, n_pairs)
 *      Array where the angles will be stored, in contiguous C order.
 *
 *  All of the arrays are assumed to be contiguous. This code will
 *  segfault if they're not.
 */
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#ifdef COMPILE_WITH_TRICLINIC
void angle_mic_triclinic(const float* xyz, const int* triplets,
               const float* box_matrix, float* out,
               const int n_frames, const int n_atoms, const int n_angles)
#else
void angle_mic(const float* xyz, const int* triplets,
               const float* box_matrix, float* out,
               const int n_frames, const int n_atoms, const int n_angles)
#endif
#else
void angle(const float* xyz, const int* triplets, float* out,
           const int n_frames, const int n_atoms, const int n_angles)
#endif
{
    std::vector<float> distances(2*n_frames);
    std::vector<float> displacements(6*n_frames);
    for (int i = 0; i < n_angles; i++) {
        int pairs[4] = {triplets[3*i+1], triplets[3*i], triplets[3*i+1], triplets[3*i+2]};
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#ifdef COMPILE_WITH_TRICLINIC
        dist_mic_triclinic(xyz, pairs, box_matrix, &distances[0], &displacements[0], n_frames, n_atoms, 2);
#else
        dist_mic(xyz, pairs, box_matrix, &distances[0], &displacements[0], n_frames, n_atoms, 2);
#endif
#else
        dist(xyz, pairs, &distances[0], &displacements[0], n_frames, n_atoms, 2);
#endif
        for (int j = 0; j < n_frames; j++) {
            fvec4 v1(displacements[6*j], displacements[6*j+1], displacements[6*j+2], 0);
            fvec4 v2(displacements[6*j+3], displacements[6*j+4], displacements[6*j+5], 0);
            float angle = (float) acos(dot3(v1, v2)/(distances[2*j]*distances[2*j+1]));
            out[n_angles*j + i] = angle;
        }
    }
}
