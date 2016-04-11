/**
 * Compute the angle between sets of four atoms in every frame
 * of xyz.
 *
 * Two versions of this function can be compiled, one of which takes an extra
 * `box_matrix` argument that uses the minimum image convention and periodic
 * boundary conditions.
 *
 * Parameters
 * ----------
 * xyz : array, shape=(n_frames, n_atoms, 3)
 *     Cartesian coordinates of the atoms in every frame, in contiguous C order.
 * quartets : array, shape=(n_quartets, 3)
 *     The specific quartet of atoms whose angle you want to compute. The
 *     angle computed will be the torsion around the bound between the
 *     middle two elements (i.e aABCD). A 2d array of indices, in C order.
 * box_matrix : array, shape=(n_frames, 3, 3)
 *     The box matrix for a each frame in the trajectory, in contiguous C order.
 * out : array, shape=(n_frames, n_pairs)
 *     Array where the angles will be stored, in contiguous C order.
 *
 * All of the arrays are assumed to be contiguous. This code will
 * segfault if they're not.
*/

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#ifdef COMPILE_WITH_TRICLINIC
void dihedral_mic_triclinic(const float* xyz, const int* quartets,
                            const float* box_matrix, float* out,
                            const int n_frames, const int n_atoms, const int n_quartets)
#else
void dihedral_mic(const float* xyz, const int* quartets,
                  const float* box_matrix, float* out,
                  const int n_frames, const int n_atoms, const int n_quartets)
#endif
#else
void dihedral(const float* xyz, const int* quartets, float* out,
              const int n_frames, const int n_atoms, const int n_quartets)
#endif
{
    std::vector<float> distances(3*n_frames);
    std::vector<float> displacements(9*n_frames);
    for (int i = 0; i < n_quartets; i++) {
        int pairs[6] = {quartets[4*i], quartets[4*i+1], quartets[4*i+1], quartets[4*i+2], quartets[4*i+2], quartets[4*i+3]};
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#ifdef COMPILE_WITH_TRICLINIC
        dist_mic_triclinic(xyz, pairs, box_matrix, &distances[0], &displacements[0], n_frames, n_atoms, 3);
#else
        dist_mic(xyz, pairs, box_matrix, &distances[0], &displacements[0], n_frames, n_atoms, 3);
#endif
#else
        dist(xyz, pairs, &distances[0], &displacements[0], n_frames, n_atoms, 3);
#endif
        for (int j = 0; j < n_frames; j++) {
            fvec4 v1(displacements[9*j], displacements[9*j+1], displacements[9*j+2], 0);
            fvec4 v2(displacements[9*j+3], displacements[9*j+4], displacements[9*j+5], 0);
            fvec4 v3(displacements[9*j+6], displacements[9*j+7], displacements[9*j+8], 0);
            fvec4 c1 = cross(v2, v3);
            fvec4 c2 = cross(v1, v2);
            float p1 = dot3(v1, c1)*distances[3*j+1];
            float p2 = dot3(c1, c2);
            out[n_quartets*j + i] = atan2f(p1, p2);
        }
    }
}
