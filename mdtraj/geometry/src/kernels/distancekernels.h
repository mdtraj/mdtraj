/**
 * Compute the distance/displacement  between pairs of atoms in every frame
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
 * pairs : array, shape=(n_pairs, 2)
 *     The specific pairs of atoms whose distance you want to compute. A 2d
 *     array of pairs, in C order.
 * box_matrix : array, shape=(n_frames, 3, 3)
 *     The box matrix for a each frame in the trajectory, in contiguous C
 * distance_out : array, shape=(n_frames, n_pairs), optional
 *     Array where the distances between pairs will be stored, in contiguous
 *     C order. If NULL is passed in, this return value will not be saved
 * displacement_out : array, shaoe=(n_frames, n_pairs, 3), optional
 *     An optional return value: if you'd also like to save the displacement
 *     vectors between the pairs, you can pass a pointer here. If
 *     displacement_out is NULL, then this variable will not be saved back
 *     to memory.
 *
 * All of the arrays are assumed to be contiguous. This code will
 * segfault if they're not.
 */

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
void dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
              float* distance_out, float* displacement_out,
              const int n_frames, const int n_atoms, const int n_pairs)

#else
void dist(const float* xyz, const int* pairs, float* distance_out,
         float* displacement_out, const int n_frames, const int n_atoms,
         const int n_pairs)
#endif
{
    bool store_displacement = (displacement_out != NULL);
    bool store_distance = (distance_out != NULL);
    for (int i = 0; i < n_frames; i++) {
        // Load the periodic box vectors.

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
        fvec4 box_size(box_matrix[0], box_matrix[4], box_matrix[8], 0);
        fvec4 inv_box_size(1.0f/box_matrix[0], 1.0f/box_matrix[4], 1.0f/box_matrix[8], 0);
#endif
        for (int j = 0; j < n_pairs; j++) {
            // Compute the displacement.
            int offset1 = 3*pairs[2*j + 0];
            fvec4 pos1(xyz[offset1], xyz[offset1+1], xyz[offset1+2], 0);
            int offset2 = 3*pairs[2*j + 1];
            fvec4 pos2(xyz[offset2], xyz[offset2+1], xyz[offset2+2], 0);
            fvec4 r12 = pos2-pos1;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
            r12 -= round(r12*inv_box_size)*box_size;
#endif

            // Store results.

            if (store_displacement) {
                float temp[4];
                r12.store(temp);
                *displacement_out = temp[0];
                displacement_out++;
                *displacement_out = temp[1];
                displacement_out++;
                *displacement_out = temp[2];
                displacement_out++;
            }
            if (store_distance) {
                *distance_out = sqrtf(dot3(r12, r12));
                distance_out++;
            }
        }

        // Advance to the next frame.

        xyz += n_atoms*3;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
        box_matrix += 9;
#endif
    }
}

/**
 * Compute the distance/displacement between pairs of atoms at t=0 and t=t for every frame
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
 * pairs : array, shape=(n_pairs, 2)
 *     The specific pairs of atoms whose distance you want to compute. A 2d
 *     array of pairs, in C order. The first atom in each pair will be considered at time
 * times : array, shape=(n_times, 2)
 *     The pairs of times whose distances are to be compared. A 2d array of
 *     times, in C order. This much match the `pairs` argument, i.e. for each
 *     pair of times and pairs, the comparison will bet made between the first
 *     atom at the first time and the second atom at the second time.
 * box_matrix : array, shape=(n_frames, 3, 3)
 *     The box matrix for a each frame in the trajectory, in contiguous C
 * distance_out : array, shape=(n_frames, n_pairs), optional
 *     Array where the distances between pairs will be stored, in contiguous
 *     C order. If NULL is passed in, this return value will not be saved
 * displacement_out : array, shaoe=(n_frames, n_pairs, 3), optional
 *     An optional return value: if you'd also like to save the displacement
 *     vectors between the pairs, you can pass a pointer here. If
 *     displacement_out is NULL, then this variable will not be saved back
 *     to memory.
 *
 * All of the arrays are assumed to be contiguous. This code will
 * segfault if they're not.
 */

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
void dist_mic_t(const float* xyz, const int* pairs, const int* times,
                const float* box_matrix, float* distance_out,
                float* displacement_out, const int n_times, const int n_atoms,
                const int n_pairs)

#else
void dist_t(const float* xyz, const int* pairs, const int* times,
            float* distance_out, float* displacement_out, const int n_times,
            const int n_atoms, const int n_pairs)
#endif
{
    bool store_displacement = (displacement_out != NULL);
    bool store_distance = (distance_out != NULL);
    for (int i = 0; i < n_times; i++) {
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
        // Based on first index of time pair
        int box_offset = times[2*i + 0] * 9;
        box_matrix += box_offset;
        fvec4 box_size(box_matrix[0], box_matrix[4], box_matrix[8], 0);
        fvec4 inv_box_size(1.0f/box_matrix[0], 1.0f/box_matrix[4], 1.0f/box_matrix[8], 0);
#endif

        for (int j = 0; j < n_pairs; j++) {
            // Find where in xyz to find the appropriate coordinates
            int time_offset1 = times[2*i + 0] * n_atoms * 3;
            int time_offset2 = times[2*i + 1] * n_atoms * 3;
            int pair_offset1 = pairs[2*j + 0] * 3;
            int pair_offset2 = pairs[2*j + 1] * 3;
            int offset1 = time_offset1 + pair_offset1;
            int offset2 = time_offset2 + pair_offset2;
            fvec4 pos1(xyz[offset1], xyz[offset1+1], xyz[offset1+2], 0);
            fvec4 pos2(xyz[offset2], xyz[offset2+1], xyz[offset2+2], 0);
            fvec4 r12 = pos2-pos1;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
            r12 -= round(r12*inv_box_size)*box_size;
#endif
    
            // Store results.
    
            if (store_displacement) {
                float temp[4];
                r12.store(temp);
                *displacement_out = temp[0];
                displacement_out++;
                *displacement_out = temp[1];
                displacement_out++;
                *displacement_out = temp[2];
                displacement_out++;
            }
            if (store_distance) {
                *distance_out = sqrtf(dot3(r12, r12));
                distance_out++;
            }
        }
// Reset box_matrix
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
        box_matrix -= box_offset;
#endif
    }
}
