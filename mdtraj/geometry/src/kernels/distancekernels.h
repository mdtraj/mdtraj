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
int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
             float* distance_out, float* displacement_out,
             const int n_frames, const int n_atoms, const int n_pairs)
#else
int dist(const float* xyz, const int* pairs, float* distance_out,
         float* displacement_out, const int n_frames, const int n_atoms,
         const int n_pairs)
#endif
{
  int i, j;
  int store_displacement = displacement_out == NULL ? 0 : 1;
  int store_distance = distance_out == NULL ? 0 : 1;
  __m128 x1, x2, r12, r12_2, s;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
  __m128 hinv[3];
  __m128 h[3];
#endif

  for (i = 0; i < n_frames; i++) {
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    loadBoxMatrix(box_matrix, &h, &hinv);
#endif

    for (j = 0; j < n_pairs; j++) {
      /* Load the two vectors whos distance we want to compute */
      /* x1 = xyz[i, pairs[j,0], 0:3] */
      /* x2 = xyz[i, pairs[j,1], 0:3] */
      x1 = load_float3(xyz + 3*pairs[2*j + 0]);
      x2 = load_float3(xyz + 3*pairs[2*j + 1]);

      /* r12 = x2 - x1 */
      r12 = _mm_sub_ps(x2, x1);

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
      r12 = minimum_image(r12, &h, &hinv);
#endif

      if (store_displacement) {
        /* store the two lower entries (x,y) in memory */
        _mm_storel_pi((__m64*)(displacement_out), r12);
        displacement_out += 2;
        /* swap high-low and then store the z entry in the memory */
        _mm_store_ss(displacement_out++, _mm_movehl_ps(r12, r12));
      }
      if (store_distance) {
        /* r12_2 = r12*r12 */
        r12_2 = _mm_mul_ps(r12, r12);
        /* horizontal add the components of d2 (last one is zero)*/
        s = _mm_hsum_ps(r12_2);
        /* sqrt our final answer */
        s = _mm_sqrt_ps(s);
        /* s now contains our answer in all four elements, because */
        /* of the way the hadd works. we only want to store one */
        /* element. */
        _mm_store_ss(distance_out++, s);
      }
    }

    /* advance to the next frame */
    xyz += n_atoms*3;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
     box_matrix += 9;
#endif
  }

  return 1;
}

