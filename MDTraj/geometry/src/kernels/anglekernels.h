
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
int angle_mic(const float* xyz, const int* triplets,
            const float* box_matrix, float* out,
            const int n_frames, const int n_atoms, const int n_angles)
#else
int angle(const float* xyz, const int* triplets, float* out,
          const int n_frames, const int n_atoms, const int n_angles)
#endif
{
  int i, j;
  __m128 r_m, r_n, r_o, u_prime, u, v_prime, v;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
  __m128 hinv[3];
  __m128 h[3];
#endif

  for (i = 0; i < n_frames; i++) {
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    loadBoxMatrix(box_matrix, &h, &hinv);
#endif

    for (j = 0; j < n_angles; j++) {
      r_m = load_float3(xyz + 3*triplets[3*j + 0]);
      r_o = load_float3(xyz + 3*triplets[3*j + 1]);
      r_n = load_float3(xyz + 3*triplets[3*j + 2]);

      u_prime = _mm_sub_ps(r_m, r_o);
      v_prime = _mm_sub_ps(r_n, r_o);

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
      u_prime = minimum_image(u_prime, &h, &hinv);
      v_prime = minimum_image(v_prime, &h, &hinv);
#endif

      /* normalize the vectors u_prime and v_prime */
      u = _mm_div_ps(u_prime, _mm_sqrt_ps(_mm_dp_ps2(u_prime, u_prime, 0x7F)));
      v = _mm_div_ps(v_prime, _mm_sqrt_ps(_mm_dp_ps2(v_prime, v_prime, 0x7F)));

      /* compute the arccos of the dot product, and store the result. */
      *(out++) = (float) acos(CLIP(_mm_cvtss_f32(_mm_dp_ps2(u, v, 0x71)), -1, 1));
    }
    /* advance to the next frame */
    xyz += n_atoms*3;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    box_matrix += 9;
#endif
  }

  return 1;
}
