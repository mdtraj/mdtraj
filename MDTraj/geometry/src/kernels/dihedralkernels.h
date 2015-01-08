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
int dihedral_mic(const float* xyz, const int* quartets,
              const float* box_matrix, float* out,
              const int n_frames, const int n_atoms, const int n_quartets)
#else
int dihedral(const float* xyz, const int* quartets, float* out,
          const int n_frames, const int n_atoms, const int n_quartets)
#endif
{
  int i, j;
  __m128 x0, x1, x2, x3, b1, b2, b3, c1, c2, p1, p2;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
  __m128 hinv[3];
  __m128 h[3];
#endif

  for (i = 0; i < n_frames; i++) {
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    loadBoxMatrix(box_matrix, &h, &hinv);
#endif

    for (j = 0; j < n_quartets; j++) {
      x0 = load_float3(xyz + 3*quartets[4*j + 0]);
      x1 = load_float3(xyz + 3*quartets[4*j + 1]);
      x2 = load_float3(xyz + 3*quartets[4*j + 2]);
      x3 = load_float3(xyz + 3*quartets[4*j + 3]);

      b1 = _mm_sub_ps(x1, x0);
      b2 = _mm_sub_ps(x2, x1);
      b3 = _mm_sub_ps(x3, x2);

#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
      b1 = minimum_image(b1, &h, &hinv);
      b2 = minimum_image(b2, &h, &hinv);
      b3 = minimum_image(b3, &h, &hinv);
#endif

      c1 = cross(b2, b3);
      c2 = cross(b1, b2);

      p1 = _mm_mul_ps(_mm_dp_ps2(b1, c1, 0x71), _mm_sqrt_ps(_mm_dp_ps2(b2, b2, 0x71)));
      p2 = _mm_dp_ps2(c1, c2, 0x71);

      *(out++) = (float) atan2(_mm_cvtss_f32(p1), _mm_cvtss_f32(p2));
    };
    xyz += n_atoms*3;
#ifdef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
    box_matrix += 9;
#endif
  }
  return 1;
}
