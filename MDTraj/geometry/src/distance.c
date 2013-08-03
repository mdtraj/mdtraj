#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>   // (Meta-header, for GCC only)

inline __m128 LoadFloat3(float* value) {
  // Load (x,y,z) into a SSE register, leaving the last entry
  // set to zero.
  __m128 x = _mm_load_ss(&value[0]);
  __m128 y = _mm_load_ss(&value[1]);
  __m128 z = _mm_load_ss(&value[2]);
  __m128 xy = _mm_movelh_ps(x, y);
  return _mm_shuffle_ps(xy, z, _MM_SHUFFLE(2, 0, 2, 0));
}

int dist(float* xyz, int* pairs, float* out,
	 int n_frames, int n_atoms, int n_pairs) {
  /* Compute the distance between pairs of atoms in every frame
     of xyz.

     Uses SSE3 horizontal add for the reduction (compile with -msse3)

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     pairs : array, shape=(n_pairs, 2)
         The specific pairs of atoms whose distance you want to compute. A 2d
         array of pairs, in C order.
     out : array, shape=(n_frames, n_pairs)
         Array where the output will be stored, in contiguous C order.

     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
   */

  int i, j;
  __m128 x1, x2, d, d2, s;

  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_pairs; j++) {
      // Load the two vectors whos distance we want to compute
      // x1 = xyz[i, pairs[j,0], 0:3]
      // x2 = xyz[i, pairs[j,1], 0:3]
      x1 = LoadFloat3(xyz + 3*pairs[2*j + 0]);
      x2 = LoadFloat3(xyz + 3*pairs[2*j + 1]);

      // d = x2 - x1
      d = _mm_sub_ps(x2, x1);
      // d2 = d*d
      d2 = _mm_mul_ps(d, d);

      // horizontal add the components of d2 with
      // two instructions. note: it's critical
      // here that the last entry of x1 and x2 was 0
      // so that d2.w = 0
      s = _mm_hadd_ps(d2, d2);
      s = _mm_hadd_ps(s, s);
      // sqrt our final answer
      s = _mm_sqrt_ps(s);
      
      // s now contains our answer in all four elements, because
      // of the way the hadd works. we only want to store one
      // element.
      _mm_store_ss(out++, s);
    }

    // advance to the next frame
    xyz += n_atoms*3;
  }

  return 1;
}

int dist_orthorhombic(float* xyz, int* atom_pairs, float* out) {
  return 1;
}
