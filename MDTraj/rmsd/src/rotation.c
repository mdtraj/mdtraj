/////////////////////////////////////////////////////////////////////////////
// MDTraj: A Python Library for Loading, Saving, and Manipulating
//         Molecular Dynamics Trajectories.
// Copyright 2012-2013 Stanford University and the Authors
//
// Authors: Robert McGibbon
// Contributors:
//
// MDTraj is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////

#include "stdio.h"
#include <assert.h>
#include "util.h"
#include "theobald_rmsd.h"
#include "rotation.h"
#ifndef INLINE
#if defined(__GNUC__)
#define INLINE __inline__
#elif defined(_MSC_VER)
#define INLINE __inline
#elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define INLINE inline
#else
#define INLINE
#endif
#endif

static INLINE __m128 _mm_add3_ps(__m128 a, __m128 b, __m128 c) {
    return _mm_add_ps(_mm_add_ps(a, b), c);
}


float rot_msd_atom_major(const int n_real_atoms, const int n_padded_atoms,
                         const float* a, const float* b, const float rot[9])
/* Apply rotation matrix `rot` to conformation `a`, and returns its
 * mean-squared-displacement with respect to conformation `b`.
 *
 * Note: `rot` need *NOT* be the optimal rotation matrix for aligning the two
 * structures, so this function does not compute the distance after an optimal
 * alignment.
 *
 * The data layout for `a` and `b` follow the same pattern as used by
 * the function msd_atom_major, in theobald_rmsd.c
 *
 * the layout in memory for a structure of 7 atoms would look like
 * this (read row-major):
 *
 *   x0 y0 z0
 *   x1 y1 z1
 *   x2 y2 x2
 *   x3 y3 x3
 *   x4 y4 x4
 *   x5 y5 x5
 *   x6 y6 x6
 *    0  0  0
 */
{
    unsigned int i;
    unsigned int n_iters = 0;
    double sum_displacement = 0.0; /* accumulate the sum-squared-displacement here */
    float sum4 = 0; /* a single sum-squared-displacement for 4 values */

    __m128 ax, ay, az, bx, by, bz, tx, ty, tz, dx, dy, dz, acculm;
    __m128 rXX = _mm_load1_ps(rot + 0);
    __m128 rXY = _mm_load1_ps(rot + 1);
    __m128 rXZ = _mm_load1_ps(rot + 2);
    __m128 rYX = _mm_load1_ps(rot + 3);
    __m128 rYY = _mm_load1_ps(rot + 4);
    __m128 rYZ = _mm_load1_ps(rot + 5);
    __m128 rZX = _mm_load1_ps(rot + 6);
    __m128 rZY = _mm_load1_ps(rot + 7);
    __m128 rZZ = _mm_load1_ps(rot + 8);
    n_iters = (n_padded_atoms >> 2);
    assert(n_padded_atoms % 4 == 0);

    for (i = 0; i < n_iters; i++) {
        /* load four atoms at a time */
        aos_deinterleaved_load(a,&ax,&ay,&az);
        aos_deinterleaved_load(b,&bx,&by,&bz);

        /* rotated coordinates of the 4 atoms */
        tx = _mm_add3_ps(_mm_mul_ps(ax, rXX), _mm_mul_ps(ay, rYX), _mm_mul_ps(az, rZX));
        ty = _mm_add3_ps(_mm_mul_ps(ax, rXY), _mm_mul_ps(ay, rYY), _mm_mul_ps(az, rZY));
        tz = _mm_add3_ps(_mm_mul_ps(ax, rXZ), _mm_mul_ps(ay, rYZ), _mm_mul_ps(az, rZZ));

        /* difference */
        dx = _mm_sub_ps(bx, tx);
        dy = _mm_sub_ps(by, ty);
        dz = _mm_sub_ps(bz, tz);

        /* sum up (bx-tx)^2 + (by-ty)^2 + (bz-tz)^2 over a block of 4 atoms */
        /* and accumulate the result into sum_displacement */
        acculm = _mm_add3_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy), _mm_mul_ps(dz, dz));
        /* horizontal add of all 4 elemenets in acculm */
        acculm = _mm_add_ps(acculm, _mm_movehl_ps(acculm, acculm));
        acculm = _mm_add_ss(acculm, _mm_shuffle_ps(acculm, acculm, 1));
        _mm_store_ss(&sum4, acculm);
        sum_displacement += sum4;
        a += 12;
        b += 12;
    }
    return (float) (sum_displacement / (double) n_real_atoms);
}
