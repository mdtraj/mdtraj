/*////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////*/

#include "stdio.h"
#include <assert.h>
#include "msvccompat.h"
#include "util_arm.h"
#include "rotation.h"
#include "msvccompat.h"


INLINE float32x4_t add3(float32x4_t a, float32x4_t b, float32x4_t c) {
    return vaddq_f32(vaddq_f32(a, b), c);
}


void rot_atom_major(const int n_atoms, float* a, const float rot[9])
{
    /* Apply rotation matrix `rot` to conformation `a`.
    */

    unsigned int k;
    unsigned int n_iters = 0;
    float x, y, z;

    float32x4x3_t axyz;
    float32x4_t tx, ty, tz;
    float32x4_t rXX = vld1q_dup_f32(rot + 0);
    float32x4_t rXY = vld1q_dup_f32(rot + 1);
    float32x4_t rXZ = vld1q_dup_f32(rot + 2);
    float32x4_t rYX = vld1q_dup_f32(rot + 3);
    float32x4_t rYY = vld1q_dup_f32(rot + 4);
    float32x4_t rYZ = vld1q_dup_f32(rot + 5);
    float32x4_t rZX = vld1q_dup_f32(rot + 6);
    float32x4_t rZY = vld1q_dup_f32(rot + 7);
    float32x4_t rZZ = vld1q_dup_f32(rot + 8);

    n_iters = n_atoms / 4;

    for (k = 0; k < n_iters; k++) {
        /* load four atoms at a time */
        axyz = aos_deinterleaved_load(a);

        tx = add3(vmulq_f32(axyz.val[0], rXX), vmulq_f32(axyz.val[1], rYX), vmulq_f32(axyz.val[2], rZX));
        ty = add3(vmulq_f32(axyz.val[0], rXY), vmulq_f32(axyz.val[1], rYY), vmulq_f32(axyz.val[2], rZY));
        tz = add3(vmulq_f32(axyz.val[0], rXZ), vmulq_f32(axyz.val[1], rYZ), vmulq_f32(axyz.val[2], rZZ));

        aos_interleaved_store(a, tx, ty, tz);

        a += 12;
    }

    /* Epilogue to process the last atoms that are past the last multiple of */
    /* four */
    for (k = 0; k < n_atoms % 4; k++) {
        x = a[3*k + 0];
        y = a[3*k + 1];
        z = a[3*k + 2];
        a[3*k + 0] = x*rot[0] + y*rot[3] + z*rot[6];
        a[3*k + 1] = x*rot[1] + y*rot[4] + z*rot[7];
        a[3*k + 2] = x*rot[2] + y*rot[5] + z*rot[8];
    }
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

    float32x4x3_t axyz, bxyz;    
    float32x4_t tx, ty, tz, dx, dy, dz, acculm;
    float32x4_t rXX = vld1q_dup_f32(rot + 0);
    float32x4_t rXY = vld1q_dup_f32(rot + 1);
    float32x4_t rXZ = vld1q_dup_f32(rot + 2);
    float32x4_t rYX = vld1q_dup_f32(rot + 3);
    float32x4_t rYY = vld1q_dup_f32(rot + 4);
    float32x4_t rYZ = vld1q_dup_f32(rot + 5);
    float32x4_t rZX = vld1q_dup_f32(rot + 6);
    float32x4_t rZY = vld1q_dup_f32(rot + 7);
    float32x4_t rZZ = vld1q_dup_f32(rot + 8);
    n_iters = (n_padded_atoms >> 2);
    assert(n_padded_atoms % 4 == 0);

    for (i = 0; i < n_iters; i++) {
        /* load four atoms at a time */
        axyz = aos_deinterleaved_load(a);
        bxyz = aos_deinterleaved_load(b);

        /* rotated coordinates of the 4 atoms */
        tx = add3(vmulq_f32(axyz.val[0], rXX), vmulq_f32(axyz.val[1], rYX), vmulq_f32(axyz.val[2], rZX));
        ty = add3(vmulq_f32(axyz.val[0], rXY), vmulq_f32(axyz.val[1], rYY), vmulq_f32(axyz.val[2], rZY));
        tz = add3(vmulq_f32(axyz.val[0], rXZ), vmulq_f32(axyz.val[1], rYZ), vmulq_f32(axyz.val[2], rZZ));

        /* difference */
        dx = vsubq_f32(bxyz.val[0], tx);
        dy = vsubq_f32(bxyz.val[1], ty);
        dz = vsubq_f32(bxyz.val[2], tz);

        /* sum up (bx-tx)^2 + (by-ty)^2 + (bz-tz)^2 over a block of 4 atoms */
        /* and accumulate the result into sum_displacement */
        acculm = add3(vmulq_f32(dx, dx), vmulq_f32(dy, dy), vmulq_f32(dz, dz));
        /* horizontal add of all 4 elemenets in acculm */    
        float32x4_t tmp = vpaddq_f32(acculm,acculm); /* tmp = a01 a23 a01 a23 */
        sum4 = vgetq_lane_f32(tmp, 0) + vgetq_lane_f32(tmp, 1);
        sum_displacement += sum4;
        a += 12;
        b += 12;
    }
    return (float) (sum_displacement / (double) n_real_atoms);
}
