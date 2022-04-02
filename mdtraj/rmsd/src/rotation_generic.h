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
#include "msvccompat.h"
#include "rotation.h"
#include "msvccompat.h"


void rot_atom_major(const int n_atoms, float* a, const float rot[9])
{
    /* Apply rotation matrix `rot` to conformation `a`.
    */
    unsigned int k;
    float x, y, z;

    for (k = 0; k < n_atoms; k++) {
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
    float ax, ay, az, bx, by, bz;
    float tx, ty, tz, dx, dy, dz, acculm;
    double sum_displacement = 0.0; /* accumulate the sum-squared-displacement here */

    for (i = 0; i < n_real_atoms; i++) {
        ax = a[3*i + 0];
        ay = a[3*i + 1];
        az = a[3*i + 2];
        bx = b[3*i + 0];
        by = b[3*i + 1];
        bz = b[3*i + 2];

        /* rotated coordinates of the 4 atoms */
        tx = ax*rot[0] + ay*rot[3] + az*rot[6];
        ty = ax*rot[1] + ay*rot[4] + az*rot[7];
        tz = ax*rot[2] + ay*rot[5] + az*rot[8];


        /* difference */
        dx = bx-tx;
        dy = by-ty;
        dz = bz-tz;

        /* sum up (bx-tx)^2 + (by-ty)^2 + (bz-tz)^2 over a block of 4 atoms */
        /* and accumulate the result into sum_displacement */
        acculm = dx*dx+dy*dy+dz*dz;
        sum_displacement += acculm;
    }
    return (float) (sum_displacement / (double) n_real_atoms);
}