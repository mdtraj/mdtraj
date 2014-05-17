/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2014- Stanford University and the Authors                   */
/*                                                                       */
/* Authors: Robert McGibbon                                              */
/* Contributors:                                                         */
/*                                                                       */
/* MDTraj is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as        */
/* published by the Free Software Foundation, either version 2.1         */
/* of the License, or (at your option) any later version.                */
/*                                                                       */
/* This library is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU Lesser General Public License for more details.                   */
/*                                                                       */
/* You should have received a copy of the GNU Lesser General Public      */
/* License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.*/
/*=======================================================================*/

#include "math.h"
#include "stdlib.h"
#include "dridkernels.h"
#include "ssetools.h"
#include "msvccompat.h"
#include "moments.h"
#include <pmmintrin.h>
#include "cbrt.h"

int
drid_moments(float* coords, int index, int* partners, int n_partners, double* moments)
{
    int i;
    float d;
    moments_t onlinemoments;
    __m128 x, y, r, r2, s;
    
    moments_clear(&onlinemoments);
    x = load_float3(&coords[3 * index]);

    for (i = 0; i < n_partners; i++) {
        y = load_float3(&coords[3 * partners[i]]);
        r = _mm_sub_ps(x, y);     /* x - y       */
        r2 = _mm_mul_ps(r, r);    /* (x - y)**2  */

        /* horizontal add the components of d2 with */
        /* two instructions. note: it's critical */
        /* here that the last entry of x1 and x2 was 0 */
        /* so that d2.w = 0 */
	s = _mm_add_ps(r2, _mm_movehl_ps(r2, r2));
        s = _mm_add_ss(s, _mm_shuffle_ps(s, s, 1));
        /* store into a regular float. I tried using _mm_rsqrt_ps, but it's not
           accurate to pass the tests */
        _mm_store_ss(&d, s);
        moments_push(&onlinemoments, 1.0 / sqrt((double) d));
    }

    moments[0] = moments_mean(&onlinemoments);
    moments[1] = sqrt(moments_second(&onlinemoments));
    moments[2] = cbrt(moments_third(&onlinemoments));

    return 1;
}
