/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2014-2016 Stanford University and the Authors               */
/*                                                                       */
/* Authors: Robert McGibbon, Peter Eastman                               */
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
#include "vectorize.h"
#include "msvccompat.h"
#include "moments.h"
#include "math_patch.h"

void drid_moments(float* coords, int32_t index, int32_t* partners, int32_t n_partners, double* moments)
{
    moments_t onlinemoments;
    moments_clear(&onlinemoments);
    fvec4 x(coords[3*index], coords[3*index+1], coords[3*index+2], 0);

    for (int i = 0; i < n_partners; i++) {
        fvec4 y(coords[3*partners[i]], coords[3*partners[i]+1], coords[3*partners[i]+2], 0);
        fvec4 r = x-y;
        float d = dot3(r, r);
        moments_push(&onlinemoments, 1.0/sqrt((double) d));
    }

    moments[0] = moments_mean(&onlinemoments);
    moments[1] = sqrt(moments_second(&onlinemoments));
    moments[2] = cbrt(moments_third(&onlinemoments));
}
