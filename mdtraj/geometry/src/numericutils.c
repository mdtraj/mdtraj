/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2012-2013 Stanford University and the Authors               */
/*                                                                       */
/* Authors: Christoph Klein                                              */
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

#include <math.h>


int histogram(const float* data, const int n_data, const float* bins,
              const int n_bins, const float min_bin, const float max_bin,
              int* out)
{
    /* Compute the histogram of a set of data. */
    int i, bin;
    float binwidth = (max_bin - min_bin) / n_bins;
    printf("%.1f %.1f %d  %.1f\n\n", max_bin, min_bin, n_bins, binwidth);

    for (i = 0; i < n_data; ++i)
    {
        printf("%.1f %.1f   ", data[i], min_bin);
        bin = floor((data[i] - min_bin)) / binwidth;
        printf("%d\n", bin);
        out[bin]++;
    }

    return 1;
}
