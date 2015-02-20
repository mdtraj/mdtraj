##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

from __future__ import print_function, division

import collections

import cython
import numpy as np

##############################################################################
# Headers
##############################################################################

cdef extern from 'numericutils.h':
    int histogram(const float* a, const float* bin_edges, float* out) nogil

##############################################################################
# Wrappers
##############################################################################

def _histogram(a, bins=10, bin_range=None):
    """histogram(a, bins=10, bin_range=None)

    Compute the histogram of a set of data.

    Parameters
    ----------
    a : array_like
        Input data. The histogram is computed over the flattened array.
    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given bin_range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    bin_range : (float, float), optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(a.min(), a.max())``.  Values outside the bin_range are
        ignored.

    Returns
    -------
    hist : array
        The values of the histogram. See `normed` and `weights` for a
        description of the possible semantics.
    bin_edges : array of dtype float
        Return the bin edges ``(length(hist)+1)``.
    """
    a = np.asarray(a)
    a = a.ravel()

    if bin_range is not None:
        mn, mx = bin_range
        if mn > mx:
            raise AttributeError(
                'max must be larger than min in bin_range parameter.')

    if not isinstance(bins, collections.Iterable):
        if np.isscalar(bins) and bins < 1:
            raise ValueError(
                '`bins` should be a positive integer.')
        if bin_range is None:
            if a.size == 0:
                # Handle empty arrays. Can't determine range, so use 0-1.
                bin_range = (0, 1)
            else:
                bin_range = (a.min(), a.max())
        mn, mx = [mi + 0.0 for mi in bin_range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins + 1, endpoint=True)
    else:
        bins = np.asarray(bins)
        if (np.diff(bins) < 0).any():
            raise AttributeError(
                'bins must increase monotonically.')

    bin_edges = [0, 1, 2, 3, 4, 5]

    cdef float[::1] ca = np.ascontiguousarray(a, dtype=np.float32)
    cdef float[::1] cbin_edges = np.ascontiguousarray(bin_edges, dtype=np.float32)
    cdef float[::1] out = np.zeros(cbin_edges.shape[0] - 1, dtype=np.float32)
    _hist(ca, cbin_edges, out)
    return np.asarray(out), np.asarray(cbin_edges)


@cython.boundscheck(False)
def _hist(float[::1] ca,
               float[::1] cbin_edges,
               float[::1] out):
    histogram(&ca[0], &cbin_edges[0], &out[0])
