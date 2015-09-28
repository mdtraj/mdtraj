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


# Portions of this code were copied from the Numpy source
# those portions are owned by the numpy developers, and released
# under the following license:

# Copyright 2005-2012, NumPy Developers.
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# Neither the name of the NumPy Developers nor the names of any contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.


from __future__ import print_function

import numpy as np

from mdtraj.geometry import _numericutils
from mdtraj.testing import eq, assert_raises


def test_simple1():
    n = 100
    v = np.random.rand(n)
    (a, b) = _numericutils._histogram(v)

    # Matches numpy
    a_np, b_np = np.histogram(v)
    eq(a, a_np)
    eq(b, b_np)

    # Check if the sum of the bins equals the number of samples.
    eq(np.sum(a, axis=0), np.int64(n))


def test_simple2():
    # Check that the bin counts are evenly spaced when the data is from a
    # linear function.
    v = np.linspace(0, 10, 100)
    (a, b) = _numericutils._histogram(v)

    # Matches numpy
    a_np, b_np = np.histogram(v)
    eq(a, a_np)
    eq(b, b_np)

    eq(a, np.ones(10) * 10)


def test_one_bin():
    # Ticket 632
    v = [1, 2, 3, 4]
    bins = [1, 2]
    a, b = _numericutils._histogram(v, bins)

    # Matches numpy
    a_np, b_np = np.histogram(v, bins)
    eq(a, a_np)
    eq(b, b_np)

    eq(a, np.array([2]))
    eq(b, np.array([1, 2]))

    assert_raises(ValueError, _numericutils._histogram, [1, 2], bins=0)

    v = [1, 2]
    bins = 1
    a, b = _numericutils._histogram(v, bins)

    # Matches numpy
    a_np, b_np = np.histogram(v, bins)
    eq(a, a_np)
    eq(b, b_np)

    eq(a, np.array([2]))
    eq(b, np.array([1., 2.]))


# def test_normed():
#     # Check that the integral of the density equals 1.
#     n = 100
#     v = rand(n)
#     a, b = histogram(v, normed=True)
#     area = np.sum(a * diff(b))
#     assert_almost_equal(area, 1)
#
#     # Check with non-constant bin widths (buggy but backwards compatible)
#     v = np.arange(10)
#     bins = [0, 1, 5, 9, 10]
#     a, b = histogram(v, bins, normed=True)
#     area = np.sum(a * diff(b))
#     assert_almost_equal(area, 1)


# def test_density():
#     # Check that the integral of the density equals 1.
#     n = 100
#     v = rand(n)
#     a, b = histogram(v, density=True)
#     area = np.sum(a * diff(b))
#     assert_almost_equal(area, 1)
#
#     # Check with non-constant bin widths
#     v = np.arange(10)
#     bins = [0, 1, 3, 6, 10]
#     a, b = histogram(v, bins, density=True)
#     assert_array_equal(a, .1)
#     assert_equal(np.sum(a*diff(b)), 1)
#
#     # Variale bin widths are especially useful to deal with
#     # infinities.
#     v = np.arange(10)
#     bins = [0, 1, 3, 6, np.inf]
#     a, b = histogram(v, bins, density=True)
#     assert_array_equal(a, [.1, .1, .1, 0.])
#
#     # Taken from a bug report from N. Becker on the numpy-discussion
#     # mailing list Aug. 6, 2010.
#     counts, dmy = np.histogram(
#         [1, 2, 3, 4], [0.5, 1.5, np.inf], density=True)
#     assert_equal(counts, [.25, 0])


def test_outliers():
    # Check that outliers are not tallied
    v = np.arange(10) + .5

    # Lower outliers
    a, b = _numericutils._histogram(v, bin_range=[0, 9])
    eq(a.sum(), np.int64(9))

    # Matches numpy
    a_np, b_np = np.histogram(v, range=[0, 9])
    eq(a, a_np)
    eq(b, b_np)

    # Upper outliers
    a, b = _numericutils._histogram(v, bin_range=[1, 10])
    eq(a.sum(), np.int64(9))

    # Matches numpy
    a_np, b_np = np.histogram(v, range=[1, 10])
    eq(a, a_np)
    eq(b, b_np)

    # # Normalization
    # a, b = histogram(a, range=[1, 9], normed=True)
    # assert_almost_equal((a * diff(b)).sum(), 1, decimal=15)
    #
    # # Weights
    # w = np.arange(10) + .5
    # h, b = histogram(a, range=[1, 9], weights=w, normed=True)
    # assert_equal((h * diff(b)).sum(), 1)
    #
    # h, b = histogram(a, bins=8, range=[1, 9], weights=w)
    # assert_equal(h, w[1:-1])


# def test_type():
#     # Check the type of the returned histogram
#     a = np.arange(10) + .5
#     h, b = histogram(a)
#     assert_(issubdtype(h.dtype, int))
#
#     h, b = histogram(a, normed=True)
#     assert_(issubdtype(h.dtype, float))
#
#     h, b = histogram(a, weights=np.ones(10, int))
#     assert_(issubdtype(h.dtype, int))
#
#     h, b = histogram(a, weights=np.ones(10, float))
#     assert_(issubdtype(h.dtype, float))


def test_f32_rounding():
    # gh-4799, check that the rounding of the edges works with float32
    v = np.array([276.318359  , -69.593948  , 21.329449], dtype=np.float32)
    bins = 100
    a, b = _numericutils._histogram(v, bins)
    eq(a.sum(), np.int64(3))

    # Matches numpy
    a_np, b_np = np.histogram(v, bins)
    eq(a, a_np)
    eq(b, b_np, decimal=5)  # TODO: Fix round off error.

    v = np.array([5005.689453, 4481.327637, 6010.369629], dtype=np.float32)
    bins = 100
    a, b = _numericutils._histogram(v, bins)
    eq(a.sum(), np.int64(3))

    # Matches numpy
    a_np, b_np = np.histogram(v, bins)
    eq(a, a_np)
    eq(b, b_np, decimal=3)  # TODO: Fix round off error.

# def test_weights():
#     v = rand(100)
#     w = np.ones(100) * 5
#     a, b = histogram(v)
#     na, nb = histogram(v, normed=True)
#     wa, wb = histogram(v, weights=w)
#     nwa, nwb = histogram(v, weights=w, normed=True)
#     assert_array_almost_equal(a * 5, wa)
#     assert_array_almost_equal(na, nwa)
#
#     # Check weights are properly applied.
#     v = np.linspace(0, 10, 10)
#     w = np.concatenate((np.zeros(5), np.ones(5)))
#     wa, wb = histogram(v, bins=np.arange(11), weights=w)
#     assert_array_almost_equal(wa, w)
#
#     # Check with integer weights
#     wa, wb = histogram([1, 2, 2, 4], bins=4, weights=[4, 3, 2, 1])
#     assert_array_equal(wa, [4, 5, 0, 1])
#     wa, wb = histogram(
#         [1, 2, 2, 4], bins=4, weights=[4, 3, 2, 1], normed=True)
#     assert_array_almost_equal(wa, np.array([4, 5, 0, 1]) / 10. / 3. * 4)
#
#     # Check weights with non-uniform bin widths
#     a, b = histogram(
#         np.arange(9), [0, 1, 3, 6, 10],
#         weights=[2, 1, 1, 1, 1, 1, 1, 1, 1], density=True)
#     assert_almost_equal(a, [.2, .1, .1, .075])


def test_empty():
    a, b = _numericutils._histogram([], bins=([0, 1]))
    eq(a, np.array([0]))
    eq(b, np.array([0, 1]))

if __name__ == '__main__':
    test_simple1()