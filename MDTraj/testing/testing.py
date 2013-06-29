# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

##############################################################################
# imports
##############################################################################

import os
import functools
import numpy as np
from numpy.testing import (assert_allclose, assert_almost_equal,
  assert_approx_equal, assert_array_almost_equal, assert_array_almost_equal_nulp,
  assert_array_equal, assert_array_less, assert_array_max_ulp, assert_equal,
  assert_raises, assert_string_equal, assert_warns)
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises
from nose import SkipTest
from pkg_resources import resource_filename


# if the system doesn't have scipy, we'd like
# this package to still work:
# we'll just redefine isspmatrix as a function that always returns
# false
try:
    from scipy.sparse import isspmatrix
except ImportError:
    isspmatrix = lambda x: False
    
    
__all__ = ['assert_allclose', 'assert_almost_equal', 'assert_approx_equal',
           'assert_array_almost_equal', 'assert_array_almost_equal_nulp',
           'assert_array_equal', 'assert_array_less', 'assert_array_max_ulp',
           'assert_equal', 'assert_raises', 'assert_string_equal', 'assert_warns',
           'get_fn', 'eq', 'assert_dict_equal', 'assert_spase_matrix_equal',
           'expected_failure', 'skip', 'ok_', 'eq_', 'raises', 'skipif']

##############################################################################
# functions
##############################################################################


def get_fn(name):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``MDTraj/testing/reference``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the reference/ folder).

    Examples
    --------
    >>> import mdtraj as md
    >>> t = md.load(get_fn('2EQQ.pdb'))
    >>> eq(t.n_frames, 20)    # this runs the assert, using the eq() func.
    """
    
    fn = resource_filename('mdtraj', os.path.join('testing/reference', name))

    if not os.path.exists(fn):
        raise ValueError('Sorry! %s does not exists. If you just '
            'added it, you\'ll have to re install' % fn)

    return fn


def eq(o1, o2, decimal=6, err_msg=''):
    """Convenience function for asserting that two objects are equal to one another

    If the objects are both arrays or sparse matrices, this method will
    dispatch to an appropriate handler, which makes it a little bit more
    useful than just calling ``assert o1 == o2`` (which wont work for numpy
    arrays -- it returns an array of bools, not a single True or False)

    Parameters
    ----------
    o1 : object
        The first object
    o2 : object
        The second object
    decimal : int
        If the two objects are floats or arrays of floats, they'll be checked for
        equality up to this decimal place.
    err_msg : str
        Custom error message
    """
    assert (type(o1) is type(o2)), 'o1 and o2 not the same type: %s %s' % (type(o1), type(o2))

    if isinstance(o1, dict):
        assert_dict_equal(o1, o1, decimal)
    elif isinstance(o1, float):
        np.testing.assert_almost_equal(o1, o2, decimal)
    elif isspmatrix(o1):
        assert_spase_matrix_equal(o1, o1, decimal)
    elif isinstance(o1, np.ndarray):
        if o1.dtype.kind == 'f' or o2.dtype.kind == 'f':
            # compare floats for almost equality
            assert_array_almost_equal(o1, o2, decimal, err_msg=err_msg)
        elif o1.dtype.type == np.core.records.record:
            # if its a record array, we need to comparse each term
            assert o1.dtype.names == o2.dtype.names
            for name in o1.dtype.names:
                eq(o1[name], o2[name], decimal=decimal, err_msg=err_msg)
        else:
            # compare everything else (ints, bools) for absolute equality
            assert_array_equal(o1, o2, err_msg=err_msg)
    # probably these are other specialized types
    # that need a special check?
    else:
        eq_(o1, o2)


def assert_dict_equal(t1, t2, decimal=6):
    """
    Assert two dicts are equal.
    This method should actually
    work for any dict of numpy arrays/objects
    """

    # make sure the keys are the same
    eq_(t1.keys(), t2.keys())

    for key, val in t1.iteritems():
        # compare numpy arrays using numpy.testing
        if isinstance(val, np.ndarray):
            if val.dtype.kind ==  'f':
                # compare floats for almost equality
                assert_array_almost_equal(val, t2[key], decimal)
            else:
                # compare everything else (ints, bools) for absolute equality
                assert_array_equal(val, t2[key])
        else:
            eq_(val, t2[key])


def assert_spase_matrix_equal(m1, m2, decimal=6):
    """Assert two scipy.sparse matrices are equal."""
    # both are sparse matricies
    assert isspmatrix(m1)
    assert isspmatrix(m1)

    # make sure they have the same format
    eq_(m1.format, m2.format)

    # even though its called assert_array_almost_equal, it will
    # work for scalars
    assert_array_almost_equal((m1 - m2).sum(), 0, decimal=decimal)


# decorator to mark tests as expected failure
def expected_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except BaseException:
            raise SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner


# decorator to skip tests
def skip(reason):
    def wrap(test):
        @functools.wraps(test)
        def inner(*args, **kwargs):
            raise SkipTest
            print "After f(*args)"
        return inner
    return wrap
