##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp
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


##############################################################################
# imports
##############################################################################

import ast
import sys

import numpy as np
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_approx_equal,
    assert_array_almost_equal,
    assert_array_almost_equal_nulp,
    assert_array_equal,
    assert_array_less,
    assert_array_max_ulp,
    assert_equal,
    assert_raises,
    assert_string_equal,
    assert_warns,
)

# if the system doesn't have scipy, we'd like
# this package to still work:
# we'll just redefine isspmatrix as a function that always returns
# false
try:
    from scipy.sparse import isspmatrix
except ImportError:

    def isspmatrix(x):
        return False


try:
    # need special logic to check for equality of pandas DataFrames.
    # but this is only relevant if the user has pandas installed
    import pandas as pd
except ImportError:
    pass

__all__ = [
    "assert_allclose",
    "assert_almost_equal",
    "assert_approx_equal",
    "assert_array_almost_equal",
    "assert_array_almost_equal_nulp",
    "assert_array_equal",
    "assert_array_less",
    "assert_array_max_ulp",
    "assert_equal",
    "assert_raises",
    "assert_string_equal",
    "assert_warns",
    "get_fn",
    "eq",
    "assert_dict_equal",
    "assert_sparse_matrix_equal",
]


def eq_(a, b, msg=None):
    """Shorthand for 'assert a == b, "%r != %r" % (a, b)

    Parameters
    ----------
    a : object
        The first object
    b : object
        The second object
    msg : string
        Optional assertion message. Why are you using this function then??
    """
    if not a == b:
        raise AssertionError(msg or f"{a!r} != {b!r}")


def get_fn(name):
    raise NotImplementedError(
        "Testing reference data is no longer included in the MDTraj package",
    )


def eq(o1, o2, decimal=6, err_msg=""):
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

    Returns
    -------
    passed : bool
        True if the tests pass. If the tests doesn't pass, since the AssertionError will be raised

    Raises
    ------
    AssertionError
        If the tests fail
    """    
    # assert type(o1) is type(o2), f"o1 and o2 not the same type: {type(o1)} {type(o2)}"
    if isinstance(o1, (np.ndarray, np.ma.MaskedArray)) and isinstance(o2, (np.ndarray, np.ma.MaskedArray)):
        pass
    else:
        assert type(o1) is type(o2), f"o1 and o2 not the same type: {type(o1)} {type(o2)}"

    if isinstance(o1, dict):
        assert_dict_equal(o1, o1, decimal)
    elif isinstance(o1, float):
        np.testing.assert_almost_equal(o1, o2, decimal)
    elif isspmatrix(o1):
        assert_sparse_matrix_equal(o1, o1, decimal)
    elif isinstance(o1, np.ndarray):
        if o1.dtype.kind == "f" or o2.dtype.kind == "f":
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
    elif "pandas" in sys.modules and isinstance(o1, pd.DataFrame):
        # pandas dataframes are basically like dictionaries of numpy arrayss
        assert_dict_equal(o1, o2, decimal=decimal)
    elif isinstance(o1, ast.AST) and isinstance(o2, ast.AST):
        eq_(ast.dump(o1), ast.dump(o2))

    # probably these are other specialized types
    # that need a special check?
    else:
        eq_(o1, o2, msg=err_msg)

    return True


def assert_dict_equal(t1, t2, decimal=6):
    """Assert two dicts are equal. This method should actually
    work for any dict of numpy arrays/objects

    Parameters
    ----------
    t1 : object
    t2 : object
    decimal : int
        Number of decimal places to check, for arrays inside the dicts
    """

    # make sure the keys are the same
    eq_(list(t1.keys()), list(t2.keys()))

    for key, val in t1.items():
        # compare numpy arrays using numpy.testing
        if isinstance(val, np.ndarray) or ("pandas" in sys.modules and isinstance(t1, pd.DataFrame)):
            if val.dtype.kind == "f":
                # compare floats for almost equality
                assert_array_almost_equal(val, t2[key], decimal)
            else:
                # compare everything else (ints, bools) for absolute equality
                assert_array_equal(val, t2[key])
        else:
            eq_(val, t2[key])


def assert_sparse_matrix_equal(m1, m2, decimal=6):
    """Assert two scipy.sparse matrices are equal.

    Parameters
    ----------
    m1 : sparse_matrix
    m2 : sparse_matrix
    decimal : int
        Number of decimal places to check.
    """
    # both are sparse matricies
    assert isspmatrix(m1)
    assert isspmatrix(m1)

    # make sure they have the same format
    eq_(m1.format, m2.format)

    # even though its called assert_array_almost_equal, it will
    # work for scalars
    assert_array_almost_equal((m1 - m2).sum(), 0, decimal=decimal)
