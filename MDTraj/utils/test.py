"""Tests for some of the utilities
"""
##############################################################################
# imports
##############################################################################

import numpy as np
from .arrays import ensure_type, TypeCastPerformanceWarning
from .unit import in_units_of
from mdtraj.testing import raises
import warnings

##############################################################################
# globals
##############################################################################

a = np.ones(10, dtype=np.float32)
b = np.ones((10,10), dtype=np.float64, order='F')

##############################################################################
# tests
##############################################################################


def test_ensure_type_1():
    ensure_type(a, np.float32, 1, '', length=10)


def test_ensure_type_2():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        val = ensure_type(a, np.float64, 1, '', length=10)

        assert val.dtype == np.float64
        assert a.dtype == np.float32  # a should not be changed
        assert len(w) == 1
        assert issubclass(w[-1].category, TypeCastPerformanceWarning)


def test_ensure_type_25():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        val = ensure_type(a, np.float64, 1, '', length=10, warn_on_cast=False)

        assert val.dtype == np.float64
        assert a.dtype == np.float32  # a should not be changed
        assert len(w) == 0  # no warning since we set warn_on_cast to False


@raises(ValueError)
def test_ensure_type_3():
    ensure_type(a, np.float32, 1, '', length=11)


def test_ensure_type_4():
    ensure_type(None, np.float64, 1, '', length=11, can_be_none=True)


@raises(ValueError)
def test_ensure_type_5():
    ensure_type(a, np.float32, 1, '', length=11, can_be_none=True)


def test_ensure_type_6():
    val = ensure_type(b, np.float64, 2, '', shape=(10,10))
    assert val.flags.c_contiguous is True


def test_ensure_type_7():
    c = ensure_type(a, np.float32, ndim=2, name='', add_newaxis_on_deficient_ndim=True)
    assert c.shape == (1, len(a))


def test_ensure_type_8():
    c = ensure_type(np.zeros((5,10)), np.float32, ndim=2, name='', shape=(None, 10))
    assert c.shape == (5, 10)


@raises(ValueError)
def test_ensure_type_9():
    c = ensure_type(np.zeros((5,11)), np.float32, ndim=2, name='', shape=(None, 10))


def test_unit_1():
    assert 1 == in_units_of(100, 'meter', 'centimeter')

def test_unit_2():
    a = in_units_of(1, 'meter**2/second', 'nanometers**2/picosecond')
    b = 1e-6
    assert abs(a-b) < 1e-10
