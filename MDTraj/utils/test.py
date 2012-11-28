import numpy as np
from arrays import ensure_type, TypeCastPerformanceWarning
from mdtraj.testing import raises
import warnings

a = np.ones(10, dtype=np.float32)
b = np.ones((10,10), dtype=np.float64, order='F')

def test_ensure_type_1():
    ensure_type(a, np.float32, 1, '', length=10)

def test_ensure_type_2():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        val = ensure_type(a, np.float64, 1, '', length=10)

        assert val.dtype == np.float64
        assert a.dtype == np.float32 # a should not be changed
        assert len(w) == 1
        assert issubclass(w[-1].category, TypeCastPerformanceWarning)

def test_ensure_type_25():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        val = ensure_type(a, np.float64, 1, '', length=10, warn_on_cast=False)

        assert val.dtype == np.float64
        assert a.dtype == np.float32 # a should not be changed
        assert len(w) == 0 # no warning since we set warn_on_cast to False

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
