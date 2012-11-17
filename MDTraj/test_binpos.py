"""
Test the cython binpos module

Note, this file cannot be located in the binpos subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""
import tempfile, os
from mdtraj import binpos, dcd, io
from mdtraj.testing import get_fn, eq
import numpy as np

# frame0.binpos was generated from frame0.dcd using
# VMD. The binpos file has one more frame than the dcd file
# because VMD also saved the PDB coordinates
fn_binpos = get_fn('frame0.binpos')
fn_dcd = get_fn('frame0.dcd')

temp = tempfile.mkstemp(suffix='.dcd')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)

def test_read_chunk1():
    xyz = binpos.read_xyz(fn_binpos)
    xyz2 = dcd.read_xyz(fn_dcd)
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')

    eq(xyz[1:], xyz2)
    eq(xyz, xyz3)

def test_read_chunk10():
    xyz = binpos.read_xyz(fn_binpos, chunk=10)
    xyz2 = dcd.read_xyz(fn_dcd)
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    
    eq(xyz[1:], xyz2)
    eq(xyz, xyz3)

def test_read_chunk1000():
    xyz = binpos.read_xyz(fn_binpos, chunk=1000)
    xyz2 = dcd.read_xyz(fn_dcd)
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    
    eq(xyz[1:], xyz2)
    eq(xyz, xyz3)

def test_write_1():
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    binpos.write_xyz(temp, xyz, force_overwrite=True)
    xyz2 = binpos.read_xyz(temp)
    
    eq(xyz, xyz2)