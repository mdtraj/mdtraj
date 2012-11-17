"""
Test the cython dcd module

Note, this file cannot be located in the dcd subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import tempfile, os
import numpy as np
from mdtraj import dcd, io
from mdtraj.testing import get_fn, eq
import warnings

fn_dcd = get_fn('frame0.dcd')
pdb = get_fn('native.pdb')

temp = tempfile.mkstemp(suffix='.dcd')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)

def test_read():
    xyz = dcd.read_xyz(fn_dcd)
    xyz2 = io.loadh(get_fn('frame0.dcd.h5'), 'xyz')

    eq(xyz, xyz2)
    
def test_write_0():
    xyz = dcd.read_xyz(fn_dcd)
    dcd.write_xyz(temp, xyz, force_overwrite=True)
    xyz2 = dcd.read_xyz(temp)

    eq(xyz, xyz2)
    
def test_write_1():
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    dcd.write_xyz(temp, xyz, force_overwrite=True)
    xyz2 = dcd.read_xyz(temp)

    eq(xyz, xyz2)