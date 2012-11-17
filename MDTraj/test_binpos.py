"""
Test the cython binpos reader

Note, this file cannot be located in the binpos subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

from mdtraj import binpos, dcd
from mdtraj.testing import get_fn, eq

fn_binpos = get_fn('frame0.binpos')
fn_dcd = get_fn('frame0.dcd')

def test_binpos_read0():
    xyz = binpos.read_xyz(fn_binpos)
    xyz2 = dcd.read_xyz(fn_dcd)
    eq(xyz[1:], xyz2)

def test_binpos_read1():
    xyz = binpos.read_xyz(fn_binpos, chunk=11)
    xyz2 = dcd.read_xyz(fn_dcd)
    eq(xyz[1:], xyz2)
    
def test_binpos_read2():
    xyz = binpos.read_xyz(fn_binpos, chunk=1000)
    xyz2 = dcd.read_xyz(fn_dcd)
    eq(xyz[1:], xyz2)
