"""
Test the cython dcd module

Note, this file cannot be located in the dcd subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import tempfile, os
from mdtraj import dcd
from mdtraj.testing import get_fn, eq
import warnings

fn_dcd = get_fn('frame0.dcd')
pdb = get_fn('native.pdb')

temp = tempfile.mkstemp(suffix='.dcd')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)

def test_dread():
    #"ReadDCD"
    from msmbuilder import Trajectory
    
    xyz = 0.1*dcd.read_xyz(fn_dcd)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = Trajectory.load_trajectory_file(fn_dcd, Conf=Trajectory.load_trajectory_file(pdb))

    eq(xyz, t['XYZList'])
    
def test_dwrite0():
    #"Read write read"
    xyz = dcd.read_xyz(fn_dcd)
    dcd.write_xyz(temp, xyz, force_overwrite=True)
    xyz2 = dcd.read_xyz(temp)

    eq(xyz, xyz2)
    
def test_dwrite1():
    #"Read, write, read(cytpes)"
    from msmbuilder import Trajectory
    
    xyz = dcd.read_xyz(fn_dcd)
    dcd.write_xyz(temp, xyz, force_overwrite=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = Trajectory.load_trajectory_file(temp, Conf=Trajectory.load_trajectory_file(pdb))

    eq(0.1*xyz, t['XYZList'])