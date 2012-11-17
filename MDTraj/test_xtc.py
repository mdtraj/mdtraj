"""
Test the cython xtc module

Note, this file cannot be located in the xtc subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import os, tempfile
import numpy as np
from mdtraj import xtc
from mdtraj.testing import get_fn, eq

fn_xtc = get_fn('frame0.xtc')
pdb = get_fn('native.pdb')

temp = tempfile.mkstemp(suffix='.xtc')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)


def test_xread0():
    from msmbuilder import Trajectory
    #"Reads the same as msmbuilder.Traj with chunk1"
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=1)
    t = Trajectory.load_trajectory_file(fn_xtc, Conf=Trajectory.load_trajectory_file(pdb))
    eq(xyz, t['XYZList'])

def test_xread1():
    from msmbuilder import Trajectory
    #"Reads the same with chunk 100"
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=100)
    t = Trajectory.load_trajectory_file(fn_xtc, Conf=Trajectory.load_trajectory_file(pdb))
    eq(xyz, t['XYZList'])
        

def test_xread2():
    from msmbuilder import Trajectory
    #"Reads the same with chunk 1000"
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=1000)
    t = Trajectory.load_trajectory_file(fn_xtc, Conf=Trajectory.load_trajectory_file(pdb))
    eq(xyz, t['XYZList'])


# xtc write
def test_xwrite0():
    xyz, time, step, box, prec = xtc.read(fn_xtc)
    xtc.write(temp, xyz=xyz, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)
    
    eq(xyz, xyz2)

# xtc write
def test_xwrite1():
    xyz, time, step, box, prec = xtc.read(fn_xtc)
    xtc.write(temp, xyz=xyz, time=time, step=step, box=box, prec=prec, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)

    eq(xyz, xyz2)
    eq(step, step2)
    eq(box, box2)
    eq(time, time2)
    eq(prec, prec2)

def test_xwrite2():
    xyz, time, step, box, prec = xtc.read(fn_xtc)
    time = np.zeros_like(time)
    
    xtc.write(temp, xyz=xyz, time=time, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)

    eq(xyz, xyz2)
    eq(time, time2)
   
def test_xwrite3():
    from msmbuilder import Trajectory
    t = Trajectory.load_trajectory_file(fn_xtc, Conf=Trajectory.load_trajectory_file(pdb))
    if os.path.exists(temp):
        os.unlink(temp)
    t.save_to_xtc(temp) 
    
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)

    eq(t['XYZList'], xyz2)
