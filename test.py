import os
import numpy as np
from mdtraj import xtc, dcd
import numpy.testing as npt
from msmbuilder import Trajectory
import warnings

pdb = './test_data/native.pdb'
fnx = './test_data/frame0.xtc'
fnd = './test_data/frame0.dcd'

tx = 'test.xtc'
td = 'test.dcd'
if os.path.exists(tx):
    os.unlink(tx)
if os.path.exists(td):
    os.unlink(td)


# xtc read

def test_xread0():
    "Reads the same as msmbuilder.Traj with chunk1"
    xyz, time, step, box, prec = xtc.read(fnx, chunk=1)
    t = Trajectory.load_trajectory_file(fnx, Conf=Trajectory.load_trajectory_file(pdb))
    npt.assert_array_almost_equal(xyz, t['XYZList'])

def test_xread1():
    "Reads the same with chunk 100"
    xyz, time, step, box, prec = xtc.read(fnx, chunk=100)
    t = Trajectory.load_trajectory_file(fnx, Conf=Trajectory.load_trajectory_file(pdb))
    npt.assert_array_almost_equal(xyz, t['XYZList'])
        

def test_xread2():
    "Reads the same with chunk 1000"
    xyz, time, step, box, prec = xtc.read(fnx, chunk=1000)
    t = Trajectory.load_trajectory_file(fnx, Conf=Trajectory.load_trajectory_file(pdb))
    npt.assert_array_almost_equal(xyz, t['XYZList'])


# xtc write
def test_xwrite0():
    xyz, time, step, box, prec = xtc.read(fnx)
    xtc.write(tx, xyz=xyz, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(tx)
    
    npt.assert_array_almost_equal(xyz, xyz2)

# xtc write
def test_xwrite1():
    xyz, time, step, box, prec = xtc.read(fnx)
    xtc.write(tx, xyz=xyz, time=time, step=step, box=box, prec=prec, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(tx)

    npt.assert_array_almost_equal(xyz, xyz2, err_msg='xyz wrong')
    npt.assert_array_almost_equal(step, step2, err_msg='step wrong')
    npt.assert_array_almost_equal(box, box2, err_msg='box wrong')
    npt.assert_array_almost_equal(time, time2, err_msg='time wrong')
    npt.assert_array_almost_equal(prec, prec2, err_msg='prec wrong')

def test_xwrite2():
    xyz, time, step, box, prec = xtc.read(fnx)
    time = np.zeros_like(time)
    
    xtc.write(tx, xyz=xyz, time=time, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(tx)

    npt.assert_array_almost_equal(xyz, xyz2, err_msg='xyz wrong')
    npt.assert_array_almost_equal(time, time2, err_msg='time wrong')
   
def test_xwrite3():
    t = Trajectory.load_trajectory_file(fnx, Conf=Trajectory.load_trajectory_file(pdb))
    if os.path.exists(tx):
        os.unlink(tx)
    t.save_to_xtc(tx) 
    
    xyz2, time2, step2, box2, prec2 = xtc.read(tx)

    npt.assert_array_almost_equal(t['XYZList'], xyz2, err_msg='xyz wrong')


# DCD METHOD (READ AND WRITE)

def test_dread():
    "ReadDCD"
    xyz = 0.1*dcd.read_xyz(fnd)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = Trajectory.load_trajectory_file(fnd, Conf=Trajectory.load_trajectory_file(pdb))
    npt.assert_array_almost_equal(xyz, t['XYZList'])
    
def test_dwrite0():
    "Read write read"
    xyz = dcd.read_xyz(fnd)
    dcd.write_xyz(td, xyz, force_overwrite=True)
    xyz2 = dcd.read_xyz(td)
    npt.assert_array_almost_equal(xyz, xyz2)
    
def test_dwrite1():
    "Read, write, read(cytpes)"
    xyz = dcd.read_xyz(fnd)
    dcd.write_xyz(td, xyz, force_overwrite=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = Trajectory.load_trajectory_file(td, Conf=Trajectory.load_trajectory_file(pdb))
    npt.assert_array_almost_equal(0.1*xyz, t['XYZList'])
