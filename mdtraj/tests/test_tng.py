import mdtraj as md
from mdtraj.formats import TNGTrajectoryFile
from mdtraj.testing.testing import *
import os
import tempfile
import numpy as np

test_fn = get_fn('tng_example.tng')
fd, temp = tempfile.mkstemp(suffix='.xtc')
os.close(fd)

def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)

def test_load_trajectory():
    """Compare a TNG file to the PDB file it was created from."""
    
    pdbtraj = md.load_pdb(get_fn('frame0.pdb'))
    tngtraj = md.load_tng(get_fn('frame0.tng'), top=pdbtraj.topology)
    eq(pdbtraj.n_frames, tngtraj.n_frames)
    eq(pdbtraj.unitcell_vectors, tngtraj.unitcell_vectors)
    eq(pdbtraj.xyz, tngtraj.xyz)

def test_write():
    """Write a TNG file, then read it back."""
    xyz = np.asarray(np.around(np.random.randn(100, 10, 3), 3), dtype=np.float32)
    time = 0.1*np.arange(0, 100, dtype=np.float32)
    box = np.asarray(np.random.randn(100,3,3), dtype=np.float32)

    with TNGTrajectoryFile(temp, 'w') as f:
        f.write(xyz, time=time, box=box)
    with TNGTrajectoryFile(temp) as f:
        xyz2, time2, box2 = f.read()
    eq(xyz, xyz2)
    eq(time, time2)
    eq(box, box2)
