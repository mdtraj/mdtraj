##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors

# Authors: Christoph Klein
# Contributors:
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


import tempfile, os

import numpy as np

import mdtraj as md
from mdtraj.formats import LAMMPSTrajectoryFile
from mdtraj.testing import eq


fd, temp = tempfile.mkstemp(suffix='.lammpstrj')
def teardown_module(module):
    """Remove the temporary file created by tests in this file
    this gets automatically called by pytest. """
    os.close(fd)
    os.unlink(temp)

def test_read_0(get_fn):
    with LAMMPSTrajectoryFile(get_fn('frame0.lammpstrj')) as f:
        xyz, _, _ = f.read()
    with LAMMPSTrajectoryFile(get_fn('frame0.lammpstrj')) as f:
        xyz3, _, _ = f.read(stride=3)
    eq(xyz[::3], xyz3)

def test_read_1(get_fn):
    reference = md.load(get_fn('frame0.dcd'), top=get_fn('native.pdb'))
    traj = md.load(get_fn('frame0.lammpstrj'), top=get_fn('native.pdb'))

    eq(reference.xyz[0], traj.xyz[0], decimal=3)

def test_read_write_0():
    xyz = 10 * np.random.randn(100, 11, 3)
    lengths = np.ones(shape=(100, 3))
    angles = np.empty(shape=(100, 3))
    angles.fill(45)

    with LAMMPSTrajectoryFile(temp, mode='w') as f:
        f.write(xyz, lengths, angles)
    with LAMMPSTrajectoryFile(temp) as f:
        xyz2, new_lengths, new_angles = f.read()

    eq(lengths, new_lengths)
    eq(angles, new_angles)
    eq(xyz, xyz2, decimal=3)

def test_mdwrite(get_fn):
    t = md.load(get_fn('frame0.dcd'), top=get_fn('native.pdb'))
    t.save(temp)

def test_multiread(get_fn):
    reference = md.load(get_fn('frame0.lammpstrj'), top=get_fn('native.pdb'))
    with LAMMPSTrajectoryFile(get_fn('frame0.lammpstrj')) as f:
        xyz0, _, _ = f.read(n_frames=1)
        xyz1, _, _ = f.read(n_frames=1)

    eq(reference.xyz[0], xyz0[0]/10)
    eq(reference.xyz[1], xyz1[0]/10)

def test_seek(get_fn):
    reference = md.load(get_fn('frame0.lammpstrj'), top=get_fn('native.pdb'))

    with LAMMPSTrajectoryFile(get_fn('frame0.lammpstrj')) as f:
        f.seek(1)
        eq(1, f.tell())
        xyz1, _, _ = f.read(n_frames=1)
        eq(reference.xyz[1], xyz1[0]/10)

        f.seek(10)
        eq(10, f.tell())
        xyz10, _, _ = f.read(n_frames=1)
        eq(reference.xyz[10], xyz10[0]/10)
        eq(11, f.tell())

        f.seek(-8, 1)
        xyz3, _, _ = f.read(n_frames=1)
        eq(reference.xyz[3], xyz3[0]/10)

        f.seek(4, 1)
        xyz8, _, _ = f.read(n_frames=1)
        eq(reference.xyz[8], xyz8[0]/10)

def test_custom(get_fn):
    t0 = md.load(get_fn('custom.lammpstrj'), top=get_fn('custom.pdb'))
    t0.save(temp)

    t1 = md.load(temp, top=get_fn('custom.pdb'))
    eq(t0.xyz, t1.xyz)
