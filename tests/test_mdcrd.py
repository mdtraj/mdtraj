##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
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
from mdtraj.formats import MDCRDTrajectoryFile
from mdtraj.testing import eq

fd, temp = tempfile.mkstemp(suffix='.mdcrd')


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.close(fd)
    os.unlink(temp)


def test_read_0(get_fn):
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz, box = f.read()
        assert box is None
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrdbox'), n_atoms=22) as f:
        xyz2, box = f.read()
        eq(box.shape, (1002, 3))
        eq(xyz, xyz2)


def test_read_1(get_fn):
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz, _ = f.read()
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz3, _ = f.read(stride=3)

    eq(xyz[::3], xyz3)


def test_read_write_0():
    xyz = 10 * np.random.randn(100, 11, 3)
    with MDCRDTrajectoryFile(temp, mode='w') as f:
        f.write(xyz)
    with MDCRDTrajectoryFile(temp, n_atoms=11) as f:
        xyz2, _ = f.read()

    eq(_, None)
    eq(xyz, xyz2, decimal=3)


def test_read_write_1():
    xyz = 10 * np.random.randn(100, 11, 3)
    box = np.random.randn(100, 3)
    with MDCRDTrajectoryFile(temp, mode='w') as f:
        f.write(xyz, box)
    with MDCRDTrajectoryFile(temp, n_atoms=11) as f:
        xyz2, box2 = f.read()

    eq(box, box2, decimal=3)
    eq(xyz, xyz2, decimal=3)


def test_read_write_2(get_fn):
    pdb = md.load(get_fn('1bpi.pdb'))
    pdb.save(temp)

    t = md.load(temp, top=pdb.topology)

    eq(t.xyz, pdb.xyz)
    eq(t.unitcell_vectors, pdb.unitcell_vectors)


def test_multiread(get_fn):
    reference = md.load(get_fn('frame0.mdcrd'), top=get_fn('native.pdb'))
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz0, box0 = f.read(n_frames=1)
        xyz1, box1 = f.read(n_frames=1)

    eq(reference.xyz[0], xyz0[0] / 10)
    eq(reference.xyz[1], xyz1[0] / 10)


def test_seek(get_fn):
    reference = md.load(get_fn('frame0.mdcrd'), top=get_fn('native.pdb'))
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        f.seek(1)
        eq(1, f.tell())
        xyz1, box1 = f.read(n_frames=1)
        eq(reference.xyz[1], xyz1[0] / 10)

        f.seek(10)
        eq(10, f.tell())
        xyz10, box10 = f.read(n_frames=1)
        eq(reference.xyz[10], xyz10[0] / 10)
        eq(11, f.tell())

        f.seek(-8, 1)
        xyz3, box3 = f.read(n_frames=1)
        eq(reference.xyz[3], xyz3[0] / 10)

        f.seek(4, 1)
        xyz8, box8 = f.read(n_frames=1)
        eq(reference.xyz[8], xyz8[0] / 10)


def test_atom_indices_0(get_fn):
    atom_indices = np.arange(10)
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz0, box0 = f.read(n_frames=1)
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz1, box1 = f.read(n_frames=1, atom_indices=atom_indices)

    eq(xyz0[:, atom_indices], xyz1)


def test_atom_indices_1(get_fn):
    atom_indices = np.arange(10)
    top = md.load(get_fn('native.pdb'))
    t0 = md.load(get_fn('frame0.mdcrd'), top=top)
    t1 = md.load(get_fn('frame0.mdcrd'), top=top, atom_indices=atom_indices)

    eq(t0.xyz[:, atom_indices], t1.xyz)
