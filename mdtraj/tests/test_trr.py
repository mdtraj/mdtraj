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

from mdtraj.formats import TRRTrajectoryFile
import os, tempfile
import numpy as np
from nose.tools import assert_raises
from mdtraj.testing import eq, get_fn
from mdtraj.formats import trr

fd, temp = tempfile.mkstemp(suffix='.trr')
os.close(fd)
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_1():
    "Write data and read it back"
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_read_stride():
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz3, time3, step3, box3, lambd3 = f.read(stride=3)
    yield lambda: eq(xyz[::3], xyz3)
    yield lambda: eq(step[::3], step3)
    yield lambda: eq(box[::3], box3)
    yield lambda: eq(time[::3], time3)


def test_read_stride_2():
    "trr read stride when n_frames is supplied (different path)"
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz3, time3, step3, box3, lambd3 = f.read(n_frames=1000, stride=3)
    yield lambda: eq(xyz[::3], xyz3)
    yield lambda: eq(step[::3], step3)
    yield lambda: eq(box[::3], box3)
    yield lambda: eq(time[::3], time3)


def test_15():
    "Write data and read it back"
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(n_frames=500)

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_read_atomindices_1():
    with TRRTrajectoryFile(temp) as f:
        xyz, time, step, box, lambd = f.read()

    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=[0,1,2])
    yield lambda: eq(xyz[:, [0,1,2]], xyz2)
    yield lambda: eq(step, step2)
    yield lambda: eq(box, box2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(time, time2)


def test_read_atomindices_2():
    with TRRTrajectoryFile(temp) as f:
        xyz, time, step, box, lambd = f.read()

    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=slice(None, None, 2))
    yield lambda: eq(xyz[:, ::2], xyz2)
    yield lambda: eq(step, step2)
    yield lambda: eq(box, box2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(time, time2)

def test_2():
    """Write data one frame at a time. This checks how the shape is dealt with,
    because each frame is deficient in shape."""
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        for i in range(len(xyz)):
            f.write(xyz=xyz[i], time=time[i], step=step[i], lambd=lambd[i])
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_ragged_1():
    # try first writing no box vectors,, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    with TRRTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz)
        assert_raises(ValueError, lambda: f.write(xyz, time, box))


def test_ragged_2():
    # try first writing no box vectors, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    with TRRTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time=time, box=box)
        assert_raises(ValueError, lambda: f.write(xyz))


def test_tell():
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek():
    reference = TRRTrajectoryFile(get_fn('frame0.trr')).read()[0]
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        eq(len(f), len(reference))
        eq(len(f.offsets), len(reference))

        eq(f.tell(), 0)
        eq(f.read(1)[0][0], reference[0])
        eq(f.tell(), 1)

        xyz = f.read(1)[0][0]
        eq(xyz, reference[1])
        eq(f.tell(), 2)

        f.seek(0)
        eq(f.tell(), 0)
        xyz = f.read(1)[0][0]
        eq(f.tell(), 1)
        eq(xyz, reference[0])

        f.seek(5)
        eq(f.read(1)[0][0], reference[5])
        eq(f.tell(), 6)

        f.seek(-5, 1)
        eq(f.tell(), 1)
        eq(f.read(1)[0][0], reference[1])


def test_get_velocities():
    """Write data with velocities and read it back"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    vel = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd, vel=vel)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2= f._read(
            n_frames=500, atom_indices=None, get_velocities=True,
            get_forces=False
        )

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(vel, vel2)
    yield lambda: eq(None, forces2)


def test_get_forces():
    """Write data with forces and read it back"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    forces = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd,
                 forces=forces)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2 = f._read(
            n_frames=500, atom_indices=None, get_velocities=False,
            get_forces=True
        )

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(None, vel2)
    yield lambda: eq(forces, forces2)


def test_get_velocities_and_forces():
    """Write data with velocities and forces, and read it back"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    vel = np.array(np.random.randn(500,50,3), dtype=np.float32)
    forces = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd,
                 vel=vel, forces=forces)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2 = f._read(
            n_frames=500, atom_indices=None, get_velocities=True,
            get_forces=True
        )

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(vel, vel2)
    yield lambda: eq(forces, forces2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_get_velocities_forces_atom_indices_1():
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    vel = np.array(np.random.randn(500,50,3), dtype=np.float32)
    forces = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd,
                 vel=vel, forces=forces)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2 = f._read(
            n_frames=500, atom_indices=[0, 1, 2],
            get_velocities=True, get_forces=True
        )

    yield lambda: eq(xyz[:, [0,1,2]], xyz2)
    yield lambda: eq(vel[:, [0,1,2]], vel2)
    yield lambda: eq(forces[:, [0,1,2]], forces2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_get_velocities_forces_atom_indices_2():
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    vel = np.array(np.random.randn(500,50,3), dtype=np.float32)
    forces = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd,
                 vel=vel, forces=forces)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2 = f._read(
            n_frames=500, atom_indices=slice(None, None, 2),
            get_velocities=True, get_forces=True
        )

    yield lambda: eq(xyz[:, ::2], xyz2)
    yield lambda: eq(vel[:, ::2], vel2)
    yield lambda: eq(forces[:, ::2], forces2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)
    pass

def test_read_velocities_do_not_exist():
    """Requesting velocities from a file that does not have them"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, box=box, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        assert_raises(RuntimeError, f._read, n_frames=500,
                      atom_indices=None, get_velocities=True,
                      get_forces=False)


def test_read_forces_do_not_exist():
    """Requesting forces from a file that does not have them"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    box = np.array(np.random.randn(500,3,3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, box=box, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        assert_raises(RuntimeError, f._read, n_frames=500,
                      atom_indices=None, get_velocities=False,
                      get_forces=True)

