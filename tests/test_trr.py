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
from mdtraj.testing import eq
from mdtraj import io
import pytest

fd, temp = tempfile.mkstemp(suffix='.trr')
os.close(fd)


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.unlink(temp)


def test_1():
    # Write data and read it back
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    assert eq(xyz, xyz2)
    assert eq(time, time2)
    assert eq(step, step2)
    assert eq(lambd, lambd2)


def test_read_stride(get_fn):
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz3, time3, step3, box3, lambd3 = f.read(stride=3)
    assert eq(xyz[::3], xyz3)
    assert eq(step[::3], step3)
    assert eq(box[::3], box3)
    assert eq(time[::3], time3)


def test_read_stride_n_frames(get_fn):
    # trr read stride when n_frames is supplied (different path)
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz3, time3, step3, box3, lambd3 = f.read(n_frames=1000, stride=3)
    assert eq(xyz[::3], xyz3)
    assert eq(step[::3], step3)
    assert eq(box[::3], box3)
    assert eq(time[::3], time3)


def test_read_stride_offsets(get_fn):
    # read xtc with stride and offsets
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    for s in (1, 2, 3, 4, 5):
        with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
            f.offsets # pre-compute byte offsets between frames
            xyz_s, time_s, step_s, box_s, lamb_s = f.read(stride=s)
        assert eq(xyz_s, xyz[::s])
        assert eq(step_s, step[::s])
        assert eq(box_s, box[::s])
        assert eq(time_s, time[::s])


def test_read_stride_n_frames_offsets(get_fn):
    # read trr with stride with n_frames and offsets
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    for s in (1, 2, 3, 4, 5):
        with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
            f.offsets # pre-compute byte offsets between frames
            xyz_s, time_s, step_s, box_s, lamb_s = f.read(n_frames=1000, stride=s)
        assert eq(xyz_s, xyz[::s], err_msg='stride=%s' % s)
        assert eq(step_s, step[::s])
        assert eq(box_s, box[::s])
        assert eq(time_s, time[::s])


def test_read_stride_switching(get_fn):
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        f.offsets  # pre-compute byte offsets between frames
        # read the first 10 frames with stride of 2
        s = 2
        n_frames = 10
        xyz_s, time_s, step_s, box_s, lamb_s = f.read(n_frames=n_frames, stride=s)
        assert eq(xyz_s, xyz[:n_frames*s:s])
        assert eq(step_s, step[:n_frames*s:s])
        assert eq(box_s, box[:n_frames*s:s])
        assert eq(time_s, time[:n_frames*s:s])
        # now read the rest with stride 3, should start from frame index 8.
        # eg. np.arange(0, n_frames*s + 1, 2)[-1] == 18
        offset = f.tell()
        assert offset == 20
        s = 3
        xyz_s, time_s, step_s, box_s, lamb_s = f.read(n_frames=None, stride=s)
        assert eq(xyz_s, xyz[offset::s])
        assert eq(step_s, step[offset::s])
        assert eq(box_s, box[offset::s])
        assert eq(time_s, time[offset::s])


def test_write_read():
    # Write data and read it back
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(n_frames=500)

    assert eq(xyz, xyz2)
    assert eq(time, time2)
    assert eq(step, step2)
    assert eq(lambd, lambd2)


def test_read_atomindices_1():
    with TRRTrajectoryFile(temp) as f:
        xyz, time, step, box, lambd = f.read()

    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=[0, 1, 2])
    assert eq(xyz[:, [0, 1, 2]], xyz2)
    assert eq(step, step2)
    assert eq(box, box2)
    assert eq(lambd, lambd2)
    assert eq(time, time2)


def test_read_atomindices_2():
    with TRRTrajectoryFile(temp) as f:
        xyz, time, step, box, lambd = f.read()

    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=slice(None, None, 2))
    assert eq(xyz[:, ::2], xyz2)
    assert eq(step, step2)
    assert eq(box, box2)
    assert eq(lambd, lambd2)
    assert eq(time, time2)


def test_deficient_shape():
    # Write data one frame at a time. This checks how the shape is dealt with,
    # because each frame is deficient in shape.
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        for i in range(len(xyz)):
            f.write(xyz=xyz[i], time=time[i], step=step[i], lambd=lambd[i])
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    assert eq(xyz, xyz2)
    assert eq(time, time2)
    assert eq(step, step2)
    assert eq(lambd, lambd2)


def test_ragged_1():
    # try first writing no box vectors,, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    with TRRTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz)
        with pytest.raises(ValueError):
            f.write(xyz, time, box)


def test_ragged_2():
    # try first writing no box vectors, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    with TRRTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time=time, box=box)
        with pytest.raises(ValueError):
            f.write(xyz)


def test_malformed():
    with open(temp, 'w') as tmpf:
        tmpf.write("foo")  # very badly malformed TRR

    with pytest.raises(IOError):
        TRRTrajectoryFile(temp)

    psutil = pytest.importorskip("psutil")
    open_files = psutil.Process().open_files()
    paths = [os.path.realpath(f.path) for f in open_files]
    assert os.path.realpath(temp) not in paths


def test_tell(get_fn):
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek(get_fn):
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
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    vel = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f._write(xyz=xyz, time=time, step=step, box=box, lambd=lambd, vel=vel)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2, vel2, forces2 = f._read(
            n_frames=500, atom_indices=None, get_velocities=True,
            get_forces=False
        )

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(lambd, lambd2)
    eq(vel, vel2)
    eq(None, forces2)


def test_get_forces():
    """Write data with forces and read it back"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    forces = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
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

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(lambd, lambd2)
    eq(None, vel2)
    eq(forces, forces2)


def test_get_velocities_and_forces():
    """Write data with velocities and forces, and read it back"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    vel = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    forces = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
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

    eq(xyz, xyz2)
    eq(vel, vel2)
    eq(forces, forces2)
    eq(time, time2)
    eq(step, step2)
    eq(lambd, lambd2)


def test_get_velocities_forces_atom_indices_1():
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    vel = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    forces = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
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

    eq(xyz[:, [0, 1, 2]], xyz2)
    eq(vel[:, [0, 1, 2]], vel2)
    eq(forces[:, [0, 1, 2]], forces2)
    eq(time, time2)
    eq(step, step2)
    eq(lambd, lambd2)


def test_get_velocities_forces_atom_indices_2():
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    vel = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    forces = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
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

    eq(xyz[:, ::2], xyz2)
    eq(vel[:, ::2], vel2)
    eq(forces[:, ::2], forces2)
    eq(time, time2)
    eq(step, step2)
    eq(lambd, lambd2)
    pass


def test_read_velocities_do_not_exist():
    """Requesting velocities from a file that does not have them"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, box=box, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        with pytest.raises(RuntimeError):
            f._read(n_frames=500, atom_indices=None, get_velocities=True, get_forces=False)


def test_read_forces_do_not_exist():
    """Requesting forces from a file that does not have them"""
    # NOTE: this is a test of a hidden API
    xyz = np.array(np.random.randn(500, 50, 3), dtype=np.float32)
    box = np.array(np.random.randn(500, 3, 3), dtype=np.float32)
    time = np.array(np.random.randn(500), dtype=np.float32)
    step = np.array(np.arange(500), dtype=np.int32)
    lambd = np.array(np.random.randn(500), dtype=np.float32)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, box=box, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        with pytest.raises(RuntimeError):
            f._read(n_frames=500, atom_indices=None, get_velocities=False, get_forces=True)
