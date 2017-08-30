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
import pytest

fd, temp = tempfile.mkstemp(suffix='.trr')
os.close(fd)


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
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
