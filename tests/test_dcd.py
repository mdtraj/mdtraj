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


import numpy as np
from mdtraj.formats import DCDTrajectoryFile
from mdtraj import io
from mdtraj.testing import eq
import pytest


def test_read(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    xyz, box_lengths, box_angles = DCDTrajectoryFile(fn_dcd).read()
    xyz2 = io.loadh(get_fn('frame0.dcd.h5'), 'xyz')

    eq(xyz, xyz2)


def test_read_2(get_fn):
    # check nframes
    fn_dcd = get_fn('frame0.dcd')
    xyz1, box_lengths1, box_angles1 = DCDTrajectoryFile(fn_dcd).read()
    xyz2, box_lengths2, box_angles2 = DCDTrajectoryFile(fn_dcd).read(10000)

    assert eq(xyz1, xyz2)
    assert eq(box_lengths1, box_lengths2)
    assert eq(box_angles1, box_angles2)


def test_read_stride(get_fn):
    # Read dcd with stride
    fn_dcd = get_fn('frame0.dcd')
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz1, box_lengths1, box_angles1 = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2, box_lengths2, box_angles2 = f.read(stride=2)

    assert eq(xyz1[::2], xyz2)
    assert eq(box_lengths1[::2], box_lengths2)
    assert eq(box_angles1[::2], box_angles2)


def test_read_stride_2(get_fn):
    # Read dcd with stride when n_frames is supplied (different code path)
    fn_dcd = get_fn('frame0.dcd')
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz1, box_lengths1, box_angles1 = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2, box_lengths2, box_angles2 = f.read(n_frames=1000, stride=2)

    assert eq(xyz1[::2], xyz2)
    assert eq(box_lengths1[::2], box_lengths2)
    assert eq(box_angles1[::2], box_angles2)


def test_read_3(get_fn):
    # Check streaming read of frames 1 at a time
    fn_dcd = get_fn('frame0.dcd')
    xyz_ref, box_lengths_ref, box_angles_ref = DCDTrajectoryFile(fn_dcd).read()

    reader = DCDTrajectoryFile(fn_dcd)
    for i in range(len(xyz_ref)):
        xyz, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)


def test_read_4(get_fn):
    # Check streaming read followed by reading the 'rest'
    fn_dcd = get_fn('frame0.dcd')
    xyz_ref, box_lengths_ref, box_angles_ref = DCDTrajectoryFile(fn_dcd).read()

    reader = DCDTrajectoryFile(fn_dcd)
    for i in range(len(xyz_ref) // 2):
        xyz, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)

    xyz_rest, box_rest, angles_rest = reader.read()
    i = len(xyz_ref) // 2
    assert eq(xyz_ref[i:], xyz_rest)
    assert eq(box_lengths_ref[i:], box_rest)
    assert eq(box_angles_ref[i:], angles_rest)

    assert len(xyz_ref) == i + len(xyz_rest)


def test_read_5(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz_ref, box_lengths_ref, box_angles_ref = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz, box_lengths, box_angles = f.read(atom_indices=[1, 2, 5])

    assert eq(xyz_ref[:, [1, 2, 5], :], xyz)


def test_read_6(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz_ref, box_lengths_ref, box_angles_ref = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz, box_lengths, box_angles = f.read(atom_indices=slice(None, None, 2))

    assert eq(xyz_ref[:, ::2, :], xyz)


def test_write_0(tmpdir, get_fn):
    fn_dcd = get_fn('frame0.dcd')
    fn = '{}/x.dcd'.format(tmpdir)
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz = f.read()[0]
    with DCDTrajectoryFile(fn, 'w') as f:
        f.write(xyz)
    with DCDTrajectoryFile(fn) as f:
        xyz2 = f.read()[0]

    eq(xyz, xyz2)


def test_write_1(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)

    with DCDTrajectoryFile(fn, 'w') as f:
        f.write(xyz)
    with DCDTrajectoryFile(fn) as f:
        xyz2 = f.read()[0]

    eq(xyz, xyz2)


def test_write_2(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((500, 3), dtype=np.float32)
    box_angles = 90 * np.ones((500, 3), dtype=np.float32)
    box_lengths[0, 0] = 10.0

    f = DCDTrajectoryFile(fn, 'w')
    f.write(xyz, box_lengths, box_angles)
    f.close()

    f = DCDTrajectoryFile(fn)
    xyz2, box_lengths2, box_angles2 = f.read()
    f.close()

    assert eq(xyz, xyz2)
    assert eq(box_lengths, box_lengths2)
    assert eq(box_angles, box_angles2)


def test_write_3(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((600, 3), dtype=np.float32)

    with DCDTrajectoryFile(fn, 'w') as f:
        with pytest.raises(ValueError):
            f.write(xyz, box_lengths)


def test_write_4(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((500, 3), dtype=np.float32)
    box_angles = 90 * np.ones((500, 3), dtype=np.float32)
    box_lengths[0, 0] = 10.0

    f = DCDTrajectoryFile(fn, 'w')
    for i in range(len(xyz)):
        f.write(xyz[i], box_lengths[i], box_angles[i])
    f.close()

    f = DCDTrajectoryFile(fn)
    xyz2, box_lengths2, box_angles2 = f.read()
    f.close()

    assert eq(xyz, xyz2)
    assert eq(box_lengths, box_lengths2)
    assert eq(box_angles, box_angles2)


def test_do_overwrite(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    with open(fn, 'w') as f:
        f.write('a')

    with DCDTrajectoryFile(fn, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10, 5, 3))


def test_dont_overwrite(tmpdir):
    fn = '{}/x.dcd'.format(tmpdir)
    with open(fn, 'w') as f:
        f.write('a')

    with pytest.raises(IOError):
        with DCDTrajectoryFile(fn, 'w', force_overwrite=False) as f:
            f.write(np.random.randn(10, 5, 3))


def test_read_closed(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    with pytest.raises(IOError):
        f = DCDTrajectoryFile(fn_dcd)
        f.close()
        f.read()


def test_write_closed(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    with pytest.raises(IOError):
        f = DCDTrajectoryFile(fn_dcd, 'w')
        f.close()
        f.write(np.random.randn(10, 3, 3))


def test_tell(get_fn):
    fn_dcd = get_fn('frame0.dcd')
    with DCDTrajectoryFile(get_fn('frame0.dcd')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek(get_fn):
    reference = DCDTrajectoryFile(get_fn('frame0.dcd')).read()[0]
    with DCDTrajectoryFile(get_fn('frame0.dcd')) as f:
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


def test_ragged_1(tmpdir):
    # try first writing no cell angles/lengths, and then adding some
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.random.randn(100, 5, 3)
    cell_lengths = np.random.randn(100, 3)
    cell_angles = np.random.randn(100, 3)

    with DCDTrajectoryFile(fn, 'w', force_overwrite=True) as f:
        f.write(xyz)
        with pytest.raises(ValueError):
            f.write(xyz, cell_lengths, cell_angles)


def test_ragged_2(tmpdir):
    # try first writing no cell angles/lengths, and then adding some
    fn = '{}/x.dcd'.format(tmpdir)
    xyz = np.random.randn(100, 5, 3)
    cell_lengths = np.random.randn(100, 3)
    cell_angles = np.random.randn(100, 3)

    # from mdtraj.formats import HDF5TrajectoryFile
    with DCDTrajectoryFile(fn, 'w', force_overwrite=True) as f:
        f.write(xyz, cell_lengths, cell_angles)
        with pytest.raises(ValueError):
            f.write(xyz)
