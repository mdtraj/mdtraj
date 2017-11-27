##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Teng Lin
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
from mdtraj.formats import DTRTrajectoryFile, DCDTrajectoryFile
from mdtraj.testing import eq
import pytest

def test_read(get_fn):
    # Test the default read and compare against reference trajectory in dcd format
    fn_dtr = get_fn('frame0.dtr')
    fn_dcd = get_fn('frame0.dcd')
    dtr_traj = DTRTrajectoryFile(fn_dtr)
    eq(len(dtr_traj), 501)
    xyz, times, cell_lens, cell_angles = dtr_traj.read()
    xyz2, cell_lens2, cell_angles2 = DCDTrajectoryFile(fn_dcd).read()
    eq(xyz, xyz2)
    eq(cell_lens, cell_lens2)
    eq(cell_angles, cell_angles2)


def test_read_1(get_fn):
    # test read with n_frame
    fn_dtr = get_fn('frame0.dtr')
    xyz, times, cell_lens, cell_angles = DTRTrajectoryFile(fn_dtr).read()
    xyz2, times2, cell_lens2, cell_angles2 = DTRTrajectoryFile(fn_dtr).read(n_frames=501)
    eq(xyz, xyz2)
    eq(times, times2)
    eq(cell_lens, cell_lens2)
    eq(cell_angles, cell_angles2)


def test_read_2(get_fn):
    # test read with atom indices
    fn_dtr = get_fn('frame0.dtr')
    indices = np.array([0, 3, 12, 4])
    xyz, times, cell_lens, cell_angles = DTRTrajectoryFile(fn_dtr).read()
    xyz2, times2, cell_lens2, cell_angles2 = DTRTrajectoryFile(fn_dtr).read(atom_indices=indices)
    eq(xyz[:, indices, :], xyz2)
    eq(times, times2)
    eq(cell_lens, cell_lens2)
    eq(cell_angles, cell_angles2)


def test_read_3(get_fn):
    # Test read with n_frames
    fn_dtr = get_fn('frame0.dtr')
    dtr_traj = DTRTrajectoryFile(fn_dtr)
    dtr_traj.seek(1)
    xyz, times, cell_lens, cell_angles = dtr_traj.read(n_frames=900)
    eq(len(xyz), 500)


def test_read_stride(get_fn):
    # Read dtr with stride
    fn_dtr = get_fn('frame0.dtr')
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz1, times1, box_lengths1, box_angles1 = f.read()
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz2, times2, box_lengths2, box_angles2 = f.read(stride=2)

    assert eq(xyz1[::2], xyz2)
    assert eq(times1[::2], times2)
    assert eq(box_lengths1[::2], box_lengths2)
    assert eq(box_angles1[::2], box_angles2)


def test_read_4(get_fn):
    # Read dtr with stride and n_frames
    fn_dtr = get_fn('frame0.dtr')
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz1, times1, box_lengths1, box_angles1 = f.read()
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz2, times2, box_lengths2, box_angles2 = f.read(n_frames=300, stride=2)

    assert eq(xyz1[::2], xyz2)
    assert eq(times1[::2], times2)
    assert eq(box_lengths1[::2], box_lengths2)
    assert eq(box_angles1[::2], box_angles2)


def test_read_5(get_fn):
    # Check streaming read of frames 1 at a time
    fn_dtr = get_fn('frame0.dtr')
    xyz_ref, times_ref, box_lengths_ref, box_angles_ref = DTRTrajectoryFile(fn_dtr).read()

    reader = DTRTrajectoryFile(fn_dtr)
    for i in range(len(xyz_ref)):
        xyz, times, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(times_ref[np.newaxis, i], times)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)


def test_read_6(get_fn):
    # Check streaming read followed by reading the 'rest'
    fn_dtr = get_fn('frame0.dtr')
    xyz_ref, times_ref, box_lengths_ref, box_angles_ref = DTRTrajectoryFile(fn_dtr).read()

    reader = DTRTrajectoryFile(fn_dtr)
    for i in range(len(xyz_ref) // 2):
        xyz, times, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(times_ref[np.newaxis, i], times)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)

    xyz_rest, times_rest, box_rest, angles_rest = reader.read()
    i = len(xyz_ref) // 2
    assert eq(xyz_ref[i:], xyz_rest)
    assert eq(times_ref[i:], times_rest)
    assert eq(box_lengths_ref[i:], box_rest)
    assert eq(box_angles_ref[i:], angles_rest)

    assert len(xyz_ref) == i + len(xyz_rest)


def test_read_7(get_fn):
    # Test two full reads
    fn_dtr = get_fn('frame0.dtr')
    reader = DTRTrajectoryFile(fn_dtr)
    xyz, times, cell_lens, cell_angles = reader.read()
    xyz, times, cell_lens, cell_angles = reader.read()
    eq(len(xyz), 0)
    eq(len(times), 0)
    eq(len(cell_lens), 0)
    eq(len(cell_angles), 0)


def test_read_8(get_fn):
    fn_dtr = get_fn('frame0.dtr')
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz_ref, times_ref, box_lengths_ref, box_angles_ref = f.read()
    with DTRTrajectoryFile(fn_dtr) as f:
        xyz, times, box_lengths, box_angles = f.read(atom_indices=slice(None, None, 2))

    assert eq(xyz_ref[:, ::2, :], xyz)


def test_write_1(tmpdir, get_fn):
    # Test write
    fn_dtr = get_fn('frame0.dtr')
    fn = '{}/x.dtr'.format(tmpdir)
    xyz, times, cell_lens, cell_angles = DTRTrajectoryFile(fn_dtr).read()
    xyz += 1
    DTRTrajectoryFile(fn, 'w').write(xyz, cell_lengths=cell_lens,
                                     cell_angles=cell_angles, times=times)

    xyz2, times2, cell_lens2, cell_angles2 = DTRTrajectoryFile(fn).read()

    eq(xyz, xyz2)
    eq(times, times2)
    eq(cell_lens, cell_lens2)
    eq(cell_angles, cell_angles2)


def test_write_2(tmpdir, get_fn):
    # Test two separate write call
    fn_dtr = get_fn('frame0.dtr')
    fn = '{}/x.dtr'.format(tmpdir)
    xyz, times, cell_lens, cell_angles = DTRTrajectoryFile(fn_dtr).read()
    writer = DTRTrajectoryFile(fn, 'w')
    writer.write(xyz, cell_lengths=cell_lens,
                 cell_angles=cell_angles, times=times)

    n_frames = len(xyz)
    times += 50.0
    writer.write(xyz, cell_lengths=cell_lens,
                 cell_angles=cell_angles, times=times)

    writer.close()

    xyz2, times2, cell_lens2, cell_angles2 = DTRTrajectoryFile(fn).read()

    eq(len(xyz2), n_frames * 2)
    eq(xyz, xyz2[n_frames:])
    eq(times, times2[n_frames:])
    eq(cell_lens, cell_lens2[n_frames:])
    eq(cell_angles, cell_angles2[n_frames:])


def test_write_3(tmpdir):
    # Test a random write operation
    xyz = np.array(np.random.uniform(low=-50, high=-50, size=(3, 17, 3)), dtype=np.float32)
    times = np.array([1, 23.0, 48.0], dtype=np.float64)
    cell_lengths = np.array(np.random.uniform(low=100, high=200, size=(3, 3)), dtype=np.float32)
    cell_angles = np.array([[90, 90, 90],
                            [80, 100, 120],
                            [120, 90, 80]],
                           dtype=np.float32)

    fn = '{}/x.dtr'.format(tmpdir)
    with DTRTrajectoryFile(fn, 'w') as f:
        f.write(xyz, cell_lengths=cell_lengths,
                cell_angles=cell_angles, times=times)
    with DTRTrajectoryFile(fn) as f:
        xyz2, times2, cell_lengths2, cell_angles2 = f.read()

    eq(xyz, xyz2)


def test_write_4(tmpdir):
    # Test write error
    xyz = np.array(np.random.uniform(low=-50, high=-50, size=(3, 17, 3)), dtype=np.float32)
    times = np.array([1, 23.0, 48.0], dtype=np.float64)
    cell_lengths = np.array(np.random.uniform(low=100, high=200, size=(3, 3)), dtype=np.float32)
    cell_angles = np.array([[90, 90, 90],
                            [80, 100, 120],
                            [120, 90, 80]],
                           dtype=np.float32)

    bad_times = np.array([21, 3.0, 48.0], dtype=np.float64)
    fn = '{}/x.dtr'.format(tmpdir)
    f = DTRTrajectoryFile(fn, 'w')
    with pytest.raises(ValueError):
        f.write(xyz, cell_lengths=cell_lengths)
    with pytest.raises(ValueError):
        f.write(xyz, cell_angles=cell_angles)
    with pytest.raises(ValueError):
        f.write(xyz, times=times)
    with pytest.raises(ValueError):
        f.write(xyz, cell_lengths=cell_lengths,
                cell_angles=cell_angles, times=bad_times)
    f.close()


def test_seek(get_fn):
    fn_dtr = get_fn('frame0.dtr')
    reference = DTRTrajectoryFile(fn_dtr).read()[0]
    with DTRTrajectoryFile(fn_dtr) as f:
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


def test_read_closed(get_fn):
    fn_dtr = get_fn('frame0.dtr')
    with pytest.raises(IOError):
        f = DTRTrajectoryFile(fn_dtr)
        f.close()
        f.read()


def test_tell(get_fn):
    fn_dtr = get_fn('frame0.dtr')
    with DTRTrajectoryFile(fn_dtr) as f:
        last = len(f)
        eq(f.tell(), 0)

        f.read(2)
        eq(f.tell(), 2)

        f.read(100)
        eq(f.tell(), 102)

        f.seek(600)
        eq(f.tell(), last)

