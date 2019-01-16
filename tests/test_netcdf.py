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

from mdtraj.formats import NetCDFTrajectoryFile
import os, tempfile
import numpy as np
import mdtraj as md
import subprocess
from mdtraj.testing import eq
import pytest

from distutils.spawn import find_executable

needs_cpptraj = pytest.mark.skipif(
    find_executable('cpptraj') is None,
    reason="This test requires cpptraj from AmberTools to be installed "
           "(http://ambermd.org)"
)

fd, temp = tempfile.mkstemp(suffix='.nc')
fd2, temp2 = tempfile.mkstemp(suffix='.nc')


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.close(fd)
    os.close(fd2)
    os.unlink(temp)
    os.unlink(temp2)


def test_read_after_close(get_fn):
    f = NetCDFTrajectoryFile(get_fn('mdcrd.nc'))
    assert eq(f.n_atoms, 223)
    assert eq(f.n_frames, 101)

    f.close()

    # should be an IOError if you read a file that's closed
    with pytest.raises(IOError):
        f.read()


def test_shape(get_fn):
    xyz, time, boxlength, boxangles = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()

    assert eq(xyz.shape, (101, 223, 3))
    assert eq(time.shape, (101,))
    assert eq(boxlength, None)
    assert eq(boxangles, None)


def test_read_chunk_1(get_fn):
    with NetCDFTrajectoryFile(get_fn('mdcrd.nc')) as f:
        a, b, c, d = f.read(10)
        e, f, g, h = f.read()

        assert eq(len(a), 10)
        assert eq(len(b), 10)

        assert eq(len(e), 101 - 10)
        assert eq(len(f), 101 - 10)

    xyz = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()[0]

    assert eq(a, xyz[0:10])
    assert eq(e, xyz[10:])


def test_read_chunk_2(get_fn):
    with NetCDFTrajectoryFile(get_fn('mdcrd.nc')) as f:
        a, b, c, d = f.read(10)
        e, f, g, h = f.read(100000000000)

        assert eq(len(a), 10)
        assert eq(len(b), 10)

        assert eq(len(e), 101 - 10)
        assert eq(len(f), 101 - 10)

    xyz = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()[0]

    assert eq(a, xyz[0:10])
    assert eq(e, xyz[10:])


def test_read_chunk_3(get_fn):
    # too big of a chunk should not be an issue
    a = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read(1000000000)
    b = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()

    eq(a[0], b[0])


def test_read_write_1():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)
    boxlengths = np.random.randn(100, 3)
    boxangles = np.random.randn(100, 3)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time, boxlengths, boxangles)

    with NetCDFTrajectoryFile(temp) as f:
        a, b, c, d = f.read()
        assert eq(a, xyz)
        assert eq(b, time)
        assert eq(c, boxlengths)
        assert eq(d, boxangles)


def test_read_write_2(get_fn):
    xyz = np.random.randn(5, 22, 3)
    time = np.random.randn(5)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time)

    with NetCDFTrajectoryFile(temp) as f:
        rcoord, rtime, rlengths, rangles = f.read()
        assert eq(rcoord, xyz)
        assert eq(rtime, time)
        assert eq(rlengths, None)
        assert eq(rangles, None)

    t = md.load(temp, top=get_fn('native.pdb'))
    eq(t.unitcell_angles, None)
    eq(t.unitcell_lengths, None)


def test_ragged_1():
    # try first writing no cell angles/lengths, and then adding some
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)
    cell_lengths = np.random.randn(100, 3)
    cell_angles = np.random.randn(100, 3)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time)
        with pytest.raises(ValueError):
            f.write(xyz, time, cell_lengths, cell_angles)


def test_ragged_2():
    # try first writing no cell angles/lengths, and then adding some
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)
    cell_lengths = np.random.randn(100, 3)
    cell_angles = np.random.randn(100, 3)

    # from mdtraj.formats import HDF5TrajectoryFile
    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time, cell_lengths, cell_angles)
        with pytest.raises(ValueError):
            f.write(xyz, time)


def test_read_write_25():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time)
        f.write(xyz, time)

    with NetCDFTrajectoryFile(temp) as f:
        a, b, c, d = f.read()
        assert eq(a[0:100], xyz)
        assert eq(b[0:100], time)
        assert eq(c, None)
        assert eq(d, None)

        assert eq(a[100:], xyz)
        assert eq(b[100:], time)
        assert eq(c, None)
        assert eq(d, None)


def test_write_3():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        # you can't supply cell_lengths without cell_angles
        with pytest.raises(ValueError):
            f.write(np.random.randn(100, 3, 3), cell_lengths=np.random.randn(100, 3))
        # or the other way around
        with pytest.raises(ValueError):
            f.write(np.random.randn(100, 3, 3), cell_angles=np.random.randn(100, 3))


def test_n_atoms():
    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(1, 11, 3))
    with NetCDFTrajectoryFile(temp) as f:
        eq(f.n_atoms, 11)


def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10, 5, 3))


def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with pytest.raises(IOError):
        with NetCDFTrajectoryFile(temp, 'w', force_overwrite=False) as f:
            f.write(np.random.randn(10, 5, 3))


def test_trajectory_save_load(get_fn):
    t = md.load(get_fn('native.pdb'))
    t.unitcell_lengths = 1 * np.ones((1, 3))
    t.unitcell_angles = 90 * np.ones((1, 3))

    t.save(temp)
    t2 = md.load(temp, top=t.topology)

    eq(t.xyz, t2.xyz)
    eq(t.unitcell_lengths, t2.unitcell_lengths)


@needs_cpptraj
def test_cpptraj(get_fn):
    trj0 = md.load(get_fn('frame0.dcd'), top=get_fn('frame0.pdb'))
    trj0.save(temp)

    top = get_fn("frame0.pdb")
    subprocess.check_call([
        'cpptraj',
        '-p', top,
        '-y', temp,
        '-x', temp2
    ])

    trj1 = md.load(temp, top=top)
    trj2 = md.load(temp2, top=top)

    np.testing.assert_array_almost_equal(trj0.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj0.unitcell_vectors,
                                         trj2.unitcell_vectors)
    np.testing.assert_array_almost_equal(trj1.unitcell_vectors,
                                         trj2.unitcell_vectors)

    np.testing.assert_array_almost_equal(trj0.time, trj1.time)
    np.testing.assert_array_almost_equal(trj0.time, trj2.time)
    np.testing.assert_array_almost_equal(trj1.time, trj2.time)
