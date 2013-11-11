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


"""
Tests for the AMBER netcdf reader/writer code
"""

from mdtraj import netcdf, NetCDFTrajectoryFile
import os, tempfile
from nose.tools import assert_raises
import numpy as np
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, raises

TestDocstrings = DocStringFormatTester(netcdf, error_on_none=True)

fd, temp = tempfile.mkstemp(suffix='.nc')
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.close(fd)
    os.unlink(temp)


def test_read_after_close():
    f = NetCDFTrajectoryFile(get_fn('mdcrd.nc'))
    yield lambda: eq(f.n_atoms, 223)
    yield lambda: eq(f.n_frames, 101)

    f.close()

    # should be an ioerror if you read a file that's closed
    assert_raises(IOError, lambda: f.read())


def test_shape():
    xyz, time, boxlength, boxangles = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()

    yield lambda: eq(xyz.shape, (101, 223, 3))
    yield lambda: eq(time.shape, (101,))
    yield lambda: eq(boxlength, None)
    yield lambda: eq(boxangles, None)


def test_read_chunk_1():
    with NetCDFTrajectoryFile(get_fn('mdcrd.nc')) as f:
        a, b, c, d = f.read(10)
        e, f, g, h = f.read()

        yield lambda: eq(len(a), 10)
        yield lambda: eq(len(b), 10)

        yield lambda: eq(len(e), 101-10)
        yield lambda: eq(len(f), 101-10)

    xyz = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()[0]

    yield lambda: eq(a, xyz[0:10])
    yield lambda: eq(e, xyz[10:])


def test_read_chunk_2():
    with NetCDFTrajectoryFile(get_fn('mdcrd.nc')) as f:
        a, b, c, d = f.read(10)
        e, f, g, h = f.read(100000000000)

        yield lambda: eq(len(a), 10)
        yield lambda: eq(len(b), 10)

        yield lambda: eq(len(e), 101-10)
        yield lambda: eq(len(f), 101-10)

    xyz = NetCDFTrajectoryFile(get_fn('mdcrd.nc')).read()[0]

    yield lambda: eq(a, xyz[0:10])
    yield lambda: eq(e, xyz[10:])


def test_read_chunk_3():
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
        yield lambda: eq(a, xyz)
        yield lambda: eq(b, time)
        yield lambda: eq(c, boxlengths)
        yield lambda: eq(d, boxangles)


def test_read_write_2():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time)

    with NetCDFTrajectoryFile(temp) as f:
        a, b, c, d = f.read()
        yield lambda: eq(a, xyz)
        yield lambda: eq(b, time)
        yield lambda: eq(c.mask, np.ma.masked_all((100,3)).mask)
        yield lambda: eq(d.mask, np.ma.masked_all((100,3)).mask)


def test_read_write_25():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(xyz, time)
        f.write(xyz, time)

    with NetCDFTrajectoryFile(temp) as f:
        a, b, c, d = f.read()
        yield lambda: eq(a[0:100], xyz)
        yield lambda: eq(b[0:100], time)
        yield lambda: eq(c.mask[0:100], np.ma.masked_all((100,3)).mask)
        yield lambda: eq(d.mask[0:100], np.ma.masked_all((100,3)).mask)

        yield lambda: eq(a[100:], xyz)
        yield lambda: eq(b[100:], time)
        yield lambda: eq(c.mask[100:], np.ma.masked_all((100,3)).mask)
        yield lambda: eq(d.mask[100:], np.ma.masked_all((100,3)).mask)

def test_write_3():
    xyz = np.random.randn(100, 3, 3)
    time = np.random.randn(100)

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        # you can't supply cell_lengths without cell_angles
        assert_raises(ValueError, lambda: f.write(np.random.randn(100, 3, 3), cell_lengths=np.random.randn(100, 3)))
        # or the other way aroun
        assert_raises(ValueError, lambda: f.write(np.random.randn(100, 3, 3), cell_angles=np.random.randn(100, 3)))


def test_n_atoms():
    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(1,11,3))
    with NetCDFTrajectoryFile(temp) as f:
        eq(f.n_atoms, 11)


def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10,5,3))

@raises(IOError)
def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with NetCDFTrajectoryFile(temp, 'w', force_overwrite=False) as f:
        f.write(np.random.randn(10,5,3))


def test_trajectory_save_load():
    t = md.load(get_fn('native.pdb'))
    t.unitcell_lengths = 1 * np.ones((1, 3))
    t.unitcell_angles = 90 * np.ones((1, 3))

    t.save(temp)
    t2 = md.load(temp, top=t.topology)

    eq(t.xyz, t2.xyz)
    eq(t.unitcell_lengths, t2.unitcell_lengths)
