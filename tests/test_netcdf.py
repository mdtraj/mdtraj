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

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest
from _pytest.monkeypatch import MonkeyPatch

import mdtraj as md
from mdtraj.formats import NetCDFTrajectoryFile
from mdtraj.testing import eq

needs_cpptraj = pytest.mark.skipif(
    shutil.which("cpptraj") is None,
    reason="This test requires cpptraj from AmberTools to be installed (http://ambermd.org)",
)

fd, temp = tempfile.mkstemp(suffix=".nc")
fd2, temp2 = tempfile.mkstemp(suffix=".nc")


class TestNetCDFNetCDF4:
    """
    This class contains all the tests that we would also want to run with scipy.
    Now a class so we can subclass it for later.
    """

    def teardown_module(self, module):
        """remove the temporary file created by tests in this file
        this gets automatically called by pytest"""
        os.close(fd)
        os.close(fd2)
        os.unlink(temp)
        os.unlink(temp2)

    def test_read_after_close(self, get_fn):
        """Default test using netCDF4"""
        f = NetCDFTrajectoryFile(get_fn("mdcrd.nc"))
        assert eq(f.n_atoms, 223)
        assert eq(f.n_frames, 101)

        f.close()

        # should be an IOError if you read a file that's closed
        with pytest.raises(IOError):
            f.read()

    def test_shape(self, get_fn):
        """Default test using netCDF4"""
        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as netcdf_file:
            xyz, time, boxlength, boxangles = netcdf_file.read()

        assert eq(xyz.shape, (101, 223, 3))
        assert eq(time.shape, (101,))
        assert eq(boxlength, None)
        assert eq(boxangles, None)

    def test_read_chunk_1(self, get_fn):
        """Default test using netCDF4"""
        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            a, b, c, d = file.read(10)
            e, f, g, h = file.read()

            assert eq(len(a), 10)
            assert eq(len(b), 10)

            assert eq(len(e), 101 - 10)
            assert eq(len(f), 101 - 10)

        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            xyz = file.read()[0]

        assert eq(a, xyz[0:10])
        assert eq(e, xyz[10:])

    def test_read_chunk_2(self, get_fn):
        """Default test using netCDF4"""

        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            a, b, c, d = file.read(10)
            e, f, g, h = file.read(100000000000)

            assert eq(len(a), 10)
            assert eq(len(b), 10)

            assert eq(len(e), 101 - 10)
            assert eq(len(f), 101 - 10)

        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            xyz = file.read()[0]

        assert eq(a, xyz[0:10])
        assert eq(e, xyz[10:])

    def test_read_chunk_3(self, get_fn):
        """Default test using netCDF4"""
        # too big of a chunk should not be an issue
        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            a = file.read(1000000000)
        with NetCDFTrajectoryFile(get_fn("mdcrd.nc")) as file:
            b = file.read()

        eq(a[0], b[0])

    def test_read_write_1(self):
        """Default test using netCDF4"""
        xyz = np.random.randn(100, 3, 3)
        time = np.random.randn(100)
        boxlengths = np.random.randn(100, 3)
        boxangles = np.random.randn(100, 3)

        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(xyz, time, boxlengths, boxangles)

        with NetCDFTrajectoryFile(temp) as f:
            a, b, c, d = f.read()
            assert eq(a, xyz)
            assert eq(b, time)
            assert eq(c, boxlengths)
            assert eq(d, boxangles)

    def test_read_write_2(self, get_fn):
        """Default test using netCDF4"""
        xyz = np.random.randn(5, 22, 3)
        time = np.random.randn(5)

        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(xyz, time)

        with NetCDFTrajectoryFile(temp) as f:
            rcoord, rtime, rlengths, rangles = f.read()
            assert eq(rcoord, xyz)
            assert eq(rtime, time)
            assert eq(rlengths, None)
            assert eq(rangles, None)

        t = md.load(temp, top=get_fn("native.pdb"))
        eq(t.unitcell_angles, None)
        eq(t.unitcell_lengths, None)

    def test_ragged_1(self):
        """Default test using netCDF4"""
        # try first writing no cell angles/lengths, and then adding some
        xyz = np.random.randn(100, 3, 3)
        time = np.random.randn(100)
        cell_lengths = np.random.randn(100, 3)
        cell_angles = np.random.randn(100, 3)

        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(xyz, time)
            with pytest.raises(ValueError):
                f.write(xyz, time, cell_lengths, cell_angles)

    def test_ragged_2(self):
        """Default test using netCDF4"""
        # try first writing no cell angles/lengths, and then adding some
        xyz = np.random.randn(100, 3, 3)
        time = np.random.randn(100)
        cell_lengths = np.random.randn(100, 3)
        cell_angles = np.random.randn(100, 3)

        # from mdtraj.formats import HDF5TrajectoryFile
        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(xyz, time, cell_lengths, cell_angles)
            with pytest.raises(ValueError):
                f.write(xyz, time)

    def test_read_write_25(self):
        """Default test using netCDF4"""
        xyz = np.random.randn(100, 3, 3)
        time = np.random.randn(100)

        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
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

    def test_write_3(self):
        """Default test using netCDF4"""
        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            # you can't supply cell_lengths without cell_angles
            with pytest.raises(ValueError):
                f.write(np.random.randn(100, 3, 3), cell_lengths=np.random.randn(100, 3))
            # or the other way around
            with pytest.raises(ValueError):
                f.write(np.random.randn(100, 3, 3), cell_angles=np.random.randn(100, 3))

    def test_n_atoms(self):
        """Default test using netCDF4"""
        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(np.random.randn(1, 11, 3))
        with NetCDFTrajectoryFile(temp) as f:
            eq(f.n_atoms, 11)

    def test_do_overwrite(self):
        """Default test using netCDF4"""
        with open(temp, "w") as f:
            f.write("a")

        with NetCDFTrajectoryFile(temp, "w", force_overwrite=True) as f:
            f.write(np.random.randn(10, 5, 3))

    def test_do_not_overwrite(self):
        """Default test using netCDF4"""
        with open(temp, "w") as f:
            f.write("a")

        with pytest.raises(IOError):
            with NetCDFTrajectoryFile(temp, "w", force_overwrite=False) as f:
                f.write(np.random.randn(10, 5, 3))

    def test_trajectory_save_load(self, get_fn):
        """Default test using netCDF4"""
        t = md.load(get_fn("native.pdb"))
        t.unitcell_lengths = 1 * np.ones((1, 3))
        t.unitcell_angles = 90 * np.ones((1, 3))

        t.save(temp)
        t2 = md.load(temp, top=t.topology)

        eq(t.xyz, t2.xyz)
        eq(t.unitcell_lengths, t2.unitcell_lengths)


class TestNetCDFScipy(TestNetCDFNetCDF4):
    """This inherits the TestNetCDFNetCDF4 class and run all tests with SciPy"""

    def setup_method(self, method):
        """Patching out netCDF4. This is the way to do it inside a class"""
        monkeypatch = MonkeyPatch()
        monkeypatch.setitem(sys.modules, "netCDF4", None)

    def teardown_method(self, method):
        """Undoing most changes, just in case."""
        monkeypatch = MonkeyPatch()
        monkeypatch.delitem(sys.modules, "netCDF4", None)


@needs_cpptraj
def test_cpptraj(get_fn):
    trj0 = md.load(get_fn("frame0.dcd"), top=get_fn("frame0.pdb"))
    trj0.save(temp)

    top = get_fn("frame0.pdb")
    subprocess.check_call(
        [
            "cpptraj",
            "-p",
            top,
            "-y",
            temp,
            "-x",
            temp2,
        ],
    )

    trj1 = md.load(temp, top=top)
    trj2 = md.load(temp2, top=top)

    np.testing.assert_array_almost_equal(trj0.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(
        trj0.unitcell_vectors,
        trj2.unitcell_vectors,
    )
    np.testing.assert_array_almost_equal(
        trj1.unitcell_vectors,
        trj2.unitcell_vectors,
    )

    np.testing.assert_array_almost_equal(trj0.time, trj1.time)
    np.testing.assert_array_almost_equal(trj0.time, trj2.time)
    np.testing.assert_array_almost_equal(trj1.time, trj2.time)
