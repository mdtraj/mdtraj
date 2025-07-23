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
Code to test the mdconvert script. These tests take about two minutes to run.

This checks all pairs for formats, converting from format x -> format y. it
also tries using striding to subsample the trajectory and atom_indices, so it
does significant integration testing of the XXXTrajectoryFile modules as well.
"""

import os
import subprocess
import sys

import numpy as np
import pytest

import mdtraj as md
from mdtraj.testing import eq

on_win = sys.platform == "win32"
on_py3 = sys.version_info >= (3, 0)


def test_index(h5traj):
    # Check that extracting a specific index works
    traj, in_fn, tmpdir = h5traj
    out_fn = f"{tmpdir}/frame4.pdb"
    subprocess.check_call(["mdconvert", in_fn, "-i", "4", "-o", out_fn])
    frame4 = md.load(out_fn)
    eq(frame4.xyz, traj[4].xyz)


def test_slice(h5traj):
    # Check that extracting a specific slice works
    traj, in_fn, tmpdir = h5traj
    out_fn = f"{tmpdir}/frame13.pdb"
    subprocess.check_call(["mdconvert", in_fn, "-i", "1:5:2", "-o", out_fn])
    frame13 = md.load(out_fn)
    eq(frame13.xyz, traj[1:5:2].xyz)


extensions = [
    "xtc",
    "dcd",
    "trr",
    "nc",
    "pdb",
    "h5",
    "lh5",
    "netcdf",
]


@pytest.fixture(params=extensions, ids=lambda x: "from-" + x)
def extension(request):
    if on_win and request.param == "lh5":
        pytest.skip("No lh5 on windows py3")
    return request.param


def test_pairwise(h5traj, extension, monkeypatch):
    """ensure that the xyz coordinates are preserved by a trip
    from python -> save in format X -> mdconvert to format Y -> python
    """

    def test_base(h5traj, extension, monkeypatch):
        traj, _, tmpdir = h5traj
        ext1 = extension

        # save one copy of traj for use as a topology file
        topology_fn = f"{tmpdir}/topology.pdb"
        traj[0].save(topology_fn)

        # save a .dat file for the atom_indices so that we can test
        # mdconvert's atom_indices flag
        atom_indices = np.array([0, 3])
        atom_indices_fn = f"{tmpdir}/atom_indices.dat"
        np.savetxt(atom_indices_fn, atom_indices, fmt="%d")

        in_fn = f"{tmpdir}/traj.{ext1}"
        traj.save(in_fn)
        working_dir = f"{tmpdir}/from-{ext1}"
        os.mkdir(working_dir)

        for ext2 in extensions:
            out_fn = f"traj.{ext2}"

            command1 = ["mdconvert", in_fn, "-o", out_fn, "-c 6"]
            if ext2 in ["pdb", "h5", "lh5"]:
                # if we're saving a pdb or h5, we need to give it a topology too
                command1 += ["-t", topology_fn]

            # TODO: test fixture
            subprocess.check_call(command1, cwd=working_dir)

            # Use the --atom_indices flag to mdconvert
            command2 = command1 + ["-a", atom_indices_fn]
            command2[3] = "subset." + out_fn  # make sure the output goes to a different file
            subprocess.check_call(command2, cwd=working_dir)

            # Use the --stride 3 flag
            command3 = command1 + ["-s 3"]
            command3[3] = "stride." + out_fn  # change the out filename, so they don't clobbed
            subprocess.check_call(command3, cwd=working_dir)

            # ensure that the xyz coordinates are preserved by a trip
            # from python -> save in format X -> mdconvert to format Y -> python
            load_kwargs_check1 = {}
            load_kwargs_check2 = {}
            if ext2 not in ["pdb", "h5", "lh5"]:
                load_kwargs_check1["top"] = traj.topology
                load_kwargs_check2["top"] = traj.topology.subset(atom_indices)

            out1 = md.load(os.path.join(working_dir, out_fn), **load_kwargs_check1)
            out2 = md.load(
                os.path.join(working_dir, "subset." + out_fn),
                **load_kwargs_check2,
            )
            out3 = md.load(
                os.path.join(working_dir, "stride." + out_fn),
                **load_kwargs_check1,
            )

            if ext1 in ["lh5"] or ext2 in ["lh5"]:
                decimal = 3
            else:
                decimal = 6
            eq(out1.xyz, traj.xyz, decimal=decimal)
            eq(out2.xyz, traj.xyz[:, atom_indices], decimal=decimal)
            eq(out3.xyz, traj.xyz[::3], decimal=decimal)

            if ext1 not in ["lh5"] and ext2 not in ["lh5"]:
                # binpos doesn't save unitcell information
                eq(out1.unitcell_vectors, traj.unitcell_vectors, decimal=2)
                eq(out2.unitcell_vectors, traj.unitcell_vectors, decimal=2)
                eq(out3.unitcell_vectors, traj.unitcell_vectors[::3], decimal=2)

            if all(e in ["xtc", "trr", "nc", "h5"] for e in [ext1, ext2]):
                # these formats contain time information
                if all(e in ["nc"] for e in [ext1, ext2]):
                    with monkeypatch.context() as m:
                        m.setitem(sys.modules, "netCDF4", None)
                        eq(out1.time, traj.time)
                        eq(out2.time, traj.time)
                        eq(out3.time, traj.time[::3])

                eq(out1.time, traj.time)
                eq(out2.time, traj.time)
                eq(out3.time, traj.time[::3])

            if ext2 in ["pdb", "h5", "lh5"]:
                # these formats contain a topology in the file that was
                # read from disk
                eq(out1.topology, traj.topology)
                eq(out2.topology, traj.topology.subset(atom_indices))
                eq(out3.topology, traj.topology)

    if extension in ("nc"):
        # If extension is nc, we need to test with both SciPy and NetCDF4.
        # First, we use pytest.monkeypatch to remove the netCDF4 import from the environment
        # Then, we make sure the extension for both SciPy and NetCDF4 tests are different (or
        # there will be SameFileError/FileExistsError failures)
        # All these changes will be reverted outside the context manager, and the netCDF4 test will run again
        with monkeypatch.context() as m:
            m.setitem(sys.modules, "netCDF4", None)
            test_base(h5traj, ".scipy.nc", monkeypatch)

    # For testing most formats and with netCDF4 (if format is nc)
    test_base(h5traj, extension, monkeypatch)
