##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
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
import tempfile
import os
import sys
import mdtraj as md
from mdtraj.formats import LH5TrajectoryFile
from mdtraj.testing import eq
import pytest

on_win = (sys.platform == 'win32')
on_py3 = (sys.version_info >= (3, 0))

# special pytest global to mark all tests in this module
pytestmark = pytest.mark.skipif(on_win and on_py3, reason='lh5 not supported on windows on python 3')

fd, temp = tempfile.mkstemp(suffix='.lh5')


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.close(fd)
    os.unlink(temp)


def test_write_coordinates():
    coordinates = np.random.randn(4, 10, 3)
    with LH5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(), coordinates, decimal=3)

    with LH5TrajectoryFile(temp) as f:
        f.seek(2)
        eq(f.read(), coordinates[2:], decimal=3)
        f.seek(0)
        eq(f.read(), coordinates[0:], decimal=3)
        f.seek(-1, 2)
        eq(f.read(), coordinates[3:], decimal=3)


def test_write_coordinates_reshape():
    coordinates = np.random.randn(10, 3)
    with LH5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(), coordinates.reshape(1, 10, 3), decimal=3)


def test_write_multiple():
    coordinates = np.random.randn(4, 10, 3)
    with LH5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)
        f.write(coordinates)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(), np.vstack((coordinates, coordinates)), decimal=3)


def test_topology(get_fn):
    top = md.load(get_fn('native.pdb')).topology

    with LH5TrajectoryFile(temp, 'w') as f:
        f.topology = top

    with LH5TrajectoryFile(temp) as f:
        assert f.topology == top


def test_read_slice_0():
    coordinates = np.random.randn(4, 10, 3)
    with LH5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(n_frames=2), coordinates[:2], decimal=3)
        eq(f.read(n_frames=2), coordinates[2:4], decimal=3)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(stride=2), coordinates[::2], decimal=3)

    with LH5TrajectoryFile(temp) as f:
        eq(f.read(stride=2, atom_indices=np.array([0, 1])), coordinates[::2, [0, 1], :], decimal=3)


def test_vsite_elements(get_fn):
    #  Test case for issue #263
    pdb_filename = get_fn('GG-tip4pew.pdb')
    trj = md.load(pdb_filename)
    trj.save_lh5(temp)

    trj2 = md.load(temp, top=pdb_filename)


def test_do_overwrite_0():
    with open(temp, 'w') as f:
        f.write('a')

    with LH5TrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10, 5, 3))


def test_do_overwrite_1():
    with open(temp, 'w') as f:
        f.write('a')

    with pytest.raises(IOError):
        with LH5TrajectoryFile(temp, 'w', force_overwrite=False) as f:
            f.write(np.random.randn(10, 5, 3))
