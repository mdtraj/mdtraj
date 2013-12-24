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
Test the cython dcd module

Note, this file cannot be located in the dcd subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import tempfile, os
import numpy as np
from mdtraj import DCDTrajectoryFile, io
from mdtraj.testing import get_fn, eq, DocStringFormatTester, raises

#TestDocstrings = DocStringFormatTester(dcd, error_on_none=True)

fn_dcd = get_fn('frame0.dcd')
pdb = get_fn('native.pdb')

fd, temp = tempfile.mkstemp(suffix='.dcd')
os.close(fd)
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_read():
    xyz, box_lengths, box_angles = DCDTrajectoryFile(fn_dcd).read()
    xyz2 = io.loadh(get_fn('frame0.dcd.h5'), 'xyz')

    eq(xyz, xyz2)


def test_read_2():
    "DCDReader: check nframes"
    xyz1, box_lengths1, box_angles1 = DCDTrajectoryFile(fn_dcd).read()
    xyz2, box_lengths2, box_angles2 = DCDTrajectoryFile(fn_dcd).read(10000)

    yield lambda: eq(xyz1, xyz2)
    yield lambda: eq(box_lengths1, box_lengths2)
    yield lambda: eq(box_angles1, box_angles2)


def test_read_stride():
    "Read dcd with stride"
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz1, box_lengths1, box_angles1 = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2, box_lengths2, box_angles2 = f.read(stride=2)

    yield lambda: eq(xyz1[::2], xyz2)
    yield lambda: eq(box_lengths1[::2], box_lengths2)
    yield lambda: eq(box_angles1[::2], box_angles2)


def test_read_stride_2():
    "Read dcd with stride when n_frames is supplied (different code path)"
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz1, box_lengths1, box_angles1 = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2, box_lengths2, box_angles2 = f.read(n_frames=1000, stride=2)

    yield lambda: eq(xyz1[::2], xyz2)
    yield lambda: eq(box_lengths1[::2], box_lengths2)
    yield lambda: eq(box_angles1[::2], box_angles2)


def test_read_3():
    "DCDReader: check streaming read of frames 1 at a time"
    xyz_ref, box_lengths_ref, box_angles_ref = DCDTrajectoryFile(fn_dcd).read()

    reader = DCDTrajectoryFile(fn_dcd)
    for i in range(len(xyz_ref)):
        xyz, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)


def test_read_4():
    "DCDReader: check streaming read followed by reading the 'rest'"
    xyz_ref, box_lengths_ref, box_angles_ref = DCDTrajectoryFile(fn_dcd).read()

    reader = DCDTrajectoryFile(fn_dcd)
    for i in range(int(len(xyz_ref)/2)):
        xyz, box_lenths, box_angles = reader.read(1)
        eq(xyz_ref[np.newaxis, i], xyz)
        eq(box_lengths_ref[np.newaxis, i], box_lenths)
        eq(box_angles_ref[np.newaxis, i], box_angles)

    xyz_rest, box_rest, angles_rest = reader.read()
    yield lambda: eq(xyz_ref[i+1:], xyz_rest)
    yield lambda: eq(box_lengths_ref[i+1:], box_rest)
    yield lambda: eq(box_angles_ref[i+1:], angles_rest)

    yield lambda: len(xyz_ref) == i + len(xyz_rest)


def test_read_5():
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz_ref, box_lengths_ref, box_angles_ref = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz, box_lengths, box_angles = f.read(atom_indices=[1,2,5])

    yield lambda: eq(xyz_ref[:, [1,2,5], :], xyz)


def test_read_6():
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz_ref, box_lengths_ref, box_angles_ref = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz, box_lengths, box_angles = f.read(atom_indices=slice(None, None, 2))

    yield lambda: eq(xyz_ref[:, ::2, :], xyz)


def test_write_0():
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz = f.read()[0]
    with DCDTrajectoryFile(temp, 'w') as f:
        f.write(xyz)
    with DCDTrajectoryFile(temp) as f:
        xyz2 = f.read()[0]

    eq(xyz, xyz2)


def test_write_1():
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)

    with DCDTrajectoryFile(temp, 'w') as f:
        f.write(xyz)
    with DCDTrajectoryFile(temp) as f:
        xyz2 = f.read()[0]

    eq(xyz, xyz2)


def test_write_2():
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((500, 3), dtype=np.float32)
    box_angles = 90 * np.ones((500, 3), dtype=np.float32)
    box_lengths[0,0] = 10.0

    f = DCDTrajectoryFile(temp, 'w')
    f.write(xyz, box_lengths, box_angles)
    f.close()

    f = DCDTrajectoryFile(temp)
    xyz2, box_lengths2, box_angles2 = f.read()
    f.close()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(box_lengths, box_lengths2)
    yield lambda: eq(box_angles, box_angles2)

@raises(ValueError)
def test_write_3():
    "checking"
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((600, 3), dtype=np.float32)

    with DCDTrajectoryFile(temp, 'w') as f:
        f.write(xyz, box_lengths)

def test_write_4():
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    box_lengths = 25 * np.ones((500, 3), dtype=np.float32)
    box_angles = 90 * np.ones((500, 3), dtype=np.float32)
    box_lengths[0,0] = 10.0

    f = DCDTrajectoryFile(temp, 'w')
    for i in range(len(xyz)):
        f.write(xyz[i], box_lengths[i], box_angles[i])
    f.close()

    f = DCDTrajectoryFile(temp)
    xyz2, box_lengths2, box_angles2 = f.read()
    f.close()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(box_lengths, box_lengths2)
    yield lambda: eq(box_angles, box_angles2)


def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with DCDTrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10,5,3))


@raises(IOError)
def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with DCDTrajectoryFile(temp, 'w', force_overwrite=False) as f:
        f.write(np.random.randn(10,5,3))


@raises(IOError)
def test_read_closed():
    f = DCDTrajectoryFile(fn_dcd)
    f.close()
    f.read()


@raises(IOError)
def test_write_closed():
    f = DCDTrajectoryFile(fn_dcd, 'w')
    f.close()
    f.write(np.random.randn(10,3,3))


def test_tell():
    with DCDTrajectoryFile(get_fn('frame0.dcd')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek():
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

