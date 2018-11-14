##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Jason Swails
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
import tempfile

import mdtraj as md
import numpy as np
import pytest
from mdtraj.formats import AmberRestartFile, AmberNetCDFRestartFile
from mdtraj.testing import eq

fd1, temp1 = tempfile.mkstemp(suffix='.rst7')
fd2, temp2 = tempfile.mkstemp(suffix='.ncrst')
os.close(fd1)
os.close(fd2)


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    if os.path.exists(temp1): os.unlink(temp1)
    if os.path.exists(temp2): os.unlink(temp2)


def test_read_after_close(get_fn):
    f = AmberNetCDFRestartFile(get_fn('ncinpcrd.rst7'))
    assert eq(f.n_atoms, 2101)
    assert eq(f.n_frames, 1)

    f.close()

    with pytest.raises(IOError):
        f.read()


def test_shape(get_fn):
    with AmberRestartFile(get_fn('inpcrd')) as f:
        xyz, time, lengths, angles = f.read()

    assert eq(xyz.shape, (1, 2101, 3))
    assert eq(time.shape, (1,))
    assert eq(lengths, np.asarray([[30.2642725] * 3]))
    assert eq(angles, np.asarray([[109.471219] * 3]))


def test_shape_2(get_fn):
    with AmberNetCDFRestartFile(get_fn('ncinpcrd.rst7')) as f:
        xyz, time, lengths, angles = f.read()

    assert eq(xyz.shape, (1, 2101, 3))
    assert eq(time.shape, (1,))
    assert eq(lengths, np.asarray([[30.2642725] * 3]))
    assert eq(angles, np.asarray([[109.471219] * 3]))


def test_read_write_1():
    xyz = np.random.randn(1, 10, 3)
    time = np.random.randn(1)
    boxlengths = np.random.randn(1, 3)
    boxangles = np.random.randn(1, 3)

    with AmberRestartFile(temp1, 'w', force_overwrite=True) as f:
        f.write(xyz, time, boxlengths, boxangles)

    with AmberRestartFile(temp1) as f:
        a, b, c, d = f.read()
        assert eq(a, xyz)
        assert eq(b, time)
        assert eq(c, boxlengths)
        assert eq(d, boxangles)


def test_read_write_2():
    xyz = np.random.randn(1, 10, 3)
    time = np.random.randn(1)
    boxlengths = np.random.randn(1, 3)
    boxangles = np.random.randn(1, 3)

    with AmberNetCDFRestartFile(temp2, 'w', force_overwrite=True) as f:
        f.write(xyz, time, boxlengths, boxangles)

    with AmberNetCDFRestartFile(temp2) as f:
        a, b, c, d = f.read()
        assert eq(a, xyz)
        assert eq(b, time)
        assert eq(c, boxlengths)
        assert eq(d, boxangles)


def test_read_write_3(get_fn):
    traj = md.load(get_fn('frame0.nc'), top=get_fn('native.pdb'))
    traj[0].save(temp1)
    assert os.path.exists(temp1)
    rsttraj = md.load(temp1, top=get_fn('native.pdb'))
    eq(rsttraj.xyz, traj[0].xyz)
    os.unlink(temp1)
    traj.save(temp1)
    for i in range(traj.n_frames):
        assert os.path.exists('%s.%03d' % (temp1, i + 1))
        os.unlink('%s.%03d' % (temp1, i + 1))


def test_read_write_4(get_fn):
    traj = md.load(get_fn('frame0.nc'), top=get_fn('native.pdb'))
    traj[0].save(temp2)
    assert os.path.exists(temp2)
    rsttraj = md.load(temp2, top=get_fn('native.pdb'))
    eq(rsttraj.xyz, traj[0].xyz)
    os.unlink(temp2)
    traj.save(temp2)
    for i in range(traj.n_frames):
        assert os.path.exists('%s.%03d' % (temp2, i + 1))
        os.unlink('%s.%03d' % (temp2, i + 1))
