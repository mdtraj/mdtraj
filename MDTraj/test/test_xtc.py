# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

"""
Test the cython xtc module

Note, this file cannot be located in the xtc subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import os, tempfile
import numpy as np
from mdtraj import xtc, io, XTCTrajectoryFile
from mdtraj.testing import get_fn, eq, DocStringFormatTester, raises

TestDocstrings = DocStringFormatTester(xtc, error_on_none=True)

fn_xtc = get_fn('frame0.xtc')
pdb = get_fn('native.pdb')

temp = tempfile.mkstemp(suffix='.xtc')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_read_chunk1():
    with XTCTrajectoryFile(fn_xtc, 'r', chunk_size_multiplier=0.5) as f:
        xyz, time, step, box = f.read()

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    yield lambda: eq(xyz, iofile['xyz'])
    yield lambda: eq(step, iofile['step'])
    yield lambda: eq(box, iofile['box'])
    yield lambda: eq(time, iofile['time'])


def test_read_stride():
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
         xyz, time, step, box = f.read(stride=3)
    yield lambda: eq(xyz, iofile['xyz'][::3])
    yield lambda: eq(step, iofile['step'][::3])
    yield lambda: eq(box, iofile['box'][::3])
    yield lambda: eq(time, iofile['time'][::3])


def test_read_atomindices_1():
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
         xyz, time, step, box = f.read(atom_indices=[0,1,2])
    yield lambda: eq(xyz, iofile['xyz'][:, [0,1,2]])
    yield lambda: eq(step, iofile['step'])
    yield lambda: eq(box, iofile['box'])
    yield lambda: eq(time, iofile['time'])

def test_read_atomindices_2():
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
         xyz, time, step, box = f.read(atom_indices=slice(None, None, 2))
    yield lambda: eq(xyz, iofile['xyz'][:, ::2])
    yield lambda: eq(step, iofile['step'])
    yield lambda: eq(box, iofile['box'])
    yield lambda: eq(time, iofile['time'])

def test_read_chunk2():
    with XTCTrajectoryFile(fn_xtc, 'r', chunk_size_multiplier=1) as f:
        xyz, time, step, box = f.read()

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    yield lambda: eq(xyz, iofile['xyz'])
    yield lambda: eq(step, iofile['step'])
    yield lambda: eq(box, iofile['box'])
    yield lambda: eq(time, iofile['time'])

def test_read_chunk3():
    with XTCTrajectoryFile(fn_xtc, chunk_size_multiplier=2) as f:
        xyz, time, step, box = f.read(n_frames=100)

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    yield lambda: eq(xyz, iofile['xyz'][:100])
    yield lambda: eq(step, iofile['step'][:100])
    yield lambda: eq(box, iofile['box'][:100])
    yield lambda: eq(time, iofile['time'][:100])

def test_write_0():
    with XTCTrajectoryFile(fn_xtc) as f:
        xyz = f.read()[0]

    f = XTCTrajectoryFile(temp, 'w')
    f.write(xyz)
    f.close()

    with XTCTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2 = f.read()
    eq(xyz, xyz2)


def test_write_1():
    xyz = np.around(np.random.randn(100, 10, 3), 3)
    time = np.random.randn(100)
    step = np.arange(100)
    box = np.random.randn(100,3,3)

    with XTCTrajectoryFile(temp, 'w') as f:
        f.write(xyz, time=time, step=step, box=box)
    with XTCTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2 = f.read()

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(box, box2)


def test_write_2():
    xyz = np.around(np.random.randn(100, 10, 3), 3)
    time = np.random.randn(100)
    step = np.arange(100)
    box = np.random.randn(100,3,3)

    with XTCTrajectoryFile(temp, 'w') as f:
        for i in range(len(xyz)):
            f.write(xyz[i], time=time[i], step=step[i], box=box[i])
    with XTCTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2 = f.read()

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(box, box2)

@raises(ValueError)
def test_write_error_0():
    xyz = np.around(np.random.randn(100, 10, 3), 3)

    with XTCTrajectoryFile(temp, 'r') as f:
        f.write(xyz)

@raises(ValueError)
def test_read_error_0():
    xyz = np.around(np.random.randn(100, 10, 3), 3)

    with XTCTrajectoryFile(temp, 'w') as f:
        f.read(xyz)

@raises(IOError)
def test_read_error_1():
    XTCTrajectoryFile('/tmp/sdfsdfsdf')

@raises(IOError)
def test_read_error_2():
    XTCTrajectoryFile(get_fn('frame0.dcd')).read()

@raises(ValueError)
def test_xtc_write_wierd_0():
    x0 = np.random.randn(100,3,3)
    x1 = np.random.randn(100,9,3)
    with XTCTrajectoryFile(temp, 'w') as f:
        f.write(x0)
        f.write(x1)

    xr = XTCTrajectoryFile(temp).read()[0]


    print xr.shape
