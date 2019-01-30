##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Gregory Bowman
# Contributors: Robert McGibbon
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
import unittest

import numpy as np
import tables

from mdtraj import io
from mdtraj.testing import eq

fd, temp = tempfile.mkstemp(suffix='.h5')
os.close(fd)


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.unlink(temp)


def test_overwrite_1():
    fid, fn = tempfile.mkstemp()
    try:
        a = np.arange(10)
        b = a + 1
        io.saveh(fn, a=a)
        io.saveh(fn, b=b)
        eq(io.loadh(fn, 'a'), a)
        eq(io.loadh(fn, 'b'), b)
    except:
        raise
    finally:
        if os.path.exists(fn):
            os.close(fid)
            os.unlink(fn)


def test_overwrite_2():
    fid, fn = tempfile.mkstemp()
    try:
        a = np.arange(10)
        b = a + 1
        io.saveh(fn, a=a)
        io.saveh(fn, a=b)
        eq(io.loadh(fn, 'a'), b)
    except:
        raise
    finally:
        if os.path.exists(fn):
            os.close(fid)
            os.unlink(fn)


class test_io(unittest.TestCase):
    def setUp(self):
        # setup() is called before very test and just creates a
        # temporary work space for reading/writing files.
        fid, self.filename1 = tempfile.mkstemp()
        fid, self.filename2 = tempfile.mkstemp()
        self.data = np.arange(10000, dtype=np.float32)

        # Write Data to an HDF5 file as a compressed CArray.
        hdf_file = tables.open_file(self.filename1, 'a')
        hdf_file.create_carray("/", "arr_0", tables.Float32Atom(),
                               self.data.shape, filters=io.COMPRESSION)
        hdf_file.root.arr_0[:] = self.data[:]
        hdf_file.flush()
        hdf_file.close()

    def test_load_1(self):
        # Load by specifying array name
        TestData = io.loadh(self.filename1, 'arr_0')
        eq(TestData, self.data)

    def test_load_2(self):
        # load using deferred=False
        TestData = io.loadh(self.filename1, deferred=False)['arr_0']
        eq(TestData, self.data)

    def test_load_3(self):
        # load using deferred=True
        deferred = io.loadh(self.filename1, deferred=True)
        eq(deferred['arr_0'], self.data)
        deferred.close()

    def test_save(self):
        # Save HDF5 to disk and load it back up
        io.saveh(self.filename2, self.data)
        TestData = io.loadh(self.filename2, 'arr_0')
        eq(TestData, self.data)

    def teardown(self):
        os.remove(self.filename1)
        os.remove(self.filename2)


class test_io_int(test_io):
    "Run the same test as the class above, but using int64 data"

    def setUp(self):
        # setup() is called before very test and just creates
        # a temporary work space for reading/writing files.
        fid, self.filename1 = tempfile.mkstemp()
        fid, self.filename2 = tempfile.mkstemp()
        self.data = np.arange(10000, dtype=np.int64)

        # Write Data to an HDF5 file as a compressed CArray.
        hdf_file = tables.open_file(self.filename1, 'a')
        hdf_file.create_carray("/", "arr_0", tables.Int64Atom(),
                               self.data.shape, filters=io.COMPRESSION)
        hdf_file.root.arr_0[:] = self.data[:]
        hdf_file.flush()
        hdf_file.close()


def test_groups():
    # Test to ensure that files are loaded correctly even if they contain
    # nested groups and stuff
    x = np.random.randn(10)
    y = np.random.randn(11)
    f = tables.open_file(temp, 'w')
    f.create_group(where='/', name='mygroup')
    f.create_array(where='/mygroup', name='myarray', obj=x)
    f.create_array(where='/', name='mya2', obj=y)
    f.close()

    assert eq(io.loadh(temp)['mygroup/myarray'], x)
    assert eq(io.loadh(temp)['mya2'], y)
    assert eq(io.loadh(temp, deferred=False)['mygroup/myarray'], x)
    assert eq(io.loadh(temp, deferred=False)['mya2'], y)
    assert eq(io.loadh(temp, 'mygroup/myarray'), x)
    assert eq(io.loadh(temp, 'mya2'), y)
