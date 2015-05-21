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

import bz2
import glob
import gzip
import os
import os.path
import tables
import tempfile
import unittest

from mdtraj.testing import eq
from mdtraj import io

import numpy as np

fd, temp = tempfile.mkstemp(suffix='.h5')
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.close(fd)
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


class test_open_maybe_zipped(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def teardown(self):
        for fn in glob.glob(os.path.join(self.tmpdir, '*')):
            os.remove(fn)
        os.rmdir(self.tmpdir)

    def test_read_gz(self):
        fn = os.path.join(self.tmpdir, 'read.gz')
        with gzip.GzipFile(fn, 'w') as f:
            f.write('COOKIE'.encode('utf-8'))
        eq(io.open_maybe_zipped(fn, 'r').read(), u'COOKIE')

    def test_write_gz(self):
        fn = os.path.join(self.tmpdir, 'write.gz')
        with io.open_maybe_zipped(fn, 'w') as f:
            f.write(u'COOKIE')
        with gzip.GzipFile(fn, 'r') as f:
            eq(f.read().decode('utf-8'), u'COOKIE')

    def test_read_bz2(self):
        fn = os.path.join(self.tmpdir, 'read.bz2')
        with bz2.BZ2File(fn, 'w') as f:
            f.write('COOKIE'.encode('utf-8'))
        eq(io.open_maybe_zipped(fn, 'r').read(), u'COOKIE')

    def test_write_bz2(self):
        fn = os.path.join(self.tmpdir, 'write.bz2')
        with io.open_maybe_zipped(fn, 'w') as f:
            f.write(u'COOKIE')
        with bz2.BZ2File(fn, 'r') as f:
            eq(f.read().decode('utf-8'), u'COOKIE')


class test_io(unittest.TestCase):
    def setUp(self):
        """setup() is called before very test and just creates a temporary work space for reading/writing files."""
        fid, self.filename1 = tempfile.mkstemp()
        fid, self.filename2 = tempfile.mkstemp()
        self.data = np.arange(10000, dtype=np.float32)

        #Write Data to an HDF5 file as a compressed CArray.
        hdfFile = tables.File(self.filename1, 'a')
        #The filter is the same used to save MSMB2 data
        hdfFile.createCArray("/", "arr_0", tables.Float32Atom(), self.data.shape, filters=io.COMPRESSION)
        hdfFile.root.arr_0[:] = self.data[:]
        hdfFile.flush()
        hdfFile.close()

    def test_load_1(self):
        "Load by specifying array name"
        TestData = io.loadh(self.filename1, 'arr_0')
        eq(TestData, self.data)

    def test_load_2(self):
        "load using deferred=False"
        TestData = io.loadh(self.filename1, deferred=False)['arr_0']
        eq(TestData, self.data)

    def test_load_2(self):
        "load using deferred=True"
        deferred = io.loadh(self.filename1, deferred=True)
        eq(deferred['arr_0'], self.data)
        deferred.close()

    def test_save(self):
        """Save HDF5 to disk and load it back up"""
        io.saveh(self.filename2, self.data)
        TestData = io.loadh(self.filename2, 'arr_0')
        eq(TestData, self.data)

    def teardown(self):
        os.remove(self.filename1)
        os.remove(self.filename2)


class test_io_int(test_io):
    "Run the same test as the class above, but using int64 data"
    def setUp(self):
        """setup() is called before very test and just creates a temporary work space for reading/writing files."""
        fid, self.filename1 = tempfile.mkstemp()
        fid, self.filename2 = tempfile.mkstemp()
        self.data = np.arange(10000, dtype=np.int64)

        #Write Data to an HDF5 file as a compressed CArray.
        hdfFile = tables.File(self.filename1, 'a')
        #The filter is the same used to save MSMB2 data
        hdfFile.createCArray("/", "arr_0", tables.Int64Atom(), self.data.shape, filters=io.COMPRESSION)
        hdfFile.root.arr_0[:] = self.data[:]
        hdfFile.flush()
        hdfFile.close()


def test_groups():
    """Test to ensure that files are loaded correctly even if they contain nested
    groups and stuff"""
    x = np.random.randn(10)
    y = np.random.randn(11)
    f = tables.openFile(temp, 'w')
    f.createGroup(where='/', name='mygroup')
    if tables.__version__ >= '3.0.0':
        f.createArray(where='/mygroup', name='myarray', obj=x)
        f.createArray(where='/', name='mya2', obj=y)
    else:
        f.createArray(where='/mygroup', name='myarray', object=x)
        f.createArray(where='/', name='mya2', object=y)
    f.close()

    yield lambda: eq(io.loadh(temp)['mygroup/myarray'], x)
    yield lambda: eq(io.loadh(temp)['mya2'], y)
    yield lambda: eq(io.loadh(temp, deferred=False)['mygroup/myarray'], x)
    yield lambda: eq(io.loadh(temp, deferred=False)['mya2'], y)
    yield lambda: eq(io.loadh(temp, 'mygroup/myarray'), x)
    yield lambda: eq(io.loadh(temp, 'mya2'), y)
