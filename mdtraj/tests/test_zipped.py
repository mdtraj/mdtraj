from __future__ import print_function, absolute_import, division

import os
import gzip
import bz2
import tempfile
import unittest
from mdtraj.testing import eq
from mdtraj.utils import open_maybe_zipped


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
        eq(open_maybe_zipped(fn, 'r').read(), u'COOKIE')

    def test_write_gz(self):
        fn = os.path.join(self.tmpdir, 'write.gz')
        with open_maybe_zipped(fn, 'w') as f:
            f.write(u'COOKIE')
        with gzip.GzipFile(fn, 'r') as f:
            eq(f.read().decode('utf-8'), u'COOKIE')

    def test_read_bz2(self):
        fn = os.path.join(self.tmpdir, 'read.bz2')
        with bz2.BZ2File(fn, 'w') as f:
            f.write('COOKIE'.encode('utf-8'))
        eq(open_maybe_zipped(fn, 'r').read(), u'COOKIE')

    def test_write_bz2(self):
        fn = os.path.join(self.tmpdir, 'write.bz2')
        with open_maybe_zipped(fn, 'w') as f:
            f.write(u'COOKIE')
        with bz2.BZ2File(fn, 'r') as f:
            eq(f.read().decode('utf-8'), u'COOKIE')
