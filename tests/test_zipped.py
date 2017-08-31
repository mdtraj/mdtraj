from __future__ import print_function, absolute_import, division

import gzip
import bz2
from mdtraj.testing import eq
from mdtraj.utils import open_maybe_zipped


def test_read_gz(tmpdir):
    fn = '{}/read.gz'.format(tmpdir)
    with gzip.GzipFile(fn, 'w') as f:
        f.write('COOKIE'.encode('utf-8'))
    eq(open_maybe_zipped(fn, 'r').read(), u'COOKIE')


def test_write_gz(tmpdir):
    fn = '{}/write.gz'.format(tmpdir)
    with open_maybe_zipped(fn, 'w') as f:
        f.write(u'COOKIE')
    with gzip.GzipFile(fn, 'r') as f:
        eq(f.read().decode('utf-8'), u'COOKIE')


def test_read_bz2(tmpdir):
    fn = '{}/read.bz2'.format(tmpdir)
    with bz2.BZ2File(fn, 'w') as f:
        f.write('COOKIE'.encode('utf-8'))
    eq(open_maybe_zipped(fn, 'r').read(), u'COOKIE')


def test_write_bz2(tmpdir):
    fn = '{}/write.bz2'.format(tmpdir)
    with open_maybe_zipped(fn, 'w') as f:
        f.write(u'COOKIE')
    with bz2.BZ2File(fn, 'r') as f:
        eq(f.read().decode('utf-8'), u'COOKIE')
