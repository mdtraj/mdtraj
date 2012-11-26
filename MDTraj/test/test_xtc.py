"""
Test the cython xtc module

Note, this file cannot be located in the xtc subdirectory, because that
directory is not a python package (it has no __init__.py) and is thus tests
there are not discovered by nose
"""

import os, tempfile
import numpy as np
from mdtraj import xtc, io
from mdtraj.testing import get_fn, eq

fn_xtc = get_fn('frame0.xtc')
pdb = get_fn('native.pdb')

temp = tempfile.mkstemp(suffix='.xtc')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)


def test_read_chunk1():
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=1)
    iofile = io.loadh(get_fn('frame0.xtc.h5'))
    eq(xyz, iofile['xyz'])
    eq(step, iofile['step'])
    eq(box, iofile['box'])
    eq(time, iofile['time'])
    eq(prec, iofile['prec'])

def test_read_chunk10():
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=10)
    iofile = io.loadh(get_fn('frame0.xtc.h5'))
    eq(xyz, iofile['xyz'])
    eq(step, iofile['step'])
    eq(box, iofile['box'])
    eq(time, iofile['time'])
    eq(prec, iofile['prec'])

def test_read_chunk1000():
    xyz, time, step, box, prec = xtc.read(fn_xtc, chunk=1000)
    iofile = io.loadh(get_fn('frame0.xtc.h5'))
    eq(xyz, iofile['xyz'])
    eq(step, iofile['step'])
    eq(box, iofile['box'])
    eq(time, iofile['time'])
    eq(prec, iofile['prec'])


def test_write_0():
    xyz = xtc.read(fn_xtc)[0]

    xtc.write(temp, xyz=xyz, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)
    eq(xyz, xyz2)


def test_write_1():
    xyz = np.around(np.random.randn(100, 10, 3), 3)

    xtc.write(temp, xyz=xyz, force_overwrite=True)
    xyz2, time2, step2, box2, prec2 = xtc.read(temp)
    eq(xyz, xyz2)



