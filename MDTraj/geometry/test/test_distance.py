import time
import itertools
import numpy as np
from mdtraj.testing import eq, skipif

from mdtraj.geometry.distance import _HAVE_OPT
from mdtraj.geometry.distance import _distance, _displacement
from mdtraj.geometry.distance import _distance_mic, _displacement_mic
if _HAVE_OPT:
    from mdtraj.geometry.distance import _opt_distance, _opt_displacement


N_FRAMES = 10
N_ATOMS = 20

xyz = np.asarray(np.random.randn(N_FRAMES, N_ATOMS, 3), dtype=np.float32)
pairs = np.array(list(itertools.combinations(range(N_ATOMS), 2)), dtype=np.int32)
box = np.ascontiguousarray(np.random.randn(N_FRAMES, 3, 3) + 2*np.eye(3,3), dtype=np.float32)

@skipif(not _HAVE_OPT)
def test_0():
    a = _distance(xyz, pairs)
    b = _opt_distance(xyz, pairs)
    eq(a, b)

@skipif(not _HAVE_OPT)
def test_1():
    a = _displacement(xyz, pairs)
    b = _opt_displacement(xyz, pairs)
    eq(a, b)

def test_15():
    a = _displacement(xyz, pairs)
    c = _distance(xyz, pairs)
    eq(c, np.sqrt(np.sum(np.square(a), axis=2)))

@skipif(not _HAVE_OPT)
def test_2():
    a = _distance_mic(xyz, pairs, box)
    b = _opt_distance(xyz, pairs, box)
    eq(a, b, decimal=3)

@skipif(not _HAVE_OPT)
def test_3():
    a = _displacement_mic(xyz, pairs, box)
    b = _opt_displacement(xyz, pairs, box)
    eq(a, b, decimal=3)

def test_35():
    a = _displacement_mic(xyz, pairs, box)
    c = _distance_mic(xyz, pairs, box)
    eq(c, np.sqrt(np.sum(np.square(a), axis=2)))

def test_4():
    # using a really big box, we should get the same results with and without
    # pbcs
    box = np.array([[100, 0, 0], [0, 200, 0], [0, 0, 300]])
    box = np.zeros((N_FRAMES, 3, 3)) + box #broadcast it out
    a = _displacement_mic(xyz, pairs, box)
    b = _displacement(xyz, pairs)
    eq(a, b, decimal=3)

def test_5():
    # simple wrap around along the z axis.
    xyz = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 2.2]]])
    box = np.eye(3,3).reshape(1,3,3)
    result = _displacement_mic(xyz, np.array([[0,1]]), box)
    eq(result, np.array([[[0, 0, 0.2]]]))
    