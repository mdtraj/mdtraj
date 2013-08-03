import time
import itertools
import numpy as np
from mdtraj.testing import eq
from mdtraj.geometry._distance import distance
from mdtraj.geometry.distance import _compute_distances_xyz


def test_0():
    n_frames = 100
    n_atoms = 500

    xyz = np.asarray(np.random.randn(n_frames,n_atoms,3), dtype=np.float32)
    pairs = np.array(list(itertools.combinations(range(n_atoms), 2)), dtype=np.int32)

    t0 = time.time()
    a = distance(xyz, pairs)
    t1 = time.time()
    b = _compute_distances_xyz(xyz, pairs)
    t2 = time.time()

    print
    print 'sse  ', t1 - t0
    print 'numpy', t2 - t1

    eq(a, b)


def test_0():
    print
    n_frames = 1
    n_atoms = 5
    n_pairs = 1

    xyz = np.asarray(np.random.randn(n_frames,n_atoms,3), dtype=np.float32)
    pairs = np.array([[0,1]], dtype=np.int32)

    X = np.random.randn(3,3)
    X = np.asarray(X + X.T, dtype=np.float32)

    print X

    out = np.empty((n_frames, n_pairs))
    for i, (a,b) in enumerate(pairs):
        s1 = np.dot(np.linalg.inv(X), xyz[0,a,:])
        s2 = np.dot(np.linalg.inv(X), xyz[0,b,:])
        s12 = s2 - s1
        s12 = s12 - np.round(s12)
        r12 = np.dot(X, s12)
        out[i] = np.sqrt(np.sum(r12 * r12))

    print out
    print '\n\n'
    print distance(xyz, pairs, X)
