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
