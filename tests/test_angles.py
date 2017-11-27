import itertools
import mdtraj as md
from mdtraj.testing import eq
import numpy as np


def test_angle_pbc_0():
    # Check that angles are calculated correctly across periodic boxes
    random = np.random.RandomState(0)
    boxLengths = [1, 2, 3]
    N = 1000

    X = np.cumsum(0.1 * random.randn(N * 3, 3), axis=0).reshape(N, 3, 3)

    center = X.mean(axis=1).reshape(N, 1, 3)
    cell = np.floor(center / boxLengths)
    X_cell = X - (cell * boxLengths)

    assert not np.all(X == X_cell)

    t0 = md.Trajectory(xyz=X, topology=None)
    t1 = md.Trajectory(xyz=X_cell, topology=None)
    t1.unitcell_vectors = np.tile(np.diag(boxLengths), [len(X), 1, 1])

    r0 = md.compute_angles(t0, [[0, 1, 2]], opt=False).reshape(-1)
    r1 = md.compute_angles(t1, [[0, 1, 2]], opt=False).reshape(-1)
    r2 = md.compute_angles(t0, [[0, 1, 2]], opt=True).reshape(-1)
    r3 = md.compute_angles(t1, [[0, 1, 2]], opt=True).reshape(-1)

    np.testing.assert_array_almost_equal(r0, r1, decimal=4)
    np.testing.assert_array_almost_equal(r2, r3, decimal=4)
    np.testing.assert_array_almost_equal(r0, r3, decimal=4)


def test_generator():
    N_FRAMES = 2
    N_ATOMS = 5
    xyz = np.asarray(np.random.randn(N_FRAMES, N_ATOMS, 3), dtype=np.float32)
    ptraj = md.Trajectory(xyz=xyz, topology=None)

    triplets = np.array(list(itertools.combinations(range(N_ATOMS), 3)), dtype=np.int32)
    triplets2 = itertools.combinations(range(N_ATOMS), 3)
    a = md.compute_angles(ptraj, triplets)
    b = md.compute_angles(ptraj, triplets2)
    eq(a, b)
