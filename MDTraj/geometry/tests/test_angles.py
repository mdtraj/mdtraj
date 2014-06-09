import mdtraj as md
import numpy as np
import matplotlib.pyplot as pp
import scipy.spatial

def test_angle_pbc_0():
    # Check that angles are calculated correctly accross periodic boxes
    random = np.random.RandomState(0)
    boxLengths = [1,2,3]
    N = 1000

    X = np.cumsum(0.1 * random.randn(N*3, 3), axis=0).reshape(N, 3, 3)

    center = X.mean(axis=1, keepdims=True)
    cell = np.floor(center / boxLengths)
    X_cell = X - (cell * boxLengths)

    assert not np.all(X == X_cell)

    t0 = md.Trajectory(xyz=X, topology=None)
    t1 = md.Trajectory(xyz=X_cell, topology=None)
    t1.unitcell_vectors = np.tile(np.diag(boxLengths), [len(X), 1, 1])

    r0 = md.compute_angles(t0, [[0,1,2]], opt=False).reshape(-1)
    r1 = md.compute_angles(t1, [[0,1,2]], opt=False).reshape(-1)
    r2 = md.compute_angles(t0, [[0,1,2]], opt=True).reshape(-1)
    r3 = md.compute_angles(t1, [[0,1,2]], opt=True).reshape(-1)
    
    print np.max(np.abs(r0-r1))
    print np.max(np.abs(r2-r3))    

    i = 458
    print r0[i], r1[i], r2[i], r3[i]
    print X[i]
    print X_cell[i]

    np.testing.assert_array_almost_equal(r0, r1, decimal=4)
    np.testing.assert_array_almost_equal(r2, r3, decimal=4)
    np.testing.assert_array_almost_equal(r0, r3, decimal=4)
