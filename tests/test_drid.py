from __future__ import print_function

import mdtraj as md
import numpy as np
import scipy.special
from mdtraj.geometry import compute_drid
from scipy.spatial.distance import euclidean, pdist, squareform

# To keep backwards compatibility
try:
    from scipy.stats import nanmean
except ImportError:
    from numpy import nanmean


def test_drid_1():
    n_frames = 1
    n_atoms = 20
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('X', chain)
    for i in range(n_atoms):
        top.add_atom('X', None, residue)

    t = md.Trajectory(xyz=np.random.RandomState(0).randn(n_frames, n_atoms, 3),
                      topology=top)
    # t contains no bonds
    got = compute_drid(t).reshape(n_frames, n_atoms, 3)

    for i in range(n_atoms):
        others = set(range(n_atoms)) - set([i])
        rd = 1 / np.array([euclidean(t.xyz[0, i], t.xyz[0, e]) for e in others])

        mean = np.mean(rd)
        second = np.mean((rd - mean) ** 2) ** (0.5)
        third = scipy.special.cbrt(np.mean((rd - mean) ** 3))

        ref = np.array([mean, second, third])
        np.testing.assert_array_almost_equal(got[0, i], ref, decimal=5)


def test_drid_2():
    n_frames = 3
    n_atoms = 11
    n_bonds = 5
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('X', chain)
    for i in range(n_atoms):
        top.add_atom('X', None, residue)

    random = np.random.RandomState(0)
    bonds = random.randint(n_atoms, size=(n_bonds, 2))
    for a, b in bonds:
        top.add_bond(top.atom(a), top.atom(b))

    t = md.Trajectory(xyz=random.randn(n_frames, n_atoms, 3), topology=top)
    got = compute_drid(t).reshape(n_frames, n_atoms, 3)

    for i in range(n_frames):
        recip = 1 / squareform(pdist(t.xyz[i]))
        recip[np.diag_indices(n=recip.shape[0])] = np.nan
        recip[bonds[:, 0], bonds[:, 1]] = np.nan
        recip[bonds[:, 1], bonds[:, 0]] = np.nan

        mean = nanmean(recip, axis=0)
        second = nanmean((recip - mean) ** 2, axis=0) ** (0.5)
        third = scipy.special.cbrt(nanmean((recip - mean) ** 3, axis=0))

        np.testing.assert_array_almost_equal(got[i, :, 0], mean, decimal=5)
        np.testing.assert_array_almost_equal(got[i, :, 1], second, decimal=5)
        np.testing.assert_array_almost_equal(got[i, :, 2], third, decimal=5)
