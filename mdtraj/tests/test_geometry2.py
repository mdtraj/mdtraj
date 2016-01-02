import numpy as np
import mdtraj as md
from mdtraj.geometry import geometry2

np.random.seed(2)


def test_1():
    print()

    xyz = np.random.randn(1, 10, 3).astype(np.float32)
    lengths = np.random.rand(1, 3).astype(np.float32)
    angles = 90 * np.random.rand(1, 3).astype(np.float32)
    t = md.Trajectory(xyz, topology=None, unitcell_lengths=lengths, unitcell_angles=angles)

    pairs = np.array([[0,1]]).astype(np.int32)

    d1 = geometry2.compute_distances_triclinic(
        xyz,
        np.asarray(t.unitcell_vectors, order='c'),
        pairs)
    d2 = md.compute_distances(t, pairs, opt=False)

    print('d1 (my      )', d1)
    print('d2 (existing)', d2)
