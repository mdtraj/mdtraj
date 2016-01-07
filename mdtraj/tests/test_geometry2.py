import numpy as np
import mdtraj as md
from mdtraj.geometry import geometry


def test_1():
    random = np.random.RandomState(5)
    xyz = random.randn(1, 10, 3).astype(np.float32)

    lengths = [[ 0.63979518, 0.98562443,  0.25909761]]
    angles = [[ 72.22472382,  78.34347534,  83.04747009]]

    t = md.Trajectory(xyz, topology=None, unitcell_lengths=lengths, unitcell_angles=angles)

    pairs = np.array([[0,1]])

    d1 = geometry.compute_distances(
        xyz,
        np.asarray(t.unitcell_vectors, order='c'),
        pairs)
    d2 = md.compute_distances(t, pairs, opt=False)



    print('d1 (my      )', d1)
    print('d2 (existing)', d2)
