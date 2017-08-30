import mdtraj as md
import numpy as np

random = np.random.RandomState(0)


def _run_one_test(periodic, triclinic):
    n_atoms = 100
    cutoff = 1.5
    box_size = np.array([3.0, 4.0, 5.0])
    xyz = random.rand(1, n_atoms, 3) * box_size
    if periodic:
        unitcell_lengths = box_size
    else:
        unitcell_lengths = None
    if triclinic:
        unitcell_angles = np.array([80.0, 90.0, 100.0])
    else:
        unitcell_angles = None
    traj = md.Trajectory(xyz=xyz, topology=None, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)
    neighbors = md.compute_neighborlist(traj, cutoff)
    pairs = np.array([(i, j) for i in range(n_atoms) for j in range(n_atoms) if i < j])
    dists = md.compute_distances(traj, pairs, True, False)[0]
    for (i, j), d in zip(pairs, dists):
        if d < cutoff:
            assert (j in neighbors[i])
            assert (i in neighbors[j])
        else:
            assert (j not in neighbors[i])
            assert (i not in neighbors[j])


def test_nonperiodic():
    _run_one_test(False, False)


def test_periodic():
    _run_one_test(True, False)


def test_triclinic():
    _run_one_test(True, True)
