##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################


import mdtraj as md
import numpy as np
from mdtraj.geometry.alignment import rmsd_qcp, compute_average_structure
from mdtraj.testing import eq

np.random.seed(52)


def test_trajectory_rmsd(get_fn):
    t = md.load(get_fn('traj.h5'))
    for parallel in [True, False]:
        calculated = md.rmsd(t, t, 0, parallel=parallel)
        reference = np.zeros(t.n_frames)
        for i in range(t.n_frames):
            reference[i] = rmsd_qcp(t.xyz[0], t.xyz[i])

        eq(calculated, reference, decimal=3)


def test_precentered_1(get_fn):
    # test rmsd against the numpy version, using the same trajectory
    # as target and reference
    t1 = md.load(get_fn('traj.h5'), stride=10)
    t2 = md.load(get_fn('traj.h5'), stride=10)
    # don't center t1, and use it without precentered
    # explicitly center t2, and use *with* precentered

    for parallel in [True, False]:
        t2.center_coordinates()
        eq(t1.n_frames, t2.n_frames)
        for i in range(t1.n_frames):
            ref = np.zeros(t1.n_frames)
            for j in range(t1.n_frames):
                ref[j] = rmsd_qcp(t1.xyz[j], t1.xyz[i])
            val1 = md.rmsd(t1, t1, i, parallel=parallel, precentered=False)
            val2 = md.rmsd(t2, t2, i, parallel=parallel, precentered=True)

            eq(ref, val1, decimal=3)
            eq(val1, val2)


def test_precentered_2(get_fn):
    # test rmsd against the numpy version, using the difference
    # trajectories as target and reference
    t1_a = md.load(get_fn('traj.h5'), stride=10)
    t2_a = md.load(get_fn('traj.h5'), stride=10)
    t1_b = md.load(get_fn('traj.h5'), stride=10)
    t2_b = md.load(get_fn('traj.h5'), stride=10)
    # don't center t1, and use it without precentered
    # explicitly center t2, and use *with* precentered

    t2_a.center_coordinates()
    t2_b.center_coordinates()

    for parallel in [True, False]:
        for i in range(t1_b.n_frames):
            ref = np.zeros(t1_a.n_frames)
            for j in range(t1_a.n_frames):
                ref[j] = rmsd_qcp(t1_a.xyz[j], t1_b.xyz[i])
            val1 = md.rmsd(t1_a, t1_b, i, parallel=parallel, precentered=False)
            val2 = md.rmsd(t2_a, t2_b, i, parallel=parallel, precentered=True)

            eq(ref, val1, decimal=3)
            eq(val1, val2, decimal=4)


def test_superpose_0(get_fn):
    t1 = md.load(get_fn('traj.h5'))
    reference_rmsd = md.rmsd(t1, t1, 0)

    t1.superpose(t1, 0)
    displ_rmsd = np.zeros(t1.n_frames)
    for i in range(t1.n_frames):
        delta = t1.xyz[i] - t1.xyz[0]
        displ_rmsd[i] = (delta ** 2.0).sum(1).mean() ** 0.5

    eq(reference_rmsd, displ_rmsd, decimal=5)


def test_superpose_1():
    # make one frame far from the origin
    reference = md.Trajectory(xyz=np.random.randn(1, 100, 3) + 100, topology=None)
    reference_xyz = reference.xyz.copy()

    for indices in [None, np.arange(90)]:
        # make another trajectory in a similar rotational state
        query = md.Trajectory(xyz=reference.xyz + 0.01 * np.random.randn(*reference.xyz.shape), topology=None)
        query.superpose(reference, 0, atom_indices=indices)
        assert eq(reference.xyz, reference_xyz)

        new_centers = np.mean(query.xyz[0], axis=1)
        assert 80 < new_centers[0] < 120
        assert 80 < new_centers[1] < 120
        assert 80 < new_centers[2] < 120


def test_superpose_2():
    t1 = md.Trajectory(xyz=np.random.randn(1, 100, 3) + 100, topology=None)
    t2 = md.Trajectory(xyz=np.random.randn(1, 100, 3) + 100, topology=None)
    t2_copy = t2.xyz.copy()

    t1.superpose(t2)
    t1.superpose(t2, atom_indices=[1, 2, 3, 4, 5, 6, 7])

    # make sure that superposing doesn't alter the reference traj
    eq(t2.xyz, t2_copy)


def test_superpose_refinds():
    # make one frame far from the origin
    normal = np.random.randn(1, 100, 3)
    normal_xyz = normal.copy()

    flipped = np.zeros_like(normal)
    flipped[:, :50, :] = normal[:, 50:, :]
    flipped[:, 50:, :] = normal[:, :50, :]
    flipped_xyz = flipped.copy()

    normal = md.Trajectory(xyz=normal, topology=None)
    flipped = md.Trajectory(xyz=flipped, topology=None)

    normal.superpose(flipped, atom_indices=np.arange(0, 50), ref_atom_indices=np.arange(50, 100))
    eq(normal.xyz, normal_xyz)

    flipped.superpose(normal, atom_indices=np.arange(50, 100), ref_atom_indices=np.arange(0, 50))
    eq(flipped.xyz, flipped_xyz)

    normal.superpose(flipped)
    assert not np.allclose(normal.xyz, normal_xyz)


def test_rmsd_atom_indices(get_fn):
    native = md.load(get_fn('native.pdb'))
    t1 = md.load(get_fn('traj.h5'))

    atom_indices = np.arange(10)
    dist1 = md.rmsd(t1, native, atom_indices=atom_indices)

    t2 = md.load(get_fn('traj.h5'))
    t2.restrict_atoms(atom_indices)
    native.restrict_atoms(atom_indices)
    dist2 = md.rmsd(t2, native)

    eq(dist1, dist2)


def test_rmsd_ref_ainds(get_fn):
    native = md.load(get_fn('native.pdb'))
    t1 = md.load(get_fn('traj.h5'))

    atom_indices = np.arange(10)
    dist1 = md.rmsd(t1, native, atom_indices=atom_indices,
                    ref_atom_indices=atom_indices)

    bad_atom_indices = np.arange(10, 20)
    t2 = md.load(get_fn('traj.h5'))
    dist2 = md.rmsd(t2, native, atom_indices=atom_indices,
                    ref_atom_indices=bad_atom_indices)

    assert np.all(dist2 > dist1)


def test_average_structure(get_fn):
    traj = md.load(get_fn('frame0.dcd'), top=get_fn('frame0.pdb'))
    average = compute_average_structure(traj.xyz)

    # The mean RMSD to the average structure should be less than to any individual frame.
    sum1 = 0
    sum2 = 0
    for i in range(traj.n_frames):
        sum1 += rmsd_qcp(traj.xyz[0], traj.xyz[i])
        sum2 += rmsd_qcp(average, traj.xyz[i])
    assert sum2 < sum1

def test_trajectory_rmsf(get_fn):
    t = md.load(get_fn('traj.h5'))
    for parallel in [True, False]:
        calculated = md.rmsf(t, t, 0, parallel=parallel)
        t.superpose(t, 0)
        avg_xyz = np.average(t.xyz, axis=0)
        reference = np.sqrt(3*np.mean((t.xyz - avg_xyz)**2, axis=(0, 2)))
        assert np.sum(np.abs(calculated)) > 0 # check trivial error
        eq(calculated, reference, decimal=3)

def test_trajectory_rmsf_aligned(get_fn):
    t = md.load(get_fn('traj.h5'))
    for parallel in [True, False]:
        # testing different set of atoms for alignment and RMSF calculation
        atom_indices = range(int(t.n_atoms/2))
        rmsf_indices = range(int(t.n_atoms/2), t.n_atoms)
        t.superpose(t, 99, atom_indices=atom_indices, parallel=False)
        calculated = md.rmsf(t, None, atom_indices=rmsf_indices, parallel=parallel)
        avg_xyz = np.average(t.xyz, axis=0)
        reference = np.sqrt(3*np.mean((t.xyz - avg_xyz)**2, axis=(0, 2)))[rmsf_indices]
        assert np.sum(np.abs(calculated)) > 0 # check trivial error
        eq(calculated, reference, decimal=3)

def test_rmsd_atom_indices_vs_ref_indices():
    n_frames = 1
    n_atoms_1 = 1
    n_atoms_2 = 2

    top_1 = md.Topology()
    top_1.add_chain()
    top_1.add_residue('RS2', top_1.chain(0))
    top_1.add_atom('A2', 'H', top_1.residue(0))

    top_2 = md.Topology()
    top_2.add_chain()
    top_2.add_residue('RS1', top_2.chain(0))
    top_2.add_atom('A1', 'H', top_2.residue(0))
    top_2.add_chain()
    top_2.add_residue('RS2', top_2.chain(1))
    top_2.add_atom('A2', 'H', top_2.residue(1))
    # here the 2nd chain in the top_2 is rmsd-compatible to the one in the top_1 so we should be able to compute rsmd between them.

    trj_1 = md.Trajectory(np.random.RandomState(0).randn(n_frames, n_atoms_1, 3), top_1)
    trj_2 = md.Trajectory(np.random.RandomState(0).randn(n_frames, n_atoms_2, 3), top_2)

    md.rmsd(trj_1, trj_2, atom_indices=[0], ref_atom_indices=[1])
    md.rmsd(trj_2, trj_1, atom_indices=[1], ref_atom_indices=[0])
    # is this don't fail then it's good no matter the result
