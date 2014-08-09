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


import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq
from mdtraj.geometry.alignment import rmsd_qcp, compute_translation_and_rotation


def test_trajectory_rmsd():
    t = md.load(get_fn('traj.h5'))
    for parallel in [True, False]:
        calculated = md.rmsd(t, t, 0, parallel=parallel)    
        reference = np.zeros(t.n_frames)
        for i in range(t.n_frames):
            reference[i] = rmsd_qcp(t.xyz[0], t.xyz[i])

        eq(calculated, reference, decimal=3)

def test_precentered_1():
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

def test_precentered_2():
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


def test_superpose_0():
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
        query = md.Trajectory(xyz=reference.xyz + 0.01*np.random.randn(*reference.xyz.shape), topology=None)
        query.superpose(reference, 0, atom_indices=indices)
        yield lambda: eq(reference.xyz, reference_xyz)

        new_centers = np.mean(query.xyz[0], axis=1)
        assert 80 < new_centers[0] < 120
        assert 80 < new_centers[1] < 120
        assert 80 < new_centers[2] < 120

def test_superpose_2():
    t1 = md.Trajectory(xyz=np.random.randn(1, 100, 3) + 100, topology=None)
    t2 = md.Trajectory(xyz=np.random.randn(1, 100, 3) + 100, topology=None)
    t2_copy = t2.xyz.copy()

    t1.superpose(t2)
    t1.superpose(t2, atom_indices=[1,2,3,4,5,6,7])

    # make sure that superposing doesn't alter the reference traj
    eq(t2.xyz, t2_copy)


def test_rmsd_atom_indices():
    native = md.load(get_fn('native.pdb'))
    t1 = md.load(get_fn('traj.h5'))

    atom_indices = np.arange(10)
    dist1 = md.rmsd(t1, native, atom_indices=atom_indices)

    t2 = md.load(get_fn('traj.h5'))
    t2.restrict_atoms(atom_indices)
    native.restrict_atoms(atom_indices)
    dist2 = md.rmsd(t2, native)
    
    eq(dist1, dist2)

if __name__ == '__main__':
    test_rmsd_atom_indices()

# def test_align_displace():
#     t = md.load(get_fn('traj.h5'))
#     t.center_coordinates()
#     pt = md.rmsd_cache(t, major='atom')
#     rmsd, rot = getMultipleAlignDisplaceRMSDs_atom_major(pt.cords, pt.cords, pt._traces(), pt._traces(), pt.cords, pt.cords, t.n_atoms, t.n_atoms, 0)
# 
#     for i in range(t.n_frames):
#         translation, rotation = compute_translation_and_rotation(t.xyz[0], t.xyz[i])
#         eq(rot[i], rotation)
#         eq(float(rmsd_qcp(t.xyz[0], t.xyz[i])), float(rmsd[i]), decimal=3)
