##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
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

import itertools

import numpy as np
import pytest

import mdtraj as md
from mdtraj.geometry.distance import (
    _displacement,
    _displacement_mic,
    compute_displacements,
    compute_distances,
    compute_distances_core,
    compute_distances_t,
    find_closest_contact,
)
from mdtraj.testing import assert_allclose, eq


class Test_Distance:
    def test_generator(self, gen_random_ptraj):
        self.pairs2 = itertools.combinations(range(self.N_ATOMS), 2)
        a = compute_distances(self.ptraj, self.pairs)
        b = compute_distances(self.ptraj, self.pairs2)
        eq(a, b)

    def test_0(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=False, opt=True)
        b = compute_distances(self.ptraj, self.pairs, periodic=False, opt=False)
        eq(a, b)

    def test_compute_distances_core_nonperiodic(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=False, opt=True)
        b = compute_distances_core(
            self.ptraj.xyz,
            self.pairs,
            unitcell_vectors=self.ptraj.unitcell_vectors,
            periodic=False,
            opt=True,
        )
        eq(a, b)

        a = compute_distances(self.ptraj, self.pairs, periodic=False, opt=False)
        b = compute_distances_core(
            self.ptraj.xyz,
            self.pairs,
            unitcell_vectors=self.ptraj.unitcell_vectors,
            periodic=False,
            opt=False,
        )
        eq(a, b)

    def test_compute_distances_core_periodic(self, gen_random_ptraj):
        # opt
        a = compute_distances(self.ptraj, self.pairs, periodic=True, opt=True)
        b = compute_distances_core(
            self.ptraj.xyz,
            self.pairs,
            unitcell_vectors=self.ptraj.unitcell_vectors,
            periodic=True,
            opt=True,
        )
        eq(a, b, decimal=3)

        # no-opt
        a = compute_distances(self.ptraj, self.pairs, periodic=True, opt=False)
        b = compute_distances_core(
            self.ptraj.xyz,
            self.pairs,
            unitcell_vectors=self.ptraj.unitcell_vectors,
            periodic=True,
            opt=False,
        )
        eq(a, b, decimal=3)

    def test_1(self, gen_random_ptraj):
        a = compute_displacements(self.ptraj, self.pairs, periodic=False, opt=True)
        b = compute_displacements(self.ptraj, self.pairs, periodic=False, opt=False)
        eq(a, b)

    def test_2(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=False, opt=False)
        b = compute_displacements(self.ptraj, self.pairs, periodic=False, opt=False)
        eq(a, np.sqrt(np.sum(np.square(b), axis=2)))

    def test_3(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=False, opt=True)
        b = compute_displacements(self.ptraj, self.pairs, periodic=False, opt=True)
        eq(a, np.sqrt(np.sum(np.square(b), axis=2)))

    def test_0p(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=True, opt=True)
        b = compute_distances(self.ptraj, self.pairs, periodic=True, opt=False)
        eq(a, b, decimal=3)

    def test_1p(self, gen_random_ptraj):
        a = compute_displacements(self.ptraj, self.pairs, periodic=True, opt=True)
        b = compute_displacements(self.ptraj, self.pairs, periodic=True, opt=False)
        eq(a, b, decimal=3)

    def test_2p(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=True, opt=False)
        b = compute_displacements(self.ptraj, self.pairs, periodic=True, opt=False)
        assert a.shape == (len(self.ptraj), len(self.pairs))
        assert b.shape == (len(self.ptraj), len(self.pairs), 3), str(b.shape)
        b = np.sqrt(np.sum(np.square(b), axis=2))
        eq(a, b, decimal=5)

    def test_3p(self, gen_random_ptraj):
        a = compute_distances(self.ptraj, self.pairs, periodic=True, opt=True)
        b = compute_displacements(self.ptraj, self.pairs, periodic=True, opt=True)
        eq(a, np.sqrt(np.sum(np.square(b), axis=2)))

    def test_4(self, gen_random_ptraj):
        # using a really big box, we should get the same results with and without
        # pbcs
        box = np.array([[100, 0, 0], [0, 200, 0], [0, 0, 300]])
        box = np.zeros((self.N_FRAMES, 3, 3)) + box  # broadcast it out
        a = _displacement_mic(self.xyz, self.pairs, box, False)
        b = _displacement(self.xyz, self.pairs)
        eq(a, b, decimal=3)

    def test_5(self, gen_random_ptraj):
        # simple wrap around along the z axis.
        xyz = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 2.2]]])
        box = np.eye(3, 3).reshape(1, 3, 3)
        result = _displacement_mic(xyz, np.array([[0, 1]]), box, True)
        eq(result, np.array([[[0, 0, 0.2]]]))

    def test_6(self, get_fn):
        ext_ref = (
            np.array(
                [
                    17.4835,
                    22.2418,
                    24.2910,
                    22.5505,
                    12.8686,
                    22.1090,
                    7.4472,
                    22.4253,
                    19.8283,
                    20.6935,
                ],
            )
            / 10
        )
        traj = md.load(get_fn("test_good.nc"), top=get_fn("test.parm7"))
        self._run_amber_traj(traj, ext_ref)

    def test_7(self, get_fn):
        ext_ref = (
            np.array(
                [
                    30.9184,
                    23.9040,
                    25.3869,
                    28.0060,
                    25.9704,
                    24.6836,
                    23.0508,
                    27.1983,
                    24.4954,
                    26.7448,
                ],
            )
            / 10
        )
        traj = md.load(get_fn("test_bad.nc"), top=get_fn("test.parm7"))
        self._run_amber_traj(traj, ext_ref)

    def _run_amber_traj(self, traj, ext_ref):
        # Test triclinic case where simple approach in Tuckerman text does not
        # always work
        distopt = md.compute_distances(traj, [[0, 9999]], opt=True)
        distslw = md.compute_distances(traj, [[0, 9999]], opt=False)
        dispopt = md.compute_displacements(traj, [[0, 9999]], opt=True)
        dispslw = md.compute_displacements(traj, [[0, 9999]], opt=False)

        eq(distopt, distslw, decimal=5)
        eq(dispopt, dispslw, decimal=5)

        assert_allclose(distopt.flatten(), ext_ref, atol=2e-5)

        # Make sure distances from displacements are the same
        eq(np.sqrt((dispopt.squeeze() ** 2).sum(axis=1)), distopt.squeeze())
        eq(np.sqrt((dispslw.squeeze() ** 2).sum(axis=1)), distslw.squeeze())
        eq(dispopt, dispslw, decimal=5)

    def test_closest_contact(self, gen_random_ptraj):
        box_size = np.array([3.0, 4.0, 5.0])
        traj = md.Trajectory(xyz=self.xyz * box_size, topology=None)
        self._verify_closest_contact(gen_random_ptraj, traj)
        traj.unitcell_lengths = np.array([box_size for i in range(self.N_FRAMES)])
        traj.unitcell_angles = np.array([[90.0, 90.0, 90.0] for i in range(self.N_FRAMES)])
        self._verify_closest_contact(gen_random_ptraj, traj=traj)
        traj.unitcell_angles = np.array([[80.0, 90.0, 100.0] for i in range(self.N_FRAMES)])
        self._verify_closest_contact(gen_random_ptraj, traj=traj)

    def _verify_closest_contact(self, gen_random_ptraj, traj):
        group1 = np.array([i for i in range(self.N_ATOMS // 2)], dtype=int)
        group2 = np.array([i for i in range(self.N_ATOMS // 2, self.N_ATOMS)], dtype=int)
        contact = find_closest_contact(traj, group1, group2)
        pairs = np.array([(i, j) for i in group1 for j in group2], dtype=int)
        dists = md.compute_distances(traj, pairs, True)[0]
        _ = md.compute_distances(traj, pairs, False)[0]
        nearest = np.argmin(dists)
        eq(float(dists[nearest]), contact[2], decimal=5)
        assert (pairs[nearest, 0] == contact[0] and pairs[nearest, 1] == contact[1]) or (
            pairs[nearest, 0] == contact[1] and pairs[nearest, 1] == contact[0]
        )

    def test_distance_nan(self):
        xyz = np.array([[1, 1, 1], [2, 1, 1], [np.nan, np.nan, np.nan]]).reshape(1, 3, 3)
        dists = md.compute_distances(md.Trajectory(xyz=xyz, topology=None), [[0, 1]])
        assert np.isfinite(dists).all()

    def test_closest_contact_nan_pos(self, gen_random_ptraj):
        box_size = np.array([3.0, 4.0, 5.0])
        xyz = self.rng.standard_normal((2, 20, 3), dtype=np.float32)
        xyz *= box_size
        # Set the last frame to nan
        xyz[-1] = np.nan
        # Slice of the last frame, so nans should not cause troubles.
        xyz = xyz[:-1]
        traj = md.Trajectory(xyz=xyz, topology=None)
        self._verify_closest_contact(gen_random_ptraj, traj)

    def test_distance_t_inputs(self, gen_random_ptraj):
        incorrect_pairs = np.array((0, self.ptraj.n_atoms + 1))
        with pytest.raises(ValueError, match="atom_pairs"):
            compute_distances_t(self.ptraj, incorrect_pairs, self.times)

        incorrect_times = np.array((0, self.ptraj.n_frames + 1))
        with pytest.raises(ValueError, match="time_pairs"):
            compute_distances_t(self.ptraj, self.pairs, incorrect_times)

    def test_distances_t(self, gen_random_ptraj):
        a = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=True, opt=True)
        b = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=True, opt=False)
        # There is a precision difference between cython and python code for ``compute_distances_t``
        # on linux-aarch64 for certain seeds (e.g. seed=3220). Pinning for a lower precision check.
        eq(a, b, decimal=4)

        c = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=False, opt=True)
        d = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=False, opt=False)
        eq(c, d)

    def test_distances_t_at_0(self, gen_random_ptraj):
        self.times = np.array([[0, 0]], dtype=np.int32)
        a = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=True, opt=True)
        b = compute_distances_t(self.ptraj, self.pairs, self.times, periodic=True, opt=False)
        c = compute_distances(self.ptraj[:1], self.pairs, periodic=True, opt=True)
        d = compute_distances(self.ptraj[:1], self.pairs, periodic=True, opt=False)
        eq(a, c)
        eq(b, d)

    def _run_amber_traj_t(self, traj, ext_ref):
        # Test triclinic case where simple approach in Tuckerman text does not
        # always work
        _ = compute_distances_t(
            traj,
            atom_pairs=[[0, 9999]],
            time_pairs=[[0, 2]],
            opt=True,
        )
        _ = compute_distances_t(
            traj,
            atom_pairs=[[0, 9999]],
            time_pairs=[[0, 2]],
            opt=False,
        )

    def test_amber_t(self, get_fn):
        ext_ref = (
            np.array(
                [
                    17.4835,
                    22.2418,
                    24.2910,
                    22.5505,
                    12.8686,
                    22.1090,
                    7.4472,
                    22.4253,
                    19.8283,
                    20.6935,
                ],
            )
            / 10
        )
        traj = md.load(get_fn("test_good.nc"), top=get_fn("test.parm7"))
        self._run_amber_traj_t(traj, ext_ref)
