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

import time
import itertools
import pytest
import numpy as np
import mdtraj as md

from mdtraj.testing import eq, assert_allclose
from mdtraj.geometry.distance import (
    compute_distances_core,
    compute_distances,
    compute_distances_t,
    compute_displacements,
    find_closest_contact,
)
from mdtraj.geometry.distance import _displacement_mic, _displacement

N_FRAMES = 20
N_ATOMS = 20

xyz = np.asarray(np.random.randn(N_FRAMES, N_ATOMS, 3), dtype=np.float32)
pairs = np.array(list(itertools.combinations(range(N_ATOMS), 2)), dtype=np.int32)
times = np.array([[i, 0] for i in range(N_FRAMES)[::2]], dtype=np.int32)

ptraj = md.Trajectory(xyz=xyz, topology=None)
ptraj.unitcell_vectors = np.ascontiguousarray(np.random.randn(N_FRAMES, 3, 3) + 2 * np.eye(3, 3), dtype=np.float32)


def test_generator():
    pairs2 = itertools.combinations(range(N_ATOMS), 2)
    a = compute_distances(ptraj, pairs)
    b = compute_distances(ptraj, pairs2)
    eq(a, b)


def test_0():
    a = compute_distances(ptraj, pairs, periodic=False, opt=True)
    b = compute_distances(ptraj, pairs, periodic=False, opt=False)
    eq(a, b)

def test_compute_distances_core_nonperiodic():

    a = compute_distances(ptraj, pairs, periodic=False, opt=True)
    b = compute_distances_core(ptraj.xyz,
                               pairs,
                               unitcell_vectors=ptraj.unitcell_vectors,
                               periodic=False,
                               opt=True,
                               )
    eq(a, b)

    a = compute_distances(ptraj, pairs, periodic=False, opt=False)
    b = compute_distances_core(ptraj.xyz,
                               pairs,
                               unitcell_vectors=ptraj.unitcell_vectors,
                               periodic=False,
                               opt=False,
                               )
    eq(a, b)

def test_compute_distances_core_periodic():

    # opt
    a = compute_distances(ptraj, pairs, periodic=True, opt=True)
    b = compute_distances_core(ptraj.xyz,
                               pairs,
                               unitcell_vectors=ptraj.unitcell_vectors,
                               periodic=True,
                               opt=True,
                               )
    eq(a, b, decimal=3)

    # no-opt
    a = compute_distances(ptraj, pairs, periodic=True, opt=False)
    b = compute_distances_core(ptraj.xyz,
                               pairs,
                               unitcell_vectors=ptraj.unitcell_vectors,
                               periodic=True,
                               opt=False,
                               )
    eq(a, b, decimal=3)

def test_1():
    a = compute_displacements(ptraj, pairs, periodic=False, opt=True)
    b = compute_displacements(ptraj, pairs, periodic=False, opt=False)
    eq(a, b)


def test_2():
    a = compute_distances(ptraj, pairs, periodic=False, opt=False)
    b = compute_displacements(ptraj, pairs, periodic=False, opt=False)
    eq(a, np.sqrt(np.sum(np.square(b), axis=2)))


def test_3():
    a = compute_distances(ptraj, pairs, periodic=False, opt=True)
    b = compute_displacements(ptraj, pairs, periodic=False, opt=True)
    eq(a, np.sqrt(np.sum(np.square(b), axis=2)))


def test_0p():
    a = compute_distances(ptraj, pairs, periodic=True, opt=True)
    b = compute_distances(ptraj, pairs, periodic=True, opt=False)
    print(a, b)
    eq(a, b, decimal=3)


def test_1p():
    a = compute_displacements(ptraj, pairs, periodic=True, opt=True)
    b = compute_displacements(ptraj, pairs, periodic=True, opt=False)
    eq(a, b, decimal=3)


def test_2p():
    a = compute_distances(ptraj, pairs, periodic=True, opt=False)
    b = compute_displacements(ptraj, pairs, periodic=True, opt=False)
    assert a.shape == (len(ptraj), len(pairs))
    assert b.shape == (len(ptraj), len(pairs), 3), str(b.shape)
    b = np.sqrt(np.sum(np.square(b), axis=2))
    eq(a, b, decimal=5)


def test_3p():
    a = compute_distances(ptraj, pairs, periodic=True, opt=True)
    b = compute_displacements(ptraj, pairs, periodic=True, opt=True)
    print(a,  b)
    eq(a, np.sqrt(np.sum(np.square(b), axis=2)))


def test_4():
    # using a really big box, we should get the same results with and without
    # pbcs
    box = np.array([[100, 0, 0], [0, 200, 0], [0, 0, 300]])
    box = np.zeros((N_FRAMES, 3, 3)) + box  # broadcast it out
    a = _displacement_mic(xyz, pairs, box, False)
    b = _displacement(xyz, pairs)
    eq(a, b, decimal=3)


def test_5():
    # simple wrap around along the z axis.
    xyz = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 2.2]]])
    box = np.eye(3, 3).reshape(1, 3, 3)
    result = _displacement_mic(xyz, np.array([[0, 1]]), box, True)
    eq(result, np.array([[[0, 0, 0.2]]]))


def test_6(get_fn):
    ext_ref = np.array([17.4835, 22.2418, 24.2910, 22.5505, 12.8686, 22.1090,
                        7.4472, 22.4253, 19.8283, 20.6935]) / 10
    traj = md.load(get_fn('test_good.nc'), top=get_fn('test.parm7'))
    _run_amber_traj(traj, ext_ref)


def test_7(get_fn):
    ext_ref = np.array([30.9184, 23.9040, 25.3869, 28.0060, 25.9704, 24.6836,
                        23.0508, 27.1983, 24.4954, 26.7448]) / 10
    traj = md.load(get_fn('test_bad.nc'), top=get_fn('test.parm7'))
    _run_amber_traj(traj, ext_ref)


def _run_amber_traj(traj, ext_ref):
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


def test_closest_contact():
    box_size = np.array([3.0, 4.0, 5.0])
    traj = md.Trajectory(xyz=xyz * box_size, topology=None)
    _verify_closest_contact(traj)
    traj.unitcell_lengths = np.array([box_size for i in range(N_FRAMES)])
    traj.unitcell_angles = np.array([[90.0, 90.0, 90.0] for i in range(N_FRAMES)])
    _verify_closest_contact(traj)
    traj.unitcell_angles = np.array([[80.0, 90.0, 100.0] for i in range(N_FRAMES)])
    _verify_closest_contact(traj)


def _verify_closest_contact(traj):
    group1 = np.array([i for i in range(N_ATOMS // 2)], dtype=int)
    group2 = np.array([i for i in range(N_ATOMS // 2, N_ATOMS)], dtype=int)
    contact = find_closest_contact(traj, group1, group2)
    pairs = np.array([(i, j) for i in group1 for j in group2], dtype=int)
    dists = md.compute_distances(traj, pairs, True)[0]
    dists2 = md.compute_distances(traj, pairs, False)[0]
    nearest = np.argmin(dists)
    eq(float(dists[nearest]), contact[2], decimal=5)
    assert ((pairs[nearest, 0] == contact[0] and pairs[nearest, 1] == contact[1]) or (
    pairs[nearest, 0] == contact[1] and pairs[nearest, 1] == contact[0]))


def test_distance_nan():
    xyz = np.array([[1, 1, 1], [2, 1, 1], [np.nan, np.nan, np.nan]]).reshape(1, 3, 3)
    dists = md.compute_distances(md.Trajectory(xyz=xyz, topology=None), [[0, 1]])
    assert np.isfinite(dists).all()


def test_closest_contact_nan_pos():
    box_size = np.array([3.0, 4.0, 5.0])
    xyz = np.asarray(np.random.randn(2, 20, 3), dtype=np.float32)
    xyz *= box_size
    # Set the last frame to nan
    xyz[-1] = np.nan
    # Slice of the last frame, so nans should not cause troubles.
    xyz = xyz[:-1]
    traj = md.Trajectory(xyz=xyz, topology=None)
    _verify_closest_contact(traj)


def test_distance_t_inputs():
    incorrect_pairs = np.array((0, ptraj.n_atoms+1))
    with pytest.raises(ValueError, match='atom_pairs'):
        compute_distances_t(ptraj, incorrect_pairs, times)

    incorrect_times = np.array((0, ptraj.n_frames+1))
    with pytest.raises(ValueError, match='time_pairs'):
        compute_distances_t(ptraj, pairs, incorrect_times)


def test_distances_t(get_fn):
    a = compute_distances_t(ptraj, pairs, times, periodic=True, opt=True)
    b = compute_distances_t(ptraj, pairs, times, periodic=True, opt=False)
    eq(a, b)
    c = compute_distances_t(ptraj, pairs, times, periodic=False, opt=True)
    d = compute_distances_t(ptraj, pairs, times, periodic=False, opt=False)
    eq(c, d)

def test_distances_t_at_0(get_fn):
    times = np.array([[0, 0]], dtype=np.int32)
    a = compute_distances_t(ptraj, pairs, times, periodic=True, opt=True)
    b = compute_distances_t(ptraj, pairs, times, periodic=True, opt=False)
    c = compute_distances(ptraj[:1], pairs, periodic=True, opt=True)
    d = compute_distances(ptraj[:1], pairs, periodic=True, opt=False)
    eq(a, c)
    eq(b, d)

def _run_amber_traj_t(traj, ext_ref):
    # Test triclinic case where simple approach in Tuckerman text does not
    # always work
    distopt = compute_distances_t(traj, atom_pairs=[[0, 9999]], time_pairs=[[0, 2]], opt=True)
    distslw = compute_distances_t(traj, atom_pairs=[[0, 9999]], time_pairs=[[0, 2]], opt=False)

def test_amber_t(get_fn):
    ext_ref = np.array([17.4835, 22.2418, 24.2910, 22.5505, 12.8686, 22.1090,
                        7.4472, 22.4253, 19.8283, 20.6935]) / 10
    traj = md.load(get_fn('test_good.nc'), top=get_fn('test.parm7'))
    _run_amber_traj_t(traj, ext_ref)
