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

import time
import itertools
import numpy as np
import mdtraj as md

from mdtraj.testing import eq, skipif
from mdtraj.geometry.distance import compute_distances, compute_displacements
from mdtraj.geometry.distance import _displacement_mic, _displacement

N_FRAMES = 20
N_ATOMS = 20

xyz = np.asarray(np.random.randn(N_FRAMES, N_ATOMS, 3), dtype=np.float32)
pairs = np.array(list(itertools.combinations(range(N_ATOMS), 2)), dtype=np.int32)

ptraj = md.Trajectory(xyz=xyz, topology=None)
ptraj.unitcell_vectors = np.ascontiguousarray(np.random.randn(N_FRAMES, 3, 3) + 2*np.eye(3,3), dtype=np.float32)

def test_0():
    a = compute_distances(ptraj, pairs, periodic=False, opt=True)
    b = compute_distances(ptraj, pairs, periodic=False, opt=False)
    eq(a, b)

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
    eq(a, np.sqrt(np.sum(np.square(b), axis=2)))


def test_4():
    # using a really big box, we should get the same results with and without
    # pbcs
    box = np.array([[100, 0, 0], [0, 200, 0], [0, 0, 300]])
    box = np.zeros((N_FRAMES, 3, 3)) + box #broadcast it out
    a = _displacement_mic(xyz, pairs, box)
    b = _displacement(xyz, pairs)
    eq(a, b, decimal=3)

def test_5():
    # simple wrap around along the z axis.
    xyz = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 2.2]]])
    box = np.eye(3,3).reshape(1,3,3)
    result = _displacement_mic(xyz, np.array([[0,1]]), box)
    eq(result, np.array([[[0, 0, 0.2]]]))

