##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christoph Klein
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

from __future__ import print_function

import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq, assert_raises, DocStringFormatTester
from mdtraj.geometry import order


"""The trajectories `2lines.pdb` and `3lines.pdb` contain several frames of,
respectively, 2 and 3 "residues" each consisting of two atoms in different
orientations relative to one another.

`2lines.pdb`
Frame 0: || in x         - S2 should be 1.0
Frame 1: || in y         - S2 should be 1.0
Frame 2: || in z         - S2 should be 1.0
Frame 3: |- in x/y       - S2 should be 0.25
Frame 4: |- in y/z       - S2 should be 0.25

`3lines.pdb`
Frame 0: ||| in x        - S2 should be 1.0
Frame 1: ||| in y        - S2 should be 1.0
Frame 2: ||| in z        - S2 should be 1.0
Frame 3: at right angles - S2 should be 0.0
"""

TRAJ1 = md.load(get_fn('1line.pdb'))
TRAJ2 = md.load(get_fn('2lines.pdb'))
TRAJ2_LINE1 = TRAJ2.atom_slice((0, 1))
TRAJ2_LINE2 = TRAJ2.atom_slice((2, 3))

TRAJ3 = md.load(get_fn('3lines.pdb'))

TestDocstrings = DocStringFormatTester(order, error_on_none=True)


def test_nematic_order():
    assert eq(np.array([1.0, 1.0, 1.0, 0.25, 0.25]),
              order.compute_nematic_order(TRAJ2))
    assert eq(np.array([1.0, 1.0, 1.0, 0.0]),
              order.compute_nematic_order(TRAJ3))


def test_director():
    directors = order._compute_director(TRAJ2_LINE1)
    assert eq(directors,
              np.array([[1, 0, 0],  # Frame 0
                        [0, 1, 0],  # Frame 1
                        [0, 0, 1],  # Frame 2
                        [1, 0, 0],  # Frame 3
                        [1, 0, 0],  # Frame 4
                        ]))
    directors = order._compute_director(TRAJ2_LINE2)
    assert eq(directors,
              np.array([[1, 0, 0],  # Frame 0
                        [0, 1, 0],  # Frame 1
                        [0, 0, 1],  # Frame 2
                        [0, 1, 0],  # Frame 3
                        [0, 0, 1],  # Frame 4
                        ]))


def test_inertia():
    assert eq(order.compute_inertia_tensor(TRAJ1),
              order._compute_inertia_tensor_slow(TRAJ1))
    assert eq(order.compute_inertia_tensor(TRAJ2),
              order._compute_inertia_tensor_slow(TRAJ2))
    assert eq(order.compute_inertia_tensor(TRAJ3),
              order._compute_inertia_tensor_slow(TRAJ3))


def test_nematic_order_args():
    assert_raises(ValueError, lambda: order.compute_nematic_order(TRAJ2, indices='dog'))
    assert_raises(ValueError, lambda: order.compute_nematic_order(TRAJ2, indices='O'))
    assert_raises(ValueError, lambda: order.compute_nematic_order(TRAJ2, indices=1))
    # Input indices with wrong "shapes".
    assert_raises(ValueError, lambda: order.compute_nematic_order(TRAJ2, indices=[[1, [2]], [2]]))
    assert_raises(ValueError, lambda: order.compute_nematic_order(TRAJ2, indices=[1, 2, 3]))

