##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors: Tim Moore
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
import pytest
from mdtraj.geometry import order
from mdtraj.testing import eq

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


# TRAJ1 = md.load(get_fn('1line.pdb'))
# TRAJ2 = md.load(get_fn('2lines.pdb'))
# TRAJ2_LINE1 = TRAJ2.atom_slice((0, 1))
# TRAJ2_LINE2 = TRAJ2.atom_slice((2, 3))

# TRAJ3 = md.load(get_fn('3lines.pdb'))

@pytest.fixture()
def traj1(get_fn):
    return md.load(get_fn('1line.pdb'))


@pytest.fixture()
def traj2(get_fn):
    return md.load(get_fn('2lines.pdb'))


@pytest.fixture()
def traj3(get_fn):
    return md.load(get_fn('3lines.pdb'))


def test_nematic_order(traj2, traj3):
    assert eq(np.array([1.0, 1.0, 1.0, 0.25, 0.25]),
              order.compute_nematic_order(traj2))
    assert eq(np.array([1.0, 1.0, 1.0, 0.0]),
              order.compute_nematic_order(traj3))


def test_director(traj2):
    traj2_line1 = traj2.atom_slice((0, 1))
    traj2_line2 = traj2.atom_slice((2, 3))
    directors = order._compute_director(traj2_line1)
    assert eq(directors,
              np.array([[1, 0, 0],  # Frame 0
                        [0, 1, 0],  # Frame 1
                        [0, 0, 1],  # Frame 2
                        [1, 0, 0],  # Frame 3
                        [1, 0, 0],  # Frame 4
                        ]))
    directors = order._compute_director(traj2_line2)
    assert eq(directors,
              np.array([[1, 0, 0],  # Frame 0
                        [0, 1, 0],  # Frame 1
                        [0, 0, 1],  # Frame 2
                        [0, 1, 0],  # Frame 3
                        [0, 0, 1],  # Frame 4
                        ]))


def test_inertia(traj1, traj2, traj3):
    assert eq(order.compute_inertia_tensor(traj1),
              order._compute_inertia_tensor_slow(traj1))
    assert eq(order.compute_inertia_tensor(traj2),
              order._compute_inertia_tensor_slow(traj2))
    assert eq(order.compute_inertia_tensor(traj3),
              order._compute_inertia_tensor_slow(traj3))

              
def test_nematic_order_args(traj2):
    with pytest.raises(ValueError):
        order.compute_nematic_order(traj2, indices='dog')
    with pytest.raises(ValueError):
        order.compute_nematic_order(traj2, indices='O')
    with pytest.raises(ValueError):
        order.compute_nematic_order(traj2, indices=1)

    # Input indices with wrong "shapes".
    with pytest.raises(ValueError):
        order.compute_nematic_order(traj2, indices=[[1, [2]], [2]])
    with pytest.raises(ValueError):
        order.compute_nematic_order(traj2, indices=[1, 2, 3])


def test_order_from_traj(get_fn):
    """Made a perfectly aligned monolayer, should have S2 = 1
    """
    traj = md.load(get_fn('alkane-monolayer.pdb'))
    indices = [list(range(1900 + x, 1900 + x + 36)) for x in range(0, 64 * 36, 36)]
    s2 = md.compute_nematic_order(traj, indices=indices)
    np.testing.assert_allclose(1.0, s2)


def test_directors_angle_from_traj(get_fn):
    """All chains in example system aligned along z axis
    """
    traj = md.load(get_fn('alkane-monolayer.pdb'))
    # only use the carbons for this calculation since they are excactly aligned
    # with the z-axis, and the hydrogens can throw things off just a bit
    indices = [list(range(1900 + x, 1900 + x + 30, 3)) for x in range(0, 64 * 36, 36)]
    directors = md.compute_directors(traj, indices=indices)
    for axis, target_angle in zip(np.eye(3), [90.0, 90.0, 0.0]):
        dotproducts = np.tensordot(directors, axis, axes=1)
        angles = np.degrees(np.arccos(dotproducts))
        angles = np.minimum(angles, 180 - angles)  # make vector point in positive z
        np.testing.assert_allclose(target_angle, angles)
