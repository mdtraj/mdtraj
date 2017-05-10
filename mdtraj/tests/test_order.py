##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Tim Moore
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
from mdtraj.testing import get_fn, eq, skipif, assert_allclose
from mdtraj.utils import six


def test_order():
    """Made a perfectly aligned monolayer, should have S2 = 1
    """
    traj = md.load(get_fn('alkane-monolayer.pdb'))
    indices = [list(range(1900+x, 1900+x+36)) for x in range(0, 64*36, 36)]
    s2 = md.compute_nematic_order(traj, indices=indices)
    assert_allclose(1.0, s2)


def test_directors():
    """All chains in example system aligned along z axis
    """
    traj = md.load(get_fn('alkane-monolayer.pdb'))
    # only use the carbons for this calculation since they are excactly aligned
    # with the z-axis, and the hydrogens can throw things off just a bit
    indices = [list(range(1900+x, 1900+x+30, 3)) for x in range(0, 64*36, 36)]
    directors = md.compute_directors(traj, indices=indices)
    for axis, target_angle in zip(np.eye(3), [90.0, 90.0, 0.0]):
        dotproducts = np.tensordot(directors, axis, axes=1)
        angles = np.degrees(np.arccos(dotproducts))
        angles = np.minimum(angles, 180-angles)  # make vector point in positive z
        assert_allclose(target_angle, angles)
