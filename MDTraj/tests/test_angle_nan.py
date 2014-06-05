##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Lee-Ping Wang
# Contributors: Robert McGibbon
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


from __future__ import print_function, division
import tempfile, os
import numpy as np
import mdtraj as md
from itertools import product
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.testing import get_fn, eq, DocStringFormatTester

def test_angle_nan():

    traj = md.load(get_fn("water2.pdb"))

    angles = compute_angles(traj, np.array([[3, 4, 0]]), opt=True)

    for i in angles[0]:
        assert (not np.isnan(i))

if __name__ == "__main__":
    test_angle_nan()
