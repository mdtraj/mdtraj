##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
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


import mdtraj as md
import numpy as np
from mdtraj.formats import ArcTrajectoryFile
from mdtraj.formats import PDBTrajectoryFile
from mdtraj.testing import eq


def test_read_0(get_fn):
    with ArcTrajectoryFile(get_fn('4waters.arc')) as f:
        xyz, leng, ang = f.read()
    with PDBTrajectoryFile(get_fn('4waters.pdb')) as f:
        xyz2 = f.positions
    eq(xyz, xyz2, decimal=3)

def test_read_arctraj(get_fn):
    traj = md.load(get_fn('nitrogen.arc'), top=get_fn('nitrogen.pdb'))
    owntop = md.load(get_fn('nitrogen.arc'))
    eq(traj.xyz, owntop.xyz)

def test_read_pbc(get_fn):
    with ArcTrajectoryFile(get_fn('thf200.arc')) as f:
       xyz, leng, ang = f.read()
    eq(xyz[0, 0], np.array([21.231166, 32.549819, 8.278454]), decimal=3)
    eq(xyz[0,-1], np.array([29.013880, 6.255754, 44.074519]), decimal=3)
    eq(leng, np.array([[45.119376, 45.119376, 45.119376]]), decimal=3)
    eq(ang, np.array([[90.0, 90.0, 90.0]]), decimal=3)
