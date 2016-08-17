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


"""
Tests of generic loading functionality.
"""

from mdtraj import load
from mdtraj.testing import get_fn

def test_load_single():
    """
    Just check for any raised errors coming from loading a single file.
    """
    load(get_fn('frame0.pdb'))

def test_load_single_list():
    """
    See if a single-element list of files is successfully loaded.
    """
    load([get_fn('frame0.pdb')])

def test_load_many_list():
    """
    See if a multi-element list of files is successfully loaded.
    """
    traj = load(2 * [get_fn('frame0.pdb')], discard_overlapping_frames=False)
    assert traj.n_frames == 2
