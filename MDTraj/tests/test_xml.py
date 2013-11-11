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

from mdtraj.testing import get_fn, eq
from mdtraj.trajectory import load


def test_0():
    t1 = load(get_fn('native2.xml'), top=get_fn('native2.pdb'))
    t2 = load(get_fn('native2.pdb'))

    t1.center_coordinates()
    t2.center_coordinates()

    yield lambda: eq(t1.xyz, t2.xyz)
    yield lambda: eq(t1.unitcell_vectors, t2.unitcell_vectors)
