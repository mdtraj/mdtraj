##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Jason Swails
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

import mdtraj as md
import numpy as np
from mdtraj.testing import get_fn, assert_warns
import warnings

def test_dummy_pdb_box_detection():
    assert_warns(UserWarning, lambda: md.load(get_fn('2koc.pdb')))
    warnings.filterwarnings('ignore', category=UserWarning)
    traj = md.load(get_fn('2koc.pdb'))
    assert traj.unitcell_lengths is None, 'Expected dummy box to be deleted'
    warnings.filterwarnings('default', category=UserWarning)
