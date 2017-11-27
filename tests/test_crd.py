##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: TJ Lane
# Contributors: Robert McGibbon and Jason Swails
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
from mdtraj.testing import eq

def test_prmtop_with_overflow_mdcrd(get_fn):
    traj = md.load(get_fn('overflow.crd'), top=get_fn('tz2.truncoct.parm7'))
    # Check that the coordinates that 'overflow' were parsed correctly
    check = np.asarray([-0.131, -221.878, 2.631]) / 10.0 # To nanometers
    eq(check, traj.xyz[-1][-7])
