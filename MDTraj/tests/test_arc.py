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


import tempfile, os
import numpy as np
import mdtraj as md
from mdtraj import ArcTrajectoryFile, arc, PDBTrajectoryFile
from mdtraj.testing import get_fn, eq, DocStringFormatTester
TestDocstrings = DocStringFormatTester(arc, error_on_none=True)

fd, temp = tempfile.mkstemp(suffix='.arc')
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.close(fd)
    os.unlink(temp)


def test_read_0():
    with ArcTrajectoryFile(get_fn('4waters.arc')) as f:
        xyz = f.read()
    with PDBTrajectoryFile(get_fn('4waters.pdb')) as f:
        xyz2 = f.positions
    eq(xyz, xyz2, decimal=3)
