##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Robert T. McGibbon
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


import tempfile, os
import numpy as np
import mdtraj as md
from mdtraj.formats import GroTrajectoryFile, gro
from mdtraj.testing import get_fn, eq, DocStringFormatTester
TestDocstrings = DocStringFormatTester(gro, error_on_none=True)

fd, temp = tempfile.mkstemp(suffix='.gro')
os.close(fd)

def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_read_write():
    for t in [md.load(get_fn('4waters.pdb')),   # no unit cell
              md.load(get_fn('native2.pdb'))]:  # unit cell
        with GroTrajectoryFile(temp, 'w') as f:
            f.write(t.xyz, t.topology, unitcell_vectors=t.unitcell_vectors)

        with GroTrajectoryFile(temp) as f:
            xyz, time, unitcell = f.read()
            top = f.topology

        eq(xyz, t.xyz, decimal=3)
        eq(list(top.atoms), list(t.top.atoms))
        if t.unitcell_vectors is not None:
            eq(unitcell, t.unitcell_vectors)

def test_load():
    for tref in [md.load(get_fn('4waters.pdb')),   # no unit cell
                 md.load(get_fn('native2.pdb'))]:  # unit cell
        with GroTrajectoryFile(temp, 'w') as f:
            f.write(tref.xyz, tref.topology, unitcell_vectors=tref.unitcell_vectors)

        t = md.load(temp)
        eq(t.xyz, tref.xyz, decimal=3)
        eq(list(t.top.atoms), list(tref.top.atoms))
        eq(t.unitcell_vectors, tref.unitcell_vectors)

