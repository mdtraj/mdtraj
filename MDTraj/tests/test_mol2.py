##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
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

import numpy as np
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj.formats import mol2
doc = DocStringFormatTester(mol2)


def test_load_mol2():
    trj = md.load(get_fn('imatinib.mol2'))
    ref_trj = md.load(get_fn('imatinib.pdb'))
    eq(trj.xyz, ref_trj.xyz)
    
    ref_top, ref_bonds = ref_trj.top.to_dataframe()
    top, bonds = trj.top.to_dataframe()
    eq(bonds, ref_bonds)

  # Third row from mol2 file copied below, used in testing.
  #       3 N1          8.5150   -0.1620    1.3310 n3        1 LIG     -0.732600

def test_mol2_dataframe():
    top, bonds = mol2.mol2_to_dataframes(get_fn("imatinib.mol2"))
    eq(top.name[2], "N1")
    eq(top.atype[2], "n3")
    eq(top.resName[2], "LIG")
    eq(float(top.charge[2]), -0.732600)
