# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj.geometry import alignment
import numpy as np

xyz1 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1.0]])

offset = 1.0 * np.random.normal(size=(3))
rotation = np.array([[1,0,0],[0,0,-1],[0,1,0]])
xyz2 = rotation.dot(xyz1.T).T + offset
xyz3 = xyz1 + np.random.normal(size=xyz1.shape)

def test_rmsd_zero():
    rmsd_kabsch = alignment.rmsd(xyz1, xyz2)
    rmsd_qcp = alignment.rmsd_qcp(xyz1, xyz2)
    eq(rmsd_kabsch, 0.0, decimal=5)
    eq(rmsd_qcp, 0.0, decimal=5)

def test_rmsd_nonzero():
    rmsd_kabsch = alignment.rmsd(xyz1, xyz3)
    rmsd_qcp = alignment.rmsd_qcp(xyz1, xyz3)
    eq(rmsd_kabsch, rmsd_qcp, decimal=5)
