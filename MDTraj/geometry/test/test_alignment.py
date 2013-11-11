##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
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


from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj.geometry import alignment
import numpy as np

xyz1 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1.0]])

offset = 1.0 * np.random.normal(size=(3))
rotation = np.array([[1,0,0],[0,0,-1],[0,1,0]])
xyz2 = rotation.dot(xyz1.T).T + offset
xyz3 = xyz1 + np.random.normal(size=xyz1.shape)

def test_rmsd_zero():
    rmsd_kabsch = alignment.rmsd_kabsch(xyz1, xyz2)
    rmsd_qcp = alignment.rmsd_qcp(xyz1, xyz2)
    eq(float(rmsd_kabsch), 0.0, decimal=5)
    eq(float(rmsd_qcp), 0.0, decimal=5)

def test_rmsd_nonzero():
    rmsd_kabsch = alignment.rmsd_kabsch(xyz1, xyz3)
    rmsd_qcp = alignment.rmsd_qcp(xyz1, xyz3)
    eq(rmsd_kabsch, rmsd_qcp, decimal=5)

def test_transform():
    T = alignment.compute_transformation(xyz2, xyz1)
    xyz2_prime = T.transform(xyz2)

    eq(xyz1, xyz2_prime)
