##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
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


import numpy as np
from mdtraj.geometry import alignment
from mdtraj.testing import eq

xyz1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1.0]])

offset = 1.0 * np.random.normal(size=(3))
rotation = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
xyz2 = rotation.dot(xyz1.T).T + offset
xyz3 = xyz1 + np.random.normal(size=xyz1.shape)
extra_points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]])
xyz4 = rotation.dot(np.concatenate([xyz1, extra_points]).T).T + offset


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


def test_transform2():
    # For testing transform on all points when only subset for align
    xyz1_len = np.shape(xyz1)[0]
    T = alignment.compute_transformation(xyz4[0:xyz1_len, :], xyz1)
    xyz4_prime = T.transform(xyz4)

    eq(xyz1, xyz4_prime[0:xyz1_len, :])
