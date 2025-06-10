##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Jeremy MG Leung
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
import pytest

from mdtraj.utils.unitcell import check_valid_unitcell_angles


class Test_Unitcell:
    @pytest.mark.parametrize(["alpha", "beta", "gamma"], [[90, 90, 90], [70, 90, 30]], ids=["orthogonal", "triclinic"])
    def test_valid(self, alpha, beta, gamma):
        """Test for valid unitcell angles, single frame"""
        unitcell_angles = np.array([alpha, beta, gamma], dtype=np.float32)

        assert check_valid_unitcell_angles(*unitcell_angles)  # True

    @pytest.mark.parametrize(
        ["alpha", "beta", "gamma"],
        [[126.96979, 54.83003, 178.20027], [44.82802, 110.52171, 65.693695]],
        ids=["over360", "close"],
    )
    def test_invalid_over360(self, alpha, beta, gamma):
        """Test for invalid unitcell angles, single frame"""
        unitcell_angles = np.array([alpha, beta, gamma], dtype=np.float32)

        assert not check_valid_unitcell_angles(*unitcell_angles)  # False

    def test_valid_array(self):
        """test for valid unitcell angles, input being an array"""
        unitcell_angles = np.array([[90.0, 90.0, 90], [70, 90, 30], [40, 20, 30]], dtype=np.float32)

        assert np.all(check_valid_unitcell_angles(*unitcell_angles.T) == [True, True, True])

    def test_invalid_over360_array(self):
        """Test for invalid unitcell angles that sum up to over 360˚, input being an array"""
        unitcell_angles = np.array([[90, 90, 90], [126.96979, 54.83003, 178.20027], [70, 90, 30]], dtype=np.float32)

        assert np.all(check_valid_unitcell_angles(*unitcell_angles.T) == [True, False, True])

    def test_invalid_angles_close_array(self):
        """Test for invalid unitcell angles where angles 1 and 3 sum up to more than angle 2"""
        unitcell_angles = np.array([[44.82802, 110.52171, 65.693695], [90, 90, 90], [70, 90, 30]], dtype=np.float64)

        assert np.all(check_valid_unitcell_angles(*unitcell_angles.T) == [False, True, True])
