##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors: Tim Moore
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
import pytest
from mdtraj.geometry import shape
from mdtraj.testing import eq

"""The trajectories `2lines.pdb` and `3lines.pdb` contain several frames of,
respectively, 2 and 3 "residues" each consisting of two atoms in different
orientations relative to one another.

`2lines.pdb`
Frame 0: || in x         - S2 should be 1.0
Frame 1: || in y         - S2 should be 1.0
Frame 2: || in z         - S2 should be 1.0
Frame 3: |- in x/y       - S2 should be 0.25
Frame 4: |- in y/z       - S2 should be 0.25

`3lines.pdb`
Frame 0: ||| in x        - S2 should be 1.0
Frame 1: ||| in y        - S2 should be 1.0
Frame 2: ||| in z        - S2 should be 1.0
Frame 3: at right angles - S2 should be 0.0
"""

@pytest.fixture()
def traj1(get_fn):
    return md.load(get_fn('1line.pdb'))


@pytest.fixture()
def traj2(get_fn):
    return md.load(get_fn('2lines.pdb'))


@pytest.fixture()
def traj3(get_fn):
    return md.load(get_fn('3lines.pdb'))

@pytest.fixture()
def traj4(get_fn):
    return md.load(get_fn('2koc.pdb'))

def test_gyration(traj1, traj2, traj3):
    assert eq(shape.compute_gyration_tensor(traj1),
              shape._compute_gyration_tensor_slow(traj1))
    assert eq(shape.compute_gyration_tensor(traj2),
              shape._compute_gyration_tensor_slow(traj2))
    assert eq(shape.compute_gyration_tensor(traj3),
              shape._compute_gyration_tensor_slow(traj3))

def test_principal_moments(traj4):
    rg_actual = md.compute_rg(traj4)

    principal_moments = shape.principal_moments(traj4)

    rg_computed = np.sqrt(principal_moments.sum(axis=1))

    assert eq(rg_actual,rg_computed)

def test_asphericity(traj4):
    b_computed = shape.asphericity(traj4)

    pm = shape.principal_moments(traj4)
    rg = md.compute_rg(traj4)
    b_actual = 1.5 * pm[:,2] - rg**2 / 2.0

    assert eq(b_actual,b_computed)

def test_shape_metrics(traj4):
    b = shape.asphericity(traj4)
    c = shape.acylindricity(traj4)
    rg = md.compute_rg(traj4)
    
    kappa_actual = (b**2 + 0.75*c**2)/(rg**4)
    kappa_computed = shape.relative_shape_anisotropy(traj4)

    assert eq(kappa_actual,kappa_computed)

    
