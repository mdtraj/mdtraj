##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Matthew Harrigan
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

import os
import itertools

import pytest
import numpy as np

import mdtraj as md

flaky_pdb_dl = pytest.mark.flaky(rerun=3, reason='github-node flaky pdb dl')


@pytest.fixture(scope="session")
def get_fn():
    test_dir = os.path.dirname(os.path.abspath(__file__))

    def _get_fn(fn):
        return f"{test_dir}/data/{fn}"

    return _get_fn


@pytest.fixture
def gen_random_ptraj(request):
    """
    Fixture for preparing a test trajectories with random coordinates 
    and unitcell

    """
    request.cls.N_FRAMES = N_FRAMES = 20
    request.cls.N_ATOMS = N_ATOMS = 20
    
    request.cls.rng = rng = np.random.default_rng()
    
    request.cls.xyz = xyz = np.asarray(rng.standard_normal((N_FRAMES, N_ATOMS, 3), dtype=np.float32))
    request.cls.pairs = pairs = np.array(list(itertools.combinations(range(N_ATOMS), 2)), dtype=np.int32)
    request.cls.times = times = np.array([[i, 0] for i in range(N_FRAMES)[::2]], dtype=np.int32)
    
    request.cls.ptraj = md.Trajectory(xyz=xyz, topology=None)
    request.cls.ptraj.unitcell_vectors = np.ascontiguousarray(
        rng.standard_normal((N_FRAMES, 3, 3)) + 2 * np.eye(3, 3),
        dtype=np.float32,
    )
