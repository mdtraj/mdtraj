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

import itertools
import os
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pytest

import mdtraj as md
from mdtraj import element
from mdtraj.utils.unitcell import check_valid_unitcell_angles


def pytest_configure(config):
    """This runs before all the tests. Untars the reference files"""
    test_path = Path(__file__).parent
    if not os.path.exists(f"{test_path}/data"):
        with tarfile.open(f"{test_path}/data.tar.gz") as tar:
            tar.extractall(path=f"{test_path}", filter="fully_trusted")


@pytest.fixture(scope="session")
def get_fn():
    test_dir = os.path.dirname(os.path.abspath(__file__))

    def _get_fn(fn):
        return f"{test_dir}/data/{fn}"

    return _get_fn


@pytest.fixture(scope="function")
def gen_random_ptraj(request):
    """Fixture for preparing a test trajectories with random coordinates
    and unitcell
    """
    request.cls.N_FRAMES = N_FRAMES = 20
    request.cls.N_ATOMS = N_ATOMS = 20

    request.cls.rng = rng = np.random.default_rng()

    request.cls.xyz = xyz = np.asarray(rng.standard_normal((N_FRAMES, N_ATOMS, 3), dtype=np.float32))
    request.cls.pairs = np.array(list(itertools.combinations(range(N_ATOMS), 2)), dtype=np.int32)
    request.cls.times = np.array([[i, 0] for i in range(N_FRAMES)[::2]], dtype=np.int32)

    request.cls.ptraj = md.Trajectory(xyz=xyz, topology=None)
    request.cls.ptraj.unitcell_vectors = np.ascontiguousarray(
        rng.standard_normal((N_FRAMES, 3, 3)) + 2 * np.eye(3, 3),
        dtype=np.float32,
    )

    while True:
        # RNG can give us unitcells with zero or imaginary volume.
        # This `while` loop guarantees that we don't get that.
        request.cls.ptraj.unitcell_vectors = np.ascontiguousarray(
            rng.standard_normal((N_FRAMES, 3, 3)) + 2 * np.eye(3, 3),
            dtype=np.float32,
        )

        if np.all([check_valid_unitcell_angles(*row) for row in request.cls.ptraj.unitcell_angles]):
            break


@pytest.fixture(scope="function")
def h5traj(tmp_path):
    xyz = np.around(np.random.randn(10, 5, 3).astype(np.float32), 2)
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    topology.add_atom("CA", element.carbon, residue)
    topology.add_atom("HG1", element.hydrogen, residue)
    topology.add_atom("SG", element.sulfur, residue)
    topology.add_atom("OD1", element.oxygen, residue)
    topology.add_atom("NE", element.nitrogen, residue)

    time = np.arange(10) ** 2
    unitcell_lengths = np.array([[1.1, 1.2, 1.3]] * 10)
    unitcell_angles = np.array([[90, 90, 95]] * 10)

    traj = md.Trajectory(
        xyz,
        topology=topology,
        time=time,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
    )

    fn = f"{tmp_path}/ref.h5"
    traj.save(fn)
    return traj, fn, str(tmp_path)


@pytest.fixture(scope="function")
def h5traj_full_metadata(tmp_path):
    # h5 trajectory with full bond metadata and formal charge
    # for use testing roundtrip of topologies
    # added because fixture above  (h5traj) is used
    # in many tests where full data is not roundtripped.
    # PR2101

    from mdtraj.core.topology import Double, Single

    xyz = np.around(np.random.randn(10, 5, 3).astype(np.float32), 2)
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    ca = topology.add_atom("CA", element.carbon, residue, formal_charge=0)
    hg1 = topology.add_atom("HG1", element.hydrogen, residue, formal_charge=1)
    sg = topology.add_atom("SG", element.sulfur, residue, formal_charge=-1)
    od1 = topology.add_atom("OD1", element.oxygen, residue, formal_charge=-1)
    ne = topology.add_atom("NE", element.nitrogen, residue, formal_charge=1)
    topology.add_bond(ca, hg1, order=1, type=Single)
    topology.add_bond(ca, sg, order=2, type=Double)
    topology.add_bond(od1, ne, order=1, type=Single)

    time = np.arange(10) ** 2
    unitcell_lengths = np.array([[1.1, 1.2, 1.3]] * 10)
    unitcell_angles = np.array([[90, 90, 95]] * 10)

    traj = md.Trajectory(
        xyz,
        topology=topology,
        time=time,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
    )

    fn = f"{tmp_path}/ref.h5"
    traj.save(fn)
    return traj, fn, str(tmp_path)


@pytest.fixture(scope="session")
def vcr_config(request):
    """This runs during pytest-recording setup. Untars the cassettes.

    To remake the cassettes, comment out this function. Then run
    ```
    pytest --record-mode=all test_dssp.py test_pdb.py
    ```

    Delete the old tar and tar up the new folder
    ```
    rm -rf cassettes.tar.gz
    tar czf cassettes.tar.gz cassettes
    ```

    Then, un-comment out this function again.
    """
    test_path = Path(__file__).parent
    with tarfile.open(f"{test_path}/cassettes.tar.gz") as tar:
        tar.extractall(path=f"{test_path}", filter="fully_trusted")

    def remove_tar():
        if os.path.exists(f"{test_path}/cassettes"):
            shutil.rmtree(f"{test_path}/cassettes")

    request.addfinalizer(remove_tar)
    return {}
