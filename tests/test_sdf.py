##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2026 Stanford University and the Authors
#
# Authors: Stefan Hervø-Hansen
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

import mdtraj as md
from mdtraj.geometry import compute_sdf


@pytest.fixture
def test_traj():
    """Construct a deterministic trajectory for SDF validation."""
    box_size = 4.0
    eps = 1e-6

    # 8 corners + 1 center
    corner_positions = np.array(
        [
            [0, 0, 0],
            [4 - eps, 0, 0],
            [0, 4 - eps, 0],
            [0, 0, 4 - eps],
            [4 - eps, 4 - eps, 0],
            [4 - eps, 0, 4 - eps],
            [0, 4 - eps, 4 - eps],
            [4 - eps, 4 - eps, 4 - eps],
        ]
    )
    center_position = np.array([[2, 2, 2]])

    n_frames = 9
    positions = np.zeros((n_frames, 9, 3))
    for frame in range(n_frames):
        positions[frame, :8] = corner_positions
        positions[frame, 8] = center_position
        positions[frame, :frame] = center_position  # move progressively to center

    # Build topology
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue("TST", chain)
    for i in range(9):
        top.add_atom(f"A{i}", element=md.element.carbon, residue=residue)

    traj = md.Trajectory(positions, top)
    traj.unitcell_lengths = np.array([[box_size, box_size, box_size]] * n_frames)
    traj.unitcell_angles = np.array([[90.0, 90.0, 90.0]] * n_frames)
    return traj


def test_sdf_center_and_corners(test_traj):
    """Test SDF normalization and voxel assignment."""
    particles = test_traj.topology.select("all")

    # Updated to unpack grid
    sdf, grid = compute_sdf(
        test_traj,
        solute_indices=particles,
        solvent_indices=particles,
        grid_spacing=1.0,
        pre_centered=True,
        return_reference=False,  # make sure we only get sdf + grid
    )

    # --- Analytical expectations ---
    V = 4.0**3
    N = 9
    rho = N / V
    voxel_volume = grid["dx"] * grid["dy"] * grid["dz"]
    n_frames = 9

    normalization = rho * voxel_volume * n_frames

    # Center voxel (2,2,2): sum(1..9) = 45
    expected_center = 45 / normalization
    np.testing.assert_allclose(
        sdf[2, 2, 2],
        expected_center,
        rtol=1e-6,
        err_msg="Center voxel does not match analytical value",
    )

    # Corner voxels
    corner_indices = [
        (0, 0, 0),
        (3, 0, 0),
        (0, 3, 0),
        (0, 0, 3),
        (3, 3, 0),
        (3, 0, 3),
        (0, 3, 3),
        (3, 3, 3),
    ]
    expected_counts = np.arange(1, 9)
    for idx, count in zip(corner_indices, expected_counts):
        expected_value = count / normalization
        np.testing.assert_allclose(
            sdf[idx],
            expected_value,
            rtol=1e-6,
            err_msg=f"Voxel {idx} incorrect",
        )


def test_sdf_shape(test_traj):
    """Ensure grid shape matches expected dimensions."""
    particles = test_traj.topology.select("all")

    sdf, grid = compute_sdf(
        test_traj,
        solute_indices=particles,
        solvent_indices=particles,
        grid_spacing=1.0,
        pre_centered=True,
        return_reference=False,
    )

    expected_shape = (grid["nx"], grid["ny"], grid["nz"])
    assert sdf.shape == expected_shape, "Grid shape is incorrect"


def test_sdf_non_negative(test_traj):
    """SDF should never contain negative values."""
    particles = test_traj.topology.select("all")

    sdf, grid = compute_sdf(
        test_traj,
        solute_indices=particles,
        solvent_indices=particles,
        grid_spacing=1.0,
        pre_centered=True,
        return_reference=False,
    )

    assert np.all(sdf >= 0), "SDF contains negative values"
