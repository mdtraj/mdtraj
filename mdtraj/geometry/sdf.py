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

from mdtraj.utils import ensure_type

__all__ = ["compute_sdf"]


def compute_sdf(
    traj,
    solute_indices,
    solvent_indices,
    grid_spacing=0.1,
    pre_centered=False,
    filename=None,
    return_reference=False,
):
    """
    Compute the spatial distribution function (SDF) of solvent around a solute.

    This function calculates a normalized 3D density map of solvent positions
    relative to a solute molecule or set of solute atoms. The SDF is computed
    on a voxel grid covering the simulation box and can optionally be written
    to a Gaussian cube file for visualization.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory containing solute and solvent atoms.
    solute_indices : array_like of int
        Atom indices of solute molecules.
    solvent_indices : array_like of int
        Atom indices of solvent molecules.
    grid_spacing : float, optional
        Grid spacing for voxelization in nanometers (default: 0.1 nm).
    pre_centered : bool, optional
        If True, assumes the solute is already centered; otherwise, the solute
        will be centered in the box using periodic imaging.
    filename : str, optional
        If provided, writes the computed SDF to a Gaussian cube file.
    return_reference : bool, optional
        If True, also return the (possibly centered) trajectory.

    Returns
    -------
    sdf : np.ndarray, shape (nx, ny, nz)
        3D array of normalized SDF values.
    grid : dict
        Dictionary containing voxel grid metadata:
            'nx', 'ny', 'nz' : int, number of voxels along each axis
            'dx', 'dy', 'dz' : float, voxel spacing in nm
            'box_lengths' : ndarray, maximum box lengths used for the grid
    traj : mdtraj.Trajectory, optional
        Returned only if `return_reference=True`; the trajectory with solute
        centered if `pre_centered=False`.

    Examples
    --------
    >>> import mdtraj as md
    >>> traj = md.load_xtc('trajectory.xtc', top='protein_in_water.pdb')
    >>> protein = traj.topology.select('protein')
    >>> water = traj.topology.select('water')
    >>> sdf, grid = compute_sdf(traj, solute=protein, solvent=water, grid_spacing=0.1)
    >>> # Write cube file
    >>> compute_sdf(traj, solute=protein, solvent=water, grid_spacing=0.1,
    ...             filename='sdf_output.cube')
    """

    # Check types
    solute_indices = ensure_type(solute_indices, dtype=int, ndim=1, name="solute_indices")
    solvent_indices = ensure_type(solvent_indices, dtype=int, ndim=1, name="solvent_indices")

    if not pre_centered:
        traj = _centering_solute(traj, solute_indices)

    sdf, grid = _sdf(traj, grid_spacing, solvent_indices)

    if filename is not None:
        _write_cube_file(filename, traj, sdf, grid, solute_indices, solvent_indices)

    if return_reference:
        return sdf, grid, traj
    else:
        return sdf, grid


def _sdf(traj, grid_spacing, solvent_indices):
    """
    Internal function: Compute the spatial distribution function (SDF) on a voxel grid.

    This function bins solvent atoms into a 3D grid and normalizes the counts
    relative to an ideal gas density. Periodic boundary conditions are applied
    using each frame's box dimensions.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory containing solvent molecules and box information.
    grid_spacing : float
        Resolution of the voxel grid in nanometers.
    solvent_indices : array_like of int
        Atom indices corresponding to solvent molecules.

    Returns
    -------
    sdf_grid : np.ndarray, shape (nx, ny, nz)
        Normalized SDF grid.
    grid : dict
        Grid metadata:
            'nx','ny','nz' : int, number of voxels along each axis
            'dx','dy','dz' : float, voxel spacing in nm
            'box_lengths'  : ndarray, maximum box lengths used for grid
    """

    # Simulation box dimensions (in nm)
    box_lengths_all = traj.unitcell_lengths
    box_lengths = box_lengths_all.max(axis=0)  # for grid size

    # Voxel grid setup
    nx = int(np.ceil(box_lengths[0] / grid_spacing))
    ny = int(np.ceil(box_lengths[1] / grid_spacing))
    nz = int(np.ceil(box_lengths[2] / grid_spacing))

    dx = box_lengths[0] / nx
    dy = box_lengths[1] / ny
    dz = box_lengths[2] / nz

    # Create a 3D grid for counting
    grid_shape = (nx, ny, nz)
    sdf_grid = np.zeros(grid_shape, dtype=np.float32)

    # Iterate over trajectory frames to compute density
    for i in range(traj.n_frames):
        solvent_coords = traj.xyz[i, solvent_indices, :]

        # Use per-frame box for wrapping
        box_lengths_frame = box_lengths_all[i]
        solvent_coords = solvent_coords % box_lengths_frame

        # Map solvent atoms to the voxel grid
        indices = np.floor(solvent_coords / np.array([dx, dy, dz])).astype(int)
        indices = np.minimum(indices, np.array(grid_shape) - 1)
        valid = (
            (indices[:, 0] >= 0)
            & (indices[:, 0] < grid_shape[0])
            & (indices[:, 1] >= 0)
            & (indices[:, 1] < grid_shape[1])
            & (indices[:, 2] >= 0)
            & (indices[:, 2] < grid_shape[2])
        )
        np.add.at(
            sdf_grid,
            (indices[valid, 0], indices[valid, 1], indices[valid, 2]),
            1,
        )

    # Normalize SDF grid by the ideal gas density and number of frames
    volumes = np.prod(box_lengths_all, axis=1)
    mean_volume = volumes.mean()
    number_density = len(solvent_indices) / mean_volume
    voxel_volume = grid_spacing**3
    sdf_grid /= number_density * voxel_volume * len(traj)

    # Return both SDF and grid info
    grid = {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dx": dx,
        "dy": dy,
        "dz": dz,
        "box_lengths": box_lengths,
    }
    return sdf_grid, grid


def _write_cube_file(filename, traj, sdf, grid, solute_indices, solvent_indices):
    """
    Internal function: Write the SDF to a Gaussian cube file for visualization.

    The SDF grid and solute atom positions are written in Bohr units.
    The output cube file is compatible with VMD, PyMOL, and other tools.

    Parameters
    ----------
    filename : str
        Output cube file path.
    traj : mdtraj.Trajectory
        Trajectory containing solute topology and coordinates.
    sdf : np.ndarray
        3D SDF grid (shape nx x ny x nz).
    grid : dict
        Dictionary containing grid info:
            'nx','ny','nz','dx','dy','dz','box_lengths'
    solute_indices : array_like of int
        Atom indices of solute molecules.
    solvent_indices : array_like of int
        Atom indices of solvent molecules (not used in cube writing).

    Returns
    -------
    None
        Writes the cube file; does not return a value.
    """

    NM_TO_BOHR = 18.8972612478289694072
    nx, ny, nz = grid["nx"], grid["ny"], grid["nz"]
    dx, dy, dz = grid["dx"], grid["dy"], grid["dz"]

    # Convert spacing from nm to Bohr
    dx_bohr = dx * NM_TO_BOHR
    dy_bohr = dy * NM_TO_BOHR
    dz_bohr = dz * NM_TO_BOHR

    # Cube origin: offset by half voxel to center grid points
    origin_bohr = np.array([0.5 * dx_bohr, 0.5 * dy_bohr, 0.5 * dz_bohr])

    # Write the SDF grid to a cube file
    with open(filename, "w") as cube_file:
        # Write cube header
        cube_file.write("Generated with MDTraj\n")
        cube_file.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")

        # Number of atoms + origin
        cube_file.write(
            f"{len(solute_indices):5} {origin_bohr[0]:12.6f} {origin_bohr[1]:12.6f} {origin_bohr[2]:12.6f}\n",
        )

        # Grid dimensions and spacing vectors
        cube_file.write(f"{nx:5} {dx_bohr:12.6f} 0.0 0.0\n")
        cube_file.write(f"{ny:5} 0.0 {dy_bohr:12.6f} 0.0\n")
        cube_file.write(f"{nz:5} 0.0 0.0 {dz_bohr:12.6f}\n")

        # Write atomic coordinates of solute in Bohr
        for idx in solute_indices:
            atom = traj.topology.atom(idx)
            atomic_number = atom.element.atomic_number if atom.element else 0
            x, y, z = traj.xyz[0, atom.index] * NM_TO_BOHR
            cube_file.write(f"{atomic_number:5} {0.0:12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        # Write the SDF grid data (flattened in C-order)
        sdf_flat = sdf.ravel(order="C")
        for i, value in enumerate(sdf_flat):
            cube_file.write(f"{value:13.5e}")
            if (i + 1) % 6 == 0:
                cube_file.write("\n")
        if len(sdf_flat) % 6 != 0:
            cube_file.write("\n")

    return None


def _centering_solute(traj, solute_indices):
    """
    Internal function: Center the solute in the simulation box using periodic imaging.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory containing the solute and solvent molecules.
    solute_indices : array_like of int
        Atom indices of solute molecules to center.

    Returns
    -------
    traj_centered : mdtraj.Trajectory
        New trajectory where the solute is centered within the periodic box.
    """

    anchor_molecules = [
        set(res.atoms) for res in traj.topology.residues if any(atom.index in solute_indices for atom in res.atoms)
    ]
    traj_centered = traj.image_molecules(anchor_molecules=anchor_molecules)
    return traj_centered
