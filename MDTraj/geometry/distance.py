# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""Pure python code to calculate bond lengths in a trajectory
"""
##############################################################################
# Imports
##############################################################################

import numpy as np

##############################################################################
# Functions
##############################################################################


def _compute_distances_xyz(xyz, atom_pairs):
    """Compute the bond angles for each frame in xyz

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.

    Returns
    -------
    distances : np.ndarray, shape=[n_frames, num_pairs], dtype=float

    Notes
    -----
    This is a reference implementation in python/numpy.
    """
    delta = np.diff(xyz[:, atom_pairs], axis=2)[:, :, 0]

    distances = (delta ** 2.).sum(-1) ** 0.5

    return distances


def compute_distances(traj, atom_pairs):
    """Compute the bond angles for each frame in traj

    Parameters
    ----------
    tray : Trajectory
        Trajectory to compute distances in
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.

    Returns
    -------
    distances : np.ndarray, shape=[n_frames, num_pairs], dtype=float

    Notes
    -----
    This is a reference implementation in python/numpy.
    """
    return _compute_distances_xyz(traj.xyz, atom_pairs)
