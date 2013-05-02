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
"""Pure python code to calculate bond angles in a trajectory
"""
##############################################################################
# Imports
##############################################################################

import numpy as np

##############################################################################
# Functions
##############################################################################


def _compute_bond_angles_xyz(xyz, angle_indices):
    """Compute the bond angles for each frame in xyz

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    angle_indices : np.ndarray, shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle

    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float

    Notes
    -----
    This is a reference single threaded implementation in python/numpy.
    """

    n_frames = xyz.shape[0]
    angles = np.zeros((n_frames, len(angle_indices)))

    for i in xrange(n_frames):
        for j, (m, o, n) in enumerate(angle_indices):
            u_prime = xyz[i, m, :] - xyz[i, o, :]
            v_prime = xyz[i, n, :] - xyz[i, o, :]
            u_norm = np.linalg.norm(u_prime)
            v_norm = np.linalg.norm(v_prime)

            angles[i, j] = np.arccos(np.dot(u_prime, v_prime) /
                                    (u_norm * v_norm))

    return angles


def compute_bond_angles(traj, angle_indices):
    """Compute the bond angles for each frame in traj

    This is a reference single threaded implementation in python/numpy

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    angle_indices : np.ndarray, shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle

    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float
    """
    return _compute_bond_angles_xyz(traj.xyz, angle_indices)
