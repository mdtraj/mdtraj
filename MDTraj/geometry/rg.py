##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A Beauchamp
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

"""Pure python code to calculate the radius of gyration of frames in a
trajectory
"""
##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np

__all__ = ['compute_rg']

##############################################################################
# Functions
##############################################################################


def _compute_rg_xyz(xyz, masses=None):
    """Compute the radius of gyration for every frame.

    Parameters
    ----------
    xyz : ndarray
        xyz coordinates
    masses : ndarray
        Transition matrix

    Returns
    -------
    Rg : ndarray
        Rg for every frame

    Notes
    -----
    If masses are none, assumes equal masses.
    """
    traj_length, num_atoms, num_dims = xyz.shape
    if not num_dims == 3:
        raise ValueError("What did you pass me?")
    if not xyz.dtype == np.float32:
        xyz = np.float32(xyz)
    if masses is None:
        masses = np.ones(num_atoms)

    weights = masses / masses.sum()

    mu = xyz.mean(1)
    centered = (xyz.transpose((1, 0, 2)) - mu).transpose((1, 0, 2))
    squared_dists = (centered ** 2).sum(2)
    Rg = (squared_dists * weights).sum(1) ** 0.5

    return Rg


def compute_rg(traj, masses=None):
    """Compute the radius of gyration for every frame.

    Parameters
    ----------
    traj : Trajectory

    masses : ndarray, optional
        array of atom masses.

    Returns
    -------
    Rg : ndarray
        Rg for every frame

    Notes
    -----
    If masses are none, assumes equal masses.
    """
    return _compute_rg_xyz(traj.xyz, masses=masses)
