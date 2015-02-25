##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Christoph Klein
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

from __future__ import print_function, division

import numpy as np


__all__ = ['compute_nematic_order']


def compute_nematic_order(traj, select='chains', periodic=True):
    """Compute the nematic ordering of a group in every frame.

    Parameters
    ----------
    traj : Trajectory
        This trajectory is assumed to only contain the atoms to be considered in
        the calculation.
    select: str, optional, default='chains'
        The group to consider. Valid groups are 'chains' and 'residues'. Atoms
        within individual members of the selected group are assumed to be
        consecutively indexed. E.g., a chain will could contain atoms [1, 2, 3]
        but not [1, 2, 4].

    Returns
    -------
    s2 : np.ndarray, shape=(traj.n_frames,), dtype=float
        Nematic order parameter values in every frame.

    """

    # xyz  shape=(n_frames, n_groups, n_atoms, 3)

    if select.lower() == 'chains':
        group = traj.top.chains
    elif select.lower() == 'residues':
        group = traj.top.residues

    indices = list()
    for member in group:
        indices.append()

    xyz = np.split(traj.xyz, indices_or_sections=indices, axis=1)

    # center_of_mass  shape=(n_frames, n_groups, n_atoms, 3)
    # inertia  shape=(n_frames, n_groups, 3, 3)
    # all_directors  shape=(n_frames, n_groups, 3)

    # Q_tensors  shape=(n_frames, 3, 3)

    # S2  shape=(n_frames, )
    s2 = list()
    for Q in Q_tensors:
        w, _ = np.linalg.eig(Q)
        s2.append(w.max())

    return np.asarray(s2)


def _compute_Q_tensor(directors):
    """Compute Q tensor of set of directors.

    Parameters
    ----------
    directors : np.ndarray, shape=(n, 3), dtype=float
        An array of directors.

    Returns
    -------
    Q : dtype=float
        The Q tensor describing the directors.

    See also
    --------
    _compute_director

    References
    ----------

    """
    normed = directors / np.sqrt((directors ** 2.0).sum(-1))[..., np.newaxis]
    Q = np.zeros(shape=(3, 3))
    for vector in normed:
        Q[0, 0] += 3.0 * vector[0] * vector[0] - 1
        Q[0, 1] += 3.0 * vector[0] * vector[1]
        Q[0, 2] += 3.0 * vector[0] * vector[2]
        Q[1, 0] += 3.0 * vector[1] * vector[0]
        Q[1, 1] += 3.0 * vector[1] * vector[1] - 1
        Q[1, 2] += 3.0 * vector[1] * vector[2]
        Q[2, 0] += 3.0 * vector[2] * vector[0]
        Q[2, 1] += 3.0 * vector[2] * vector[1]
        Q[2, 2] += 3.0 * vector[2] * vector[2] - 1
    Q /= (2.0 * normed.shape[0])
    return Q


def _compute_director(I):
    """Compute characteristic vector describing a group's orientation.

    Parameters
    ----------
    I : np.ndarray, shape=(3, 3), dtype=float
        Moment of inertia tensor.

    Returns
    -------
    director :  np.ndarray, shape=(3,), dtype=float
        Characteristic vector.

    References
    ----------

    """
    w, v = np.linalg.eig(I)
    director = v[:, np.argmin(w)]
    return director


def compute_inerita_tensor(xyz, masses):
    """ """
    I = np.zeros(shape=(3, 3))
    com =
    for i, coord0 in enumerate(xyz):
        mass = masses[i]
        coord = coord0 - com
        I[0, 0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
        I[1, 1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
        I[2, 2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
        I[0, 1] -= mass * coord[0] * coord[1]
        I[0, 2] -= mass * coord[0] * coord[2]
        I[1, 2] -= mass * coord[1] * coord[2]
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    return I
