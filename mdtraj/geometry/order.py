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
import pdb

import numpy as np

from mdtraj.geometry.distance import compute_center_of_mass
from mdtraj.utils import ensure_type


__all__ = ['compute_nematic_order']


def compute_nematic_order(traj, indices='chains'):
    """Compute the nematic ordering of a group in every frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute ordering in.
    indices: str or array-like, optional, default='chains'
        The group to consider. Users can pass their own indices or
        Recognized groups are 'chains' and 'residues'.

    Returns
    -------
    s2 : np.ndarray, shape=(traj.n_frames,), dtype=float
        Nematic order parameter values in every frame.

    """

    if indices.lower() == 'chains':
        group = list(traj.top.chains)
    elif indices.lower() == 'residues':
        group = list(traj.top.residues)
    else:
        raise ValueError('Invalid selection: {0}'.format(indices))
    indices = [[at.index for at in compound.atoms]
               for compound in group]

    all_directors = np.zeros(shape=(traj.n_frames, len(indices), 3),
                             dtype=np.float64)
    for i, ids in enumerate(indices):
        sub_traj = traj.atom_slice(ids)
        director = compute_director(sub_traj)
        all_directors[:, i, :] = director
    return all_directors
    Q = compute_Q_tensor(all_directors)

    print(Q)
    s2 = np.zeros(shape=(traj.n_frames,), dtype=np.float64)
    # TODO: vectorize
    for n in range(traj.n_frames):
        w, _ = np.linalg.eig(Q)
        s2[n] = w.max()
    return s2


def compute_Q_tensor(all_directors):
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
    # center_of_mass  shape=(n_frames, n_groups, n_atoms, 3)
    # inertia  shape=(n_frames, n_groups, 3, 3)
    # all_directors  shape=(n_frames, n_groups, 3)
    # Q_tensors  shape=(n_frames, 3, 3)

    all_directors = ensure_type(all_directors, dtype=np.float64, ndim=3,
                                name='directors', shape=(None, None, 3))
    Q = np.zeros(shape=(all_directors.shape[0], 3, 3), dtype=np.float64)

    # TODO: vectorize
    for n, directors in enumerate(all_directors):
        normed = directors / np.sqrt((directors ** 2.0).sum(-1))[..., np.newaxis]
        for vector in normed:
            Q[n, 0, 0] += 3.0 * vector[0] * vector[0] - 1
            Q[n, 0, 1] += 3.0 * vector[0] * vector[1]
            Q[n, 0, 2] += 3.0 * vector[0] * vector[2]
            Q[n, 1, 0] += 3.0 * vector[1] * vector[0]
            Q[n, 1, 1] += 3.0 * vector[1] * vector[1] - 1
            Q[n, 1, 2] += 3.0 * vector[1] * vector[2]
            Q[n, 2, 0] += 3.0 * vector[2] * vector[0]
            Q[n, 2, 1] += 3.0 * vector[2] * vector[1]
            Q[n, 2, 2] += 3.0 * vector[2] * vector[2] - 1
    Q /= (2.0 * normed.shape[0])
    return Q


def compute_director(traj):
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
    inertia_tensor = compute_inertia_tensor(traj)
    directors = np.zeros(shape=(traj.n_frames, 3), dtype=np.float64)
    for n, I in enumerate(inertia_tensor):
        w, v = np.linalg.eig(I)
        directors[n] = v[:, np.argmin(w)]
    return directors


def compute_inertia_tensor(traj):
    """Compute the inertia tensor of a group of atoms. """
    I = np.zeros(shape=(traj.n_frames, 3, 3), dtype=np.float64)
    com = compute_center_of_mass(traj)
    masses = np.array([a.element.mass for a in traj.top.atoms])

    # TODO: vectorize
    for n, xyz in enumerate(traj.xyz):
        for i, coord0 in enumerate(xyz):
            mass = masses[i]
            coord = coord0 - com[n]
            I[n, 0, 0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
            I[n, 1, 1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
            I[n, 2, 2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
            I[n, 0, 1] -= mass * coord[0] * coord[1]
            I[n, 0, 2] -= mass * coord[0] * coord[2]
            I[n, 1, 2] -= mass * coord[1] * coord[2]
        I[n, 1, 0] = I[n, 0, 1]
        I[n, 2, 0] = I[n, 0, 2]
        I[n, 2, 1] = I[n, 1, 2]
    return I
