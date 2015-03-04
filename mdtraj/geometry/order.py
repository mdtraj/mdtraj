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
from mdtraj.utils.six import string_types


__all__ = ['compute_nematic_order']  #, 'compute_inertia_tensor']  # Expose this?


def compute_nematic_order(traj, indices='chains'):
    """Compute the nematic ordering of a group in every frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute ordering in.
    indices: str or list of lists, optional, default='chains'
        The group to consider. Users can pass their own indices as a list of
        lists with the "shape" (n_compounds, len(each_compound)).
        Recognized string keywords are 'chains' and 'residues'.

    Returns
    -------
    S2 : np.ndarray, shape=(traj.n_frames,), dtype=float64
        Nematic order parameter values in every frame.

    References
    ----------
    .. [1] Allen, M. P.; Tildesley , D. J. (1987), "Computer Simulation of
           Liquids", Ch. 11.5
    .. [2] http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
    .. [3] http://cmt.dur.ac.uk/sjc/thesis_dlc/node19.html

    """
    if isinstance(indices, string_types):
        if indices.lower() == 'chains':
            group = list(traj.top.chains)
        elif indices.lower() == 'residues':
            group = list(traj.top.residues)
        else:
            raise ValueError('Invalid selection: {0}'.format(indices))

        indices = np.array([[at.index for at in compound.atoms]
                               for compound in group], dtype=np.int32)
    else:
        # TODO: Clean way to ensure that indices is a list of lists of ints?
        # This may be easier than I'm thinking but `ensure_type` won't work
        # since sub-lists of variable lengths should be allowed.
        # E.g. [[1, 2], [3, 4, 5], [6, 7]] should be valid.
        if isinstance(indices, (list, tuple)):
            for sublist in indices:
                if not isinstance(indices[0], (list, tuple)):
                    raise ValueError('Invalid selection: {0}'.format(indices))
                for index in sublist:
                    if not isinstance(index, int):
                        raise ValueError('Invalid selection: {0}'.format(indices))
        else:
            raise ValueError('Invalid selection: {0}'.format(indices))

    # Compute the directors for each compound for each frame.
    all_directors = np.zeros(shape=(traj.n_frames, len(indices), 3),
                             dtype=np.float64)
    for i, ids in enumerate(indices):
        sub_traj = traj.atom_slice(ids)
        director = _compute_nematic_director(sub_traj)
        all_directors[:, i, :] = director

    # From the directors, compute the Q-tensor and nematic order parameter, S2.
    Q_ab = _compute_Q_tensor(all_directors)
    w = np.linalg.eigvals(Q_ab)
    S2 = w.max(axis=1)
    return S2


def _compute_Q_tensor(all_directors):
    """Compute the Q-tensor for a set of directors.

    For each frame,
        Q_{ab} = 1/(2N) sum_{i_molecules} (3 * e_{ia} * e_{ib} - d_{ab})    [1]

    Parameters
    ----------
    directors : np.ndarray, shape=(n_frames, n_compounds, 3), dtype=float64
        An array of directors describing each compound's orientation over time.

    Returns
    -------
    Q_ab : np.ndarray, shape=(traj.n_frames, 3, 3), dtype=float64
        The Q-tensors describing the directors for each frame.

    See also
    --------
    _compute_nematic_director

    References
    ----------
    .. [1] Allen, M. P.; Tildesley , D. J. (1987), "Computer Simulation of
           Liquids", p. 305, Eq. 11.19

    """

    all_directors = ensure_type(all_directors, dtype=np.float64, ndim=3,
                                name='directors', shape=(None, None, 3))
    Q_ab = np.zeros(shape=(all_directors.shape[0], 3, 3), dtype=np.float64)

    # TODO: Look into vectorizing.
    for n, directors in enumerate(all_directors):
        normed = directors / np.sqrt((directors ** 2.0).sum(-1))[..., np.newaxis]
        for vector in normed:
            Q_ab[n, 0, 0] += 3.0 * vector[0] * vector[0] - 1
            Q_ab[n, 0, 1] += 3.0 * vector[0] * vector[1]
            Q_ab[n, 0, 2] += 3.0 * vector[0] * vector[2]
            Q_ab[n, 1, 0] += 3.0 * vector[1] * vector[0]
            Q_ab[n, 1, 1] += 3.0 * vector[1] * vector[1] - 1
            Q_ab[n, 1, 2] += 3.0 * vector[1] * vector[2]
            Q_ab[n, 2, 0] += 3.0 * vector[2] * vector[0]
            Q_ab[n, 2, 1] += 3.0 * vector[2] * vector[1]
            Q_ab[n, 2, 2] += 3.0 * vector[2] * vector[2] - 1
    Q_ab /= (2.0 * normed.shape[0])
    return Q_ab


def _compute_nematic_director(traj):
    """Compute the characteristic vector describing a trajectory's orientation.

    In this definition, the long molecular axis is found from the inertia
    tensor, I_{ab], and is taken to be the eigenvector associated with the
    smallest eigenvalue of I_{ab}.

    See [1] for brief summary and discussion on other methods to obtain the
    director.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute orientation of.

    Returns
    -------
    directors :  np.ndarray, shape=(traj.n_frames, 3), dtype=float64
        Characteristic vectors describing the trajectory for each frame.

    See also
    --------
    compute_inertia_tensor

    References
    ----------
    .. [1] http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html

    """
    inertia_tensor = compute_inertia_tensor(traj)
    # TODO: Is there a cleaner way to do this broadcasting? Closer to this which
    # does not work:    v[:, :, np.argmin(w, axis=1)]
    w, v = np.linalg.eig(inertia_tensor)
    return np.array([v[:, :, x][i] for i, x in enumerate(np.argmin(w, axis=1))])


def compute_inertia_tensor(traj):
    """Compute the inertia tensor of a trajectory.

    For each frame,
        I_{ab} = sum_{i_atoms} [m_i * (r_i^2 * d_{ab} - r_{ia} * r_{ib})]

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute inertia tensor of.

    Returns
    -------
    I_ab:  np.ndarray, shape=(traj.n_frames, 3, 3), dtype=float64
        Inertia tensors for each frame.

    """
    I_ab = np.zeros(shape=(traj.n_frames, 3, 3), dtype=np.float64)
    com = compute_center_of_mass(traj)
    masses = np.array([atom.element.mass for atom in traj.top.atoms])

    # TODO: Look into vectorizing.
    for n, xyz in enumerate(traj.xyz):
        for i, coord in enumerate(xyz):
            mass = masses[i]
            coord = coord - com[n]
            I_ab[n, 0, 0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
            I_ab[n, 1, 1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
            I_ab[n, 2, 2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
            I_ab[n, 0, 1] -= mass * coord[0] * coord[1]
            I_ab[n, 0, 2] -= mass * coord[0] * coord[2]
            I_ab[n, 1, 2] -= mass * coord[1] * coord[2]
        I_ab[n, 1, 0] = I_ab[n, 0, 1]
        I_ab[n, 2, 0] = I_ab[n, 0, 2]
        I_ab[n, 2, 1] = I_ab[n, 1, 2]
    return I_ab
