##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors: Tim Moore
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

from mdtraj.geometry.distance import compute_center_of_mass
from mdtraj.utils import ensure_type

__all__ = ["compute_nematic_order", "compute_inertia_tensor", "compute_directors"]


def compute_nematic_order(traj, indices="chains"):
    """Compute the nematic order parameter of a group in every frame.

    The nematic order parameter describes the orientational order of a system
    with a value between 0 and 1. A completely isotropic system, such as a
    liquid, has no preferred direction and a nematic order parameter of 0. An
    anisotropic system, such as many liquid crystals, monolayers or bilayers,
    have a preferred orientation and will have a positive order parameter where
    a value of 1 signifies perfect ordering.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute ordering in.
    indices : {'chains', 'residues', list of lists}, optional, default='chains'
        The group to consider. Users can pass their own indices as a list of
        lists with the "shape" (n_compounds, len(each_compound)).
        Recognized string keywords are 'chains' and 'residues'.

    Returns
    -------
    S2 : np.ndarray, shape=(traj.n_frames,), dtype=float64
        Nematic order parameter values in every frame.

    See also
    --------
    compute_directors

    References
    ----------
    .. [1] Allen, M. P.; Tildesley , D. J. (1987), "Computer Simulation of
           Liquids", Ch. 11.5
    .. [2] http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
    .. [3] http://cmt.dur.ac.uk/sjc/thesis_dlc/node19.html

    Examples
    --------
    Ordering of chains in an alkylsilane monolayer of C10H31-Si(OH)2-:

    >>> import mdtraj as md
    >>> from mdtraj.testing import get_fn
    >>> traj = md.load(get_fn('monolayer.xtc'), top=get_fn('monolayer.pdb'))
    >>> # Each of the 100 chains contains 36 atoms.
    >>> chain_indices = [[n+x for x in range(36)] for n in range(0, 3600, 36)]
    >>> S2 = md.compute_nematic_order(traj, indices=chain_indices)

    The chains were attached to a flat, crystalline substrate and are thus
    highly ordered with a mean S2 of ~0.996.

    >>> traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    >>> water_indices = [[n+x for x in range(3)] for n in range(0, 3600, 3)]
    >>> S2 = md.compute_nematic_order(traj, indices=water_indices)

    This water box is essentially isotropic and has a mean S2 of ~0.042.

    """
    # Compute the directors for each compound for each frame.
    all_directors = compute_directors(traj, indices)

    # From the directors, compute the Q-tensor and nematic order parameter, S2.
    Q_ab = _compute_Q_tensor(all_directors)

    w = np.linalg.eigvals(Q_ab)
    S2 = w.max(axis=1)

    return S2


def compute_directors(traj, indices="chains"):
    """Compute the characteristic vector describing the orientation of each group

    In this definition, the long molecular axis is found from the inertia
    tensor, :math:`I_{ab}`, and is taken to be the eigenvector associated with the
    smallest eigenvalue of :math:`I_{ab}`.

    See [1] for brief summary and discussion on other methods to obtain the
    director.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute orientation of.
    indices : {'chains', 'residues', list of lists}, optional, default='chains'
        The group to consider. Users can pass their own indices as a list of
        lists with the "shape" (n_compounds, len(each_compound)).
        Recognized string keywords are 'chains' and 'residues'.

    Returns
    -------
    directors : np.ndarray, shape=(traj.n_frames, len(indices), 3), dtype=float64
        Characteristic vectors describing the trajectory for each frame.

    See also
    --------
    compute_nematic_order
    compute_inertia_tensor

    Notes
    -----
    Since there is no preferred orientation of the director, the director
    :math:`n` has the same meaning as :math:`-n`.
    Therefore, care should be taken to ensure the director is pointing in
    the direction you think it is, e.g., by contraining it to a hemisphere that
    makes physical sense.


    References
    ----------
    .. [1] http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html

    """
    indices = _get_indices(traj, indices)
    all_directors = np.empty(
        shape=(traj.n_frames, len(indices), 3),
        dtype=np.float64,
    )
    for i, ids in enumerate(indices):
        sub_traj = traj.atom_slice(ids)
        director = _compute_director(sub_traj)
        all_directors[:, i, :] = director
    return all_directors


def compute_inertia_tensor(traj):
    """Compute the inertia tensor of a trajectory.

    For each frame,

    .. math::

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
    center_of_mass = np.expand_dims(compute_center_of_mass(traj), axis=1)
    xyz = traj.xyz - center_of_mass
    masses = np.array([atom.element.mass for atom in traj.top.atoms])

    eyes = np.empty(shape=(traj.n_frames, 3, 3), dtype=np.float64)
    eyes[:] = np.eye(3)
    A = np.einsum("i, kij->k", masses, xyz**2).reshape(traj.n_frames, 1, 1)
    B = np.einsum("ij..., ...jk->...ki", masses[:, np.newaxis] * xyz.T, xyz)
    return A * eyes - B


def _get_indices(traj, indices):
    if isinstance(indices, str):
        if indices.lower() == "chains":
            group = list(traj.top.chains)
        elif indices.lower() == "residues":
            group = list(traj.top.residues)
        else:
            raise ValueError(f"Invalid selection: {indices}")
        indices = [[at.index for at in compound.atoms] for compound in group]
    else:
        # TODO: Clean way to ensure that indices is a list of lists of ints?
        # This may be easier than I'm thinking but `ensure_type` won't work
        # since sub-lists of variable lengths should be allowed.
        # E.g. [[1, 2], [3, 4, 5], [6, 7]] should be valid.
        if isinstance(indices, (list, tuple)):
            for sublist in indices:
                if not isinstance(sublist, (list, tuple)):
                    raise ValueError(f"Invalid selection: {indices}")
                for index in sublist:
                    if not isinstance(index, int):
                        raise ValueError(
                            f"Indices must be integers: {indices}",
                        )
        else:
            raise ValueError(f"Invalid selection: {indices}")
    return indices


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
    _compute_director

    References
    ----------
    .. [1] Allen, M. P.; Tildesley , D. J. (1987), "Computer Simulation of
           Liquids", p. 305, Eq. 11.19

    """

    all_directors = ensure_type(
        all_directors,
        dtype=np.float64,
        ndim=3,
        name="directors",
        shape=(None, None, 3),
    )

    normed = all_directors / np.linalg.norm(all_directors, axis=2)[..., np.newaxis]

    Q_ab = np.zeros(shape=(all_directors.shape[0], 3, 3), dtype=np.float64)

    for n, normed_vectors in enumerate(normed):
        for vector in normed_vectors:
            Q_ab[n, 0, 0] += 3.0 * vector[0] * vector[0] - 1
            Q_ab[n, 0, 1] += 3.0 * vector[0] * vector[1]
            Q_ab[n, 0, 2] += 3.0 * vector[0] * vector[2]
            Q_ab[n, 1, 0] += 3.0 * vector[1] * vector[0]
            Q_ab[n, 1, 1] += 3.0 * vector[1] * vector[1] - 1
            Q_ab[n, 1, 2] += 3.0 * vector[1] * vector[2]
            Q_ab[n, 2, 0] += 3.0 * vector[2] * vector[0]
            Q_ab[n, 2, 1] += 3.0 * vector[2] * vector[1]
            Q_ab[n, 2, 2] += 3.0 * vector[2] * vector[2] - 1

    Q_ab /= 2.0 * all_directors.shape[1]

    return Q_ab


def _compute_director(traj):
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

    # Only works with numpy >= 1.8.
    # TODO: Is there a cleaner way to do this broadcasting? Closer to this which
    # does not work:    v[:, :, np.argmin(w, axis=1)]
    w, v = np.linalg.eig(inertia_tensor)
    directors = np.array([v[:, :, x][i] for i, x in enumerate(np.argmin(w, axis=1))])

    return directors


# Pure python reference implementations for testing #
def _compute_inertia_tensor_slow(traj):
    """Compute the inertia tensor of a trajectory."""
    center_of_mass = np.expand_dims(compute_center_of_mass(traj), axis=1)
    centered_xyz = traj.xyz - center_of_mass
    masses = np.array([atom.element.mass for atom in traj.top.atoms])

    I_ab = np.zeros(shape=(traj.n_frames, 3, 3), dtype=np.float64)
    for n, xyz in enumerate(centered_xyz):
        for i, r in enumerate(xyz):
            mass = masses[i]
            I_ab[n, 0, 0] += mass * (r[1] * r[1] + r[2] * r[2])
            I_ab[n, 1, 1] += mass * (r[0] * r[0] + r[2] * r[2])
            I_ab[n, 2, 2] += mass * (r[0] * r[0] + r[1] * r[1])
            I_ab[n, 0, 1] -= mass * r[0] * r[1]
            I_ab[n, 0, 2] -= mass * r[0] * r[2]
            I_ab[n, 1, 2] -= mass * r[1] * r[2]
        I_ab[n, 1, 0] = I_ab[n, 0, 1]
        I_ab[n, 2, 0] = I_ab[n, 0, 2]
        I_ab[n, 2, 1] = I_ab[n, 1, 2]
    return I_ab
