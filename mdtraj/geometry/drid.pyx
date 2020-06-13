##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2014- Stanford University and the Authors
#
# Authors: Robert McGibbon
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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division, absolute_import
import numpy as np
from mdtraj.utils import ensure_type

cimport cython
from cython.parallel import prange

__all__ = ['compute_drid']

cdef extern from "dridkernels.h":
    ctypedef signed int int32_t
    void drid_moments(float* coords, int32_t index, int32_t * partners,
                      int32_t n_partners, double* moments) nogil


##############################################################################
# Code
##############################################################################


def compute_drid(traj, atom_indices=None):
    """compute_drid(traj, atom_indices=None)

    Distribution of reciprocal interatomic distances (DRID) representation
    of an MD trajectory

    Parameters
    ----------
    traj : md.Trajectory
        A trajectory
    atom_indices : np.ndarray, dtype=int, shape=(n_atom_indices,) or None
        The indices (zero-based) of the atoms to use in the DRID representation.
        If None, all atoms will be used.

    Returns
    -------
    X : np.ndarray, shape=(n_frames, n_atom_indices*3), dtype=np.double
        A rotationally invariant vector representation of each frame in the
        simulation. The DRID vector contains the mean, second, and third
        central moments of the reciprocal interatomic distances (with 1-2
        interactions excluded).

    References
    ----------
    .. [1] Zhou, Caflisch; Distribution of Reciprocal of Interatomic Distances: A Fast Structural Metric. JCTC 2012 doi:10.1021/ct3003145
    """
    if traj.topology is None:
        raise ValueError('traj must have a topology defined')
    if atom_indices is None:
        atom_indices = np.arange(traj.n_atoms, dtype=np.int32)
    else:
        atom_indices = ensure_type(np.asarray(atom_indices), dtype=np.int32, ndim=1, name='atom_indices', warn_on_cast=False)
        if not np.all((atom_indices >= 0) * (atom_indices < traj.xyz.shape[1]) * (atom_indices < traj.xyz.shape[1])):
            raise ValueError("atom_indices must be valid positive indices")

    # atom_indices[i] = j implies that inverse_atom_indices[j] = i
    inverse_atom_indices = dict((j,i) for i, j in enumerate(atom_indices))

    # DEBUG
    # print(atom_indices)
    # print(inverse_atom_indices)

    bonds = [set() for _ in atom_indices]
    for a, b in traj.topology.bonds:
        if a.index in inverse_atom_indices and b.index in inverse_atom_indices:
            bonds[inverse_atom_indices[a.index]].add(b.index)
            bonds[inverse_atom_indices[b.index]].add(a.index)

    partners_l = []
    set_atom_indices = set(atom_indices)
    for i, j in enumerate(atom_indices):
        partners_l.append(set_atom_indices - bonds[i] - set([j]))
    n_partners = np.array([len(p) for p in partners_l], dtype=np.int32)

    partners = np.zeros((len(atom_indices), np.max(n_partners)), dtype=np.int32)
    partners.fill(-1)
    for i in range(len(atom_indices)):
        partners[i,:n_partners[i]] = np.array(sorted(partners_l[i]))

    return _drid(traj.xyz, atom_indices, partners, n_partners)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _drid(float[:, :, ::1] xyz, int32_t[::1] atom_indices,
           int32_t[:, ::1] partners, int32_t[::1] n_partners):

    cdef int i, j, n_frames, n_atoms, n_dims, n_atom_indices
    cdef double[:, :, ::1] result

    n_frames = xyz.shape[0]
    n_atoms = xyz.shape[1]
    n_dims = xyz.shape[2]
    n_atom_indices = len(atom_indices)
    result = np.zeros((n_frames, n_atom_indices, 3))
    assert n_dims == 3

    for i in range(n_frames):
        for j in prange(n_atom_indices, nogil=True):
            drid_moments(&xyz[i, 0, 0], atom_indices[j], &partners[j,0],
                         n_partners[j], &result[i, j, 0])

    return np.array(result, copy=False).reshape(n_frames, n_atom_indices*3)
