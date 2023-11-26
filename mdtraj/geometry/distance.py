##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp, Jason Swails
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
from mdtraj.utils import ensure_type
from mdtraj.utils.six.moves import range
from mdtraj.utils.unitcell import box_vectors_to_lengths_and_angles
from . import _geometry

__all__ = ['compute_distances_core', 'compute_distances',
           'compute_distances_t', 'compute_displacements',
           'compute_center_of_mass', 'compute_center_of_geometry',
           'find_closest_contact']

def compute_distances_core(
        positions,
        atom_pairs,
        unitcell_vectors = None,
        periodic=True,
        opt = True,
):

    """Compute the distances between pairs of atoms in each frame.

    Parameters
    ----------
    positions : np.ndarray of shape=(n_frames, n_atoms, 3), dtype=float
        The positions of all atoms for a given trajectory.

    atom_pairs : np.ndarray of shape=(num_pairs, 2), dtype=int
        Each row gives the indices of two atoms involved in the interaction.

    unitcell_vectors : None or np.ndarray of shape(n_frames, 3 x 3), default=None
        A numpy array that specifies the box vectors for all frames for a trajectory.

    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    opt : bool, default=True
        Use an optimized native library to calculate distances. Our optimized
        SSE minimum image convention calculation implementation is over 1000x
        faster than the naive numpy implementation.

    Returns
    -------

    distances : np.ndarray, shape=(n_frames, num_pairs), dtype=float
        The distance, in each frame, between each pair of atoms.
    """

    xyz = ensure_type(positions, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    pairs = ensure_type(atom_pairs, dtype=np.int32, ndim=2, name='atom_pairs', shape=(None, 2), warn_on_cast=False)
    if not np.all(np.logical_and(pairs < positions.shape[1], pairs >= 0)):
        raise ValueError('atom_pairs must be between 0 and %d' % traj.n_atoms)

    if len(pairs) == 0:
        return np.zeros((len(xyz), 0), dtype=np.float32)

    if periodic and (unitcell_vectors is not None):

        box = ensure_type(
            unitcell_vectors,
            dtype=np.float32,
            ndim=3,
            name='unitcell_vectors',
            shape=(len(xyz), 3, 3),
            warn_on_cast=False,
        )

        # convert to angles
        unitcell_angles = []
        for fr_unitcell_vectors in unitcell_vectors:
            _, _, _, alpha, beta, gamma = box_vectors_to_lengths_and_angles(
                fr_unitcell_vectors[0],
                fr_unitcell_vectors[1],
                fr_unitcell_vectors[2],
            )
            unitcell_angles.append(np.array([alpha, beta, gamma]))

        orthogonal = np.allclose(np.array(unitcell_angles), 90)

        if opt:

            out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)
            _geometry._dist_mic(xyz, pairs, box.transpose(0, 2, 1).copy(), out, orthogonal)
            return out

        else:

            return _distance_mic(xyz, pairs, box.transpose(0, 2, 1), orthogonal)

    # either there are no unitcell vectors or they dont want to use them
    if opt:
        out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)
        _geometry._dist(xyz, pairs, out)
        return out
    else:
        return _distance(xyz, pairs)


def compute_distances(traj, atom_pairs, periodic=True, opt=True):
    """Compute the distances between pairs of atoms in each frame.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    atom_pairs : np.ndarray, shape=(num_pairs, 2), dtype=int
        Each row gives the indices of two atoms involved in the interaction.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Our optimized
        SSE minimum image convention calculation implementation is over 1000x
        faster than the naive numpy implementation.

    Returns
    -------
    distances : np.ndarray, shape=(n_frames, num_pairs), dtype=float
        The distance, in each frame, between each pair of atoms.
    """

    return compute_distances_core(
        traj.xyz,
        atom_pairs,
        unitcell_vectors=traj.unitcell_vectors,
        periodic=periodic,
        opt=opt,
    )

def compute_distances_t(traj, atom_pairs, time_pairs, periodic=True, opt=True):
    """Compute the distances between pairs of atoms at pairs of times.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    atom_pairs : np.ndarray, shape=(num_atom_pairs, 2), dtype=int
        Each row gives the indices of two atoms involved in the interaction.
    time_pairs : nd.array, shape=(num_times, 2), dtype=int
        Each row gives the indices of two frames.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Our optimized
        SSE minimum image convention calculation implementation is over 1000x
        faster than the naive numpy implementation.

    Returns
    -------
    distances : np.ndarray, shape=(num_times, num_atom_pairs), dtype=float
        The distance between each pair of atoms at t=t1 and t=t2.
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    pairs = ensure_type(atom_pairs, dtype=np.int32, ndim=2, name='atom_pairs', shape=(None, 2), warn_on_cast=False)
    times = ensure_type(time_pairs, dtype=np.int32, ndim=2, name='time_pairs', shape=(None, 2), warn_on_cast=False)

    if not np.all(np.logical_and(pairs < traj.n_atoms, pairs >= 0)):
        raise ValueError('atom_pairs must be between 0 and %d' % traj.n_atoms)

    if not np.all(np.logical_and(times < traj.n_frames, times >= 0)):
        raise ValueError('time_pairs must be between 0 and %d' % traj.n_frames)

    if len(pairs) == 0:
        return np.zeros((len(xyz), 0), dtype=np.float32)

    if periodic and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3),
                          warn_on_cast=False)
        orthogonal = np.allclose(traj.unitcell_angles, 90)
        if opt:
            out = np.empty((times.shape[0], pairs.shape[0]), dtype=np.float32)
            _geometry._dist_mic_t(xyz, pairs, times, box.transpose(0, 2, 1).copy(), out, orthogonal)
            out = out.reshape((times.shape[0], pairs.shape[0]))
            return out
        else:
            return _distance_mic_t(xyz, pairs, times, box.transpose(0, 2, 1), orthogonal)

    # either there are no unitcell vectors or they dont want to use them
    if opt:
        out = np.empty((times.shape[0], pairs.shape[0]), dtype=np.float32)
        _geometry._dist_t(xyz, pairs, times, out)
        return out
    else:
        return _distance_t(xyz, pairs, times)


def compute_displacements(traj, atom_pairs, periodic=True, opt=True):
    """Compute the displacement vector between pairs of atoms in each frame of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute distances in
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Our
        optimized minimum image convention calculation implementation is
        over 1000x faster than the naive numpy implementation.

    Returns
    -------
    displacements : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32
         The displacememt vector, in each frame, between each pair of atoms.
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    pairs = ensure_type(np.asarray(atom_pairs), dtype=np.int32, ndim=2, name='atom_pairs', shape=(None, 2), warn_on_cast=False)
    if not np.all(np.logical_and(pairs < traj.n_atoms, pairs >= 0)):
        raise ValueError('atom_pairs must be between 0 and %d' % traj.n_atoms)
    if len(pairs) == 0:  # If pairs is an empty slice of an array
        return np.zeros((len(xyz), 0, 3), dtype=np.float32)

    if periodic and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3),
                          warn_on_cast=False)
        orthogonal = np.allclose(traj.unitcell_angles, 90)
        if opt:
            out = np.empty((xyz.shape[0], pairs.shape[0], 3), dtype=np.float32)
            _geometry._dist_mic_displacement(xyz, pairs, box.transpose(0, 2, 1).copy(), out, orthogonal)
            return out
        else:
            return _displacement_mic(xyz, pairs, box.transpose(0, 2, 1), orthogonal)

    # either there are no unitcell vectors or they dont want to use them
    if opt:
        out = np.empty((xyz.shape[0], pairs.shape[0], 3), dtype=np.float32)
        _geometry._dist_displacement(xyz, pairs, out)
        return out
    return _displacement(xyz, pairs)


def compute_center_of_mass(traj, select=None):
    """Compute the center of mass for each frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of mass for
    select : str, optional, default=all
        a mdtraj.Topology selection string that
        defines the set of atoms of which to calculate
        the center of mass, the default is all atoms

    Returns
    -------
    com : np.ndarray, shape=(n_frames, 3)
         Coordinates of the center of mass for each frame
    """
    com = np.empty((traj.n_frames, 3))

    if select is None:
        masses = np.array([a.element.mass for a in traj.top.atoms])
        masses /= masses.sum()

        xyz = traj.xyz

    else:
        atoms_of_interest = traj.topology.select(select)

        masses = np.array([traj.top.atom(i).element.mass for i in atoms_of_interest])
        masses /= masses.sum()

        xyz = traj.xyz[:, atoms_of_interest]

    for i, x in enumerate(xyz):
        com[i, :] = x.astype('float64').T.dot(masses)
    return com


def compute_center_of_geometry(traj):
    """Compute the center of geometry for each frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of geometry for

    Returns
    -------
    com : np.ndarray, shape=(n_frames, 3)
         Coordinates of the center of geometry for each frame

    """

    centers = np.zeros((traj.n_frames, 3))

    for i, x in enumerate(traj.xyz):
        centers[i, :] = x.astype('float64').T.mean(axis=1)
    return centers


def find_closest_contact(traj, group1, group2, frame=0, periodic=True):
    """Find the closest contact between two groups of atoms.

    Given a frame of a Trajectory and two groups of atoms, identify the pair of
    atoms (one from each group) that form the closest contact between the two groups.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    group1 : np.ndarray, shape=(num_atoms), dtype=int
        The indices of atoms in the first group.
    group2 : np.ndarray, shape=(num_atoms), dtype=int
        The indices of atoms in the second group.
    frame : int, default=0
        The frame of the Trajectory to take positions from
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    Returns
    -------
    result : tuple (int, int, float)
         The indices of the two atoms forming the closest contact, and the distance between them.
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)[frame]
    atoms1 = ensure_type(group1, dtype=np.int32, ndim=1, name='group1', warn_on_cast=False)
    atoms2 = ensure_type(group2, dtype=np.int32, ndim=1, name='group2', warn_on_cast=False)
    if periodic and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(traj.xyz), 3, 3),
                          warn_on_cast=False)[frame]
    else:
        box = None
    return _geometry._find_closest_contact(xyz, atoms1, atoms2, box)


##############################################################################
# pure python implementation of the core routines
##############################################################################


def _distance(xyz, pairs):
    "Distance between pairs of points in each frame"
    delta = np.diff(xyz[:, pairs], axis=2)[:, :, 0]
    return (delta ** 2.).sum(-1) ** 0.5


def _distance_t(xyz, pairs, times):
    "Distance between pairs of points in specified frames"
    frame1 = xyz[:, pairs[:,0]][times[:,0]]
    frame2 = xyz[:, pairs[:,1]][times[:,1]]
    out = np.linalg.norm(frame1-frame2, axis=2)

    return out


def _displacement(xyz, pairs):
    "Displacement vector between pairs of points in each frame"
    value = np.diff(xyz[:, pairs], axis=2)[:, :, 0]
    assert value.shape == (xyz.shape[0], pairs.shape[0], 3), 'v.shape %s, xyz.shape %s, pairs.shape %s' % (str(value.shape), str(xyz.shape), str(pairs.shape))
    return value


def _reduce_box_vectors(vectors):
    """Make sure box vectors are in reduced form."""
    (bv1, bv2, bv3) = vectors
    bv3 -= bv2*round(bv3[1]/bv2[1]);
    bv3 -= bv1*round(bv3[0]/bv1[0]);
    bv2 -= bv1*round(bv2[0]/bv1[0]);
    return (bv1, bv2, bv3)


def _distance_mic(xyz, pairs, box_vectors, orthogonal):
    """Distance between pairs of points in each frame under the minimum image
    convention for periodic boundary conditions.

    The computation follows scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

    This is a slow pure python implementation, mostly for testing.
    """
    out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)
    for i in range(len(xyz)):
        bv1, bv2, bv3 = _reduce_box_vectors(box_vectors[i].T)

        for j, (a,b) in enumerate(pairs):
            r12 = xyz[i,b,:] - xyz[i,a,:]
            r12 -= bv3*round(r12[2]/bv3[2]);
            r12 -= bv2*round(r12[1]/bv2[1]);
            r12 -= bv1*round(r12[0]/bv1[0]);
            dist = np.linalg.norm(r12)
            if not orthogonal:
                for ii in range(-1, 2):
                    v1 = bv1*ii
                    for jj in range(-1, 2):
                        v12 = bv2*jj + v1
                        for kk in range(-1, 2):
                            new_r12 = r12 + v12 + bv3*kk
                            dist = min(dist, np.linalg.norm(new_r12))
            out[i, j] = dist
    return out


def _distance_mic_t(xyz, pairs, times, box_vectors, orthogonal):
    """Distance between pairs of points between specified frames under the minimum image
    convention for periodic boundary conditions.

    The computation is modified from scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

    This is a slow pure python implementation, mostly for testing.
    """
    out = np.empty((times.shape[0], pairs.shape[0]), dtype=np.float32)
    for i, (a,b) in enumerate(times):
        bv1, bv2, bv3 = _reduce_box_vectors(box_vectors[a].T)
        for j, (c,d) in enumerate(pairs):
            r12 = xyz[a,c] - xyz[b,d]
            r12 -= bv3*round(r12[2]/bv3[2]);
            r12 -= bv2*round(r12[1]/bv2[1]);
            r12 -= bv1*round(r12[0]/bv1[0]);
            dist = np.linalg.norm(r12)
            if not orthogonal:
                for ii in range(-1, 2):
                    v1 = bv1*ii
                    for jj in range(-1, 2):
                        v12 = bv2*jj + v1
                        for kk in range(-1, 2):
                            new_r12 = r12 + v12 + bv3*kk
                            dist = min(dist, np.linalg.norm(new_r12))
            out[i, j] = dist

    return out


def _displacement_mic(xyz, pairs, box_vectors, orthogonal):
    """Displacement vector between pairs of points in each frame under the
    minimum image convention for periodic boundary conditions.

    The computation follows scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

    This is a very slow pure python implementation, mostly for testing.
    """
    out = np.empty((xyz.shape[0], pairs.shape[0], 3), dtype=np.float32)
    for i in range(len(xyz)):
        bv1, bv2, bv3 = _reduce_box_vectors(box_vectors[i].T)
        hinv = np.linalg.inv(np.array([bv1, bv2, bv3]).T)

        for j, (a,b) in enumerate(pairs):
            r12 = xyz[i,b,:] - xyz[i,a,:]
            r12 -= bv3*round(r12[2]/bv3[2]);
            r12 -= bv2*round(r12[1]/bv2[1]);
            r12 -= bv1*round(r12[0]/bv1[0]);
            min_disp = r12
            dist2 = (r12*r12).sum()
            if not orthogonal:
                for ii in range(-1, 2):
                    v1 = bv1*ii
                    for jj in range(-1, 2):
                        v12 = bv2*jj+v1
                        for kk in range(-1, 2):
                            tmp = r12 + v12 + bv3*kk
                            new_dist2 = (tmp*tmp).sum()
                            if new_dist2 < dist2:
                                dist2 = new_dist2
                                min_disp = tmp
            out[i, j] = min_disp

    return out
