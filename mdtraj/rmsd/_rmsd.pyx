##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon, Kyle Beauchamp
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

import cython
import warnings
import numpy as np
from mdtraj.utils import ensure_type

from cpython cimport bool
from cython.parallel cimport prange

##############################################################################
# External Declarations
##############################################################################

cdef extern float msd_axis_major(int nrealatoms, int npaddedatoms, int rowstride,
    float* aT, float* bT, float G_a, float G_b) nogil
cdef extern float msd_atom_major(int nrealatoms, int npaddedatoms,  float* a,
    float* b, float G_a, float G_b, int computeRot, float rot[9]) nogil
cdef extern float rot_msd_atom_major(const int n_real_atoms,
    const int n_padded_atoms, const float* a, const float* b, const float rot[9]) nogil
cdef extern float rot_atom_major(const int n_atoms, float* a, const float rot[9]) nogil
cdef extern from "center.h":
    void inplace_center_and_trace_atom_major(float* coords, float* traces,
    const int n_frames, const int n_atoms) nogil
cdef extern from "math.h":
    float sqrtf(float x) nogil


##############################################################################
# External (Public) Functions
##############################################################################


@cython.boundscheck(False)
def rmsd(target, reference, int frame=0, atom_indices=None,
         ref_atom_indices=None, bool parallel=True, bool precentered=False):
    """rmsd(target, reference, frame=0, atom_indices=None, parallel=True, precentered=False)

    Compute RMSD of all conformations in target to a reference conformation.
    Note, this will center the conformations in place.

    Parameters
    ----------
    target : md.Trajectory
        For each conformation in this trajectory, compute the RMSD to
        a particular 'reference' conformation in another trajectory
        object.
    reference : md.Trajectory
        The object containing the reference conformation to measure distances
        to.
    frame : int, default=0
        The index of the conformation in `reference` to measure
        distances to.
    atom_indices : array_like, or None
        The indices of the atoms to use in the RMSD calculation. If not
        supplied, all atoms will be used.
    ref_atom_indices : array_like, or None
        Use these indices for the reference trajectory. If not supplied,
        the atom indices will be the same as those for target.
    parallel : bool
        Use OpenMP to calculate each of the RMSDs in parallel over
        multiple cores.
    precentered : bool, default=False
        Assume that the conformations are already centered at the origin, and that
        the "rmsd_traces" have been computed, as is done by
        `Trajectory.center_coordinates`. The "rmsd_traces" are intermediate
        calculations needed for the RMSD calculation which can be computed
        independently on each trajectory. Note that this has the potential to
        be unsafe; if you use Trajectory.center_coordinates and then modify
        the trajectory's coordinates, the center and traces will be out of
        date and the RMSDs will be incorrect.

    Examples
    --------
    >>> import mdtraj as md                                      # doctest: +SKIP
    >>> rmsds = md.rmsd(trajectory, trajectory, 0)               # doctest: +SKIP
    >>> print rmsds                                              # doctest: +SKIP
    array([ 0.0,  0.03076187,  0.02549562, ...,  0.06230228,
        0.00666826,  0.24364147])

    The calculation is slightly faster if you precenter the trajectory

    >>> trajectory.center_coordinates()
    >>> rmsds = md.rmsd(trajectory, trajectory, 0, precentered=True)

    See Also
    --------
    Trajectory.center_coordinates

    Notes
    -----
    This function uses OpenMP to parallelize the calculation across
    multiple cores. To control the number of threads launched by OpenMP,
    you can set the environment variable ``OMP_NUM_THREADS``.

    Returns
    -------
    rmsds : np.ndarray, shape=(target.n_frames,)
        A 1-D numpy array of the optimal root-mean-square deviations from
        the `frame`-th conformation in reference to each of the conformations
        in target.
    """
    # import time
    cdef bool atom_indices_is_none = False

    if atom_indices is None:
        atom_indices_is_none = True
        atom_indices = slice(None)
    else:
        atom_indices = ensure_type(np.asarray(atom_indices), dtype=int, ndim=1, name='atom_indices')
        if not np.all((atom_indices >= 0) *
                              (atom_indices < target.xyz.shape[1])):
            raise ValueError("atom_indices must be valid positive indices")

    if ref_atom_indices is None:
        ref_atom_indices = atom_indices
    else:
        if len(ref_atom_indices) != len(atom_indices):
            raise ValueError("atom_indices and ref_atom_indices must have same number of atom indices. "
                             "found %d and %d." % (len(atom_indices), len(ref_atom_indices)))

    if not isinstance(ref_atom_indices, slice):
        ref_atom_indices = ensure_type(np.asarray(ref_atom_indices), dtype=int, ndim=1, name='ref_atom_indices')
        if not np.all((ref_atom_indices >= 0) *
                              (ref_atom_indices < reference.xyz.shape[1])):
            raise ValueError("ref_atom_indices must be valid positive indices")

    # Error checks
    assert (target.xyz.ndim == 3) and (reference.xyz.ndim == 3) and (target.xyz.shape[2]) == 3 and (reference.xyz.shape[2] == 3)
    if not ((target.xyz.shape[1]  == reference.xyz.shape[1]) or
            (len(ref_atom_indices) == len(atom_indices))):
        raise ValueError("Input trajectories must have same number of atoms. "
                         "found %d and %d." % (target.xyz.shape[1], reference.xyz.shape[1]))
    if frame >= reference.xyz.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: reference has "
                         "only %d frames." % (frame, reference.xyz.shape[0]))

    # static declarations
    cdef int i
    cdef float msd, ref_g
    cdef float[:, :, :] target_xyz
    cdef float[:, :] ref_xyz_frame
    cdef float[:] target_g
    cdef int target_n_frames = target.xyz.shape[0]
    cdef int n_atoms = target.xyz.shape[1] if np.all(atom_indices == slice(None)) else len(atom_indices)

    # make sure *every* frame in target_xyz is in proper c-major order
    target_xyz = np.asarray(target.xyz[:, atom_indices, :], order='C', dtype=np.float32)
    # only extract the `frame`-th conformation from ref_xyz
    ref_xyz_frame = np.asarray(reference.xyz[frame, ref_atom_indices, :], order='C', dtype=np.float32)

    # t0 = time.time()
    if precentered and (reference._rmsd_traces is not None) and (target._rmsd_traces is not None) and atom_indices_is_none:
        target_g = np.asarray(target._rmsd_traces, order='C', dtype=np.float32)
        ref_g = reference._rmsd_traces[frame]
    else:
        if precentered:
            warnings.warn(
                'in rmsd(), precentered is ignored when atom_indices != None',
                RuntimeWarning)
        target_g = np.empty(target_n_frames, dtype=np.float32)
        inplace_center_and_trace_atom_major(&target_xyz[0,0,0], &target_g[0], target_n_frames, n_atoms)
        inplace_center_and_trace_atom_major(&ref_xyz_frame[0, 0], &ref_g, 1, n_atoms)

    # t1 = time.time()

    cdef float[:] distances = np.zeros(target_n_frames, dtype=np.float32)
    if parallel:
        for i in prange(target_n_frames, nogil=True):
            msd = msd_atom_major(n_atoms, n_atoms, &target_xyz[i, 0, 0], &ref_xyz_frame[0, 0], target_g[i], ref_g, 0, NULL)
            distances[i] = sqrtf(msd)
    else:
        for i in range(target_n_frames):
            msd = msd_atom_major(n_atoms, n_atoms, &target_xyz[i, 0, 0], &ref_xyz_frame[0, 0], target_g[i], ref_g, 0, NULL)
            distances[i] = sqrtf(msd)

    # t2 = time.time()
    # print 'rmsd: %s, centering: %s' % (t2-t1, t1-t0)
    return np.array(distances, copy=False)


@cython.boundscheck(False)
def rmsf(target, reference, int frame=0, atom_indices=None,
         ref_atom_indices=None, bool parallel=True, bool precentered=False):
    """rmsf(target, reference, frame=0, atom_indices=None, parallel=True, precentered=False)
    
    Compute RMSF of atom positions in target trajectory. This will center target conformations in place.
    
    Parameters
    ----------
    target : md.Trajectory
        Compute the RMSF of atom positions in target trajectory. 
    reference : md.Trajectory, or None
        A trajectory with the same number of atoms as in the target to be used
        as a reference. If None, the average xyz positions of each atom will be used.
    frame : int, default=0
        The index of the conformation in `reference` to measure
        distances to.
    atom_indices : array_like, or None
        The indices of the atoms to use in the RMSF calculation. If not
        supplied, all atoms will be used. If the atoms used in alignment
        and the atoms used in RMSF calculations are different, the trajectory
        must be aligned first manually.
    ref_atom_indices : array_like, or None
        Use these indices for the reference trajectory. If not supplied,
        the atom indices will be the same as those for target.
    parallel : bool
        Use OpenMP to calculate each of the RMSFs in parallel over
        multiple cores.
    precentered : bool, default=False
        Assume that the conformations are already centered at the origin, and that
        the "rmsd_traces" have been computed, as is done by
        `Trajectory.center_coordinates`. The "rmsd_traces" are intermediate
        calculations needed for the RMSF calculation which can be computed
        independently on each trajectory. Note that this has the potential to
        be unsafe; if you use Trajectory.center_coordinates and then modify
        the trajectory's coordinates, the center and traces will be out of
        date and the RMSFs will be incorrect.
    
    Examples
    --------
    >>> import mdtraj as md                                      # doctest: +SKIP
    >>> rmsf = md.rmsf(trajectory, trajectory, 0)                # doctest: +SKIP
    >>> print rmsf                                               # doctest: +SKIP
    array([ 0.0,  0.03076187,  0.02549562, ...,  0.06230228,
        0.00666826,  0.24364147])
    
    The calculation is slightly faster if you precenter the trajectory
    >>> trajectory.center_coordinates()
    >>> rmsf = md.rmsf(trajectory, trajectory, 0, precentered=True)
    
    If the atoms used in alignment and RMSF calculations are different,
    align the trajectory before using this method.
    >>> trajectory.superpose(atom_indices=atom_indices)
    >>> rmsf = md.rmsf(trajectory, None, atom_indices=rmsf_atom_indices)
    
    See Also
    --------
    Trajectory.center_coordinates
    
    Notes
    -----
    This function uses OpenMP to parallelize the calculation across
    multiple cores. To control the number of threads launched by OpenMP,
    you can set the environment variable ``OMP_NUM_THREADS``.
    
    Returns
    -------
    rmsf : np.ndarray, shape=(atom_indices,)
        A 1-D numpy array of the optimal root-mean-square fluctuations of the
       selected atoms.
    """
    # import time
    cdef bool atom_indices_is_none = False
    cdef bool trajectory_prealigned = False

    if reference is None:
        trajectory_prealigned = True
        reference = target

    if atom_indices is None:
        atom_indices_is_none = True
        atom_indices = slice(None)
    else:
        atom_indices = ensure_type(np.asarray(atom_indices), dtype=int, ndim=1, name='atom_indices')
        if not np.all((atom_indices >= 0) *
                              (atom_indices < target.xyz.shape[1]) *
                              (atom_indices < reference.xyz.shape[1])):
            raise ValueError("atom_indices must be valid positive indices")

    if ref_atom_indices is None:
        ref_atom_indices = atom_indices
    else:
        if len(ref_atom_indices) != len(atom_indices):
            raise ValueError("atom_indices and ref_atom_indices must have same number of atom indices. "
                             "found %d and %d." % (len(atom_indices), len(ref_atom_indices)))

    if not isinstance(ref_atom_indices, slice):
        ref_atom_indices = ensure_type(np.asarray(ref_atom_indices), dtype=int, ndim=1, name='ref_atom_indices')
        if not np.all((ref_atom_indices >= 0) *
                              (ref_atom_indices < target.xyz.shape[1]) *
                              (ref_atom_indices < reference.xyz.shape[1])):
            raise ValueError("ref_atom_indices must be valid positive indices")

    # Error checks
    assert (target.xyz.ndim == 3) and (reference.xyz.ndim == 3) and (target.xyz.shape[2]) == 3 and (reference.xyz.shape[2] == 3)
    if not ((target.xyz.shape[1]  == reference.xyz.shape[1]) or
            (len(ref_atom_indices) == len(atom_indices))):
        raise ValueError("Input trajectories must have same number of atoms. "
                         "found %d and %d." % (target.xyz.shape[1], reference.xyz.shape[1]))
    if frame >= reference.xyz.shape[0]:
        raise ValueError("Cannot calculate RMSF of frame %d: reference has "
                         "only %d frames." % (frame, reference.xyz.shape[0]))

    # static declarations
    cdef int i, j
    cdef float msd, ref_g
    cdef float[:, :, :] target_xyz
    cdef float[:, :] ref_xyz_frame
    cdef float[:, :, :] target_displaced_xyz
    cdef float[:, :] avg_xyz_frame
    cdef float[:] target_g
    cdef int target_n_frames = target.xyz.shape[0]
    cdef int n_atoms = target.xyz.shape[1] if np.all(atom_indices == slice(None)) else len(atom_indices)

    # make sure *every* frame in target_xyz is in proper c-major order
    target_xyz = np.asarray(target.xyz[:, atom_indices, :], order='C', dtype=np.float32)
    # only extract the `frame`-th conformation from ref_xyz
    ref_xyz_frame = np.asarray(reference.xyz[frame, ref_atom_indices, :], order='C', dtype=np.float32)
    avg_xyz_frame = np.zeros_like(reference.xyz[frame, ref_atom_indices, :], order='C', dtype=np.float32)

    # t0 = time.time()
    if precentered and (reference._rmsd_traces is not None) and (target._rmsd_traces is not None) and atom_indices_is_none:
        target_g = np.asarray(target._rmsd_traces, order='C', dtype=np.float32)
        ref_g = reference._rmsd_traces[frame]
    else:
        if precentered:
            warnings.warn(
                'in rmsd(), precentered is ignored when atom_indices != None',
                RuntimeWarning)
        target_g = np.empty(target_n_frames, dtype=np.float32)
        inplace_center_and_trace_atom_major(&target_xyz[0,0,0], &target_g[0], target_n_frames, n_atoms)
        inplace_center_and_trace_atom_major(&ref_xyz_frame[0, 0], &ref_g, 1, n_atoms)

    # t1 = time.time()

    cdef float[:] fluctuations = np.zeros(n_atoms, dtype=np.float32)
    # we could get away with using less memory here: we only need 1 rotation
    # matrix in serial mode and 1 rotation matrix per thread in parallel mode.
    cdef float[:, :, ::1] rot = np.zeros((target_n_frames, 3, 3), dtype=np.float32)
    target_displaced_xyz = np.array(target.xyz[:, atom_indices, :], copy=True)

    if trajectory_prealigned == False:
        if parallel:
            for i in prange(target_n_frames, nogil=True):
                msd_atom_major(n_atoms, n_atoms, &target_xyz[i, 0, 0], &ref_xyz_frame[0, 0], ref_g, target_g[i], 1, &rot[i, 0, 0])
                rot_atom_major(n_atoms, &target_displaced_xyz[i, 0, 0], &rot[i, 0, 0])
        else:
            for i in range(target_n_frames):
                msd_atom_major(n_atoms, n_atoms, &target_xyz[i, 0, 0], &ref_xyz_frame[0, 0], ref_g, target_g[i], 1, &rot[i, 0, 0])
                rot_atom_major(n_atoms, &target_displaced_xyz[i, 0, 0], &rot[i, 0, 0])

    if parallel:
        for j in prange(n_atoms, nogil=True):
            for i in range(target_n_frames):
                avg_xyz_frame[j, 0] += target_displaced_xyz[i, j, 0] / target_n_frames
                avg_xyz_frame[j, 1] += target_displaced_xyz[i, j, 1] / target_n_frames
                avg_xyz_frame[j, 2] += target_displaced_xyz[i, j, 2] / target_n_frames
    else:
        for i in range(target_n_frames):
            for j in range(n_atoms):
                avg_xyz_frame[j, 0] += target_displaced_xyz[i, j, 0] / target_n_frames
                avg_xyz_frame[j, 1] += target_displaced_xyz[i, j, 1] / target_n_frames
                avg_xyz_frame[j, 2] += target_displaced_xyz[i, j, 2] / target_n_frames

    if parallel:
        for j in prange(n_atoms, nogil=True):
            for i in range(target_n_frames):
                fluctuations[j] += ((target_displaced_xyz[i, j, 0] - avg_xyz_frame[j, 0])**2 +
                                    (target_displaced_xyz[i, j, 1] - avg_xyz_frame[j, 1])**2 +
                                    (target_displaced_xyz[i, j, 2] - avg_xyz_frame[j, 2])**2) / target_n_frames
    else:
        for i in range(target_n_frames):
            for j in range(n_atoms):
                fluctuations[j] += ((target_displaced_xyz[i, j, 0] - avg_xyz_frame[j, 0])**2 +
                                    (target_displaced_xyz[i, j, 1] - avg_xyz_frame[j, 1])**2 +
                                    (target_displaced_xyz[i, j, 2] - avg_xyz_frame[j, 2])**2) / target_n_frames

    for j in range(n_atoms):
        fluctuations[j] = sqrtf(fluctuations[j])

    # t2 = time.time()
    # print 'rmsd: %s, centering: %s' % (t2-t1, t1-t0)
    return np.array(fluctuations, copy=False)


def _center_inplace_atom_major(float[:, :, ::1] xyz not None):
    assert xyz.shape[2] == 3
    cdef float[:] traces = np.empty(xyz.shape[0], dtype=np.float32)
    inplace_center_and_trace_atom_major(&xyz[0,0,0], &traces[0], xyz.shape[0], xyz.shape[1])
    return np.array(traces, copy=False)



##############################################################################
# Private Functions
##############################################################################


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_axis_major(float[:, :, ::1] xyz1 not None,
                                float[:, :, ::1] xyz2 not None, float[::1] g1,
                                float[::1] g2, int frame, bool parallel=True):
    """getMultipleRMSDs_axis_major(xyz1, xyz2, g1, g, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in axis-major orders

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, 3, n_atoms), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, 3, n_atoms), dtype=float32
        Coordinates of frames to iterate over
    g1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    frame : int
        Index of the desired reference frame in xyz1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDS between xyz1[frame] and all of xyz2
    """

    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms = xyz1.shape[2]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[1] == 3 and xyz2.shape[1] == 3
    if not (xyz1.shape[2]  == xyz2.shape[2]):
        raise ValueError("Input arrays must have same number of atoms. "
                         "found %d and %d." % (xyz1.shape[2], xyz2.shape[2]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef float[:] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_axis_major(n_atoms, n_atoms, n_atoms,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_axis_major(n_atoms, n_atoms, n_atoms,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)

    return np.array(distances, copy=False)


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_atom_major(float[:, :, ::1] xyz1 not None, float[:, :, ::1] xyz2 not None,
                                float[::1] g1 not None, float[::1] g2 not None, int frame,
                                bool parallel=True):
    """getMultipleRMSDs_atom_major(xyz1, xyz2, g1, g2, n_atoms, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in atom-major orders

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        Coordinates of frames to iterate over
    g1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    frame : int
        Index of the desired reference frame in xyz1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDS between xyz1[frame] and all of xyz2
    """

    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms = xyz1.shape[1]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[2] == 3 and xyz2.shape[2] == 3
    if not (xyz1.shape[1]  == xyz2.shape[1]):
        raise ValueError("Input arrays must have same number of atoms. "
                         "found %d and %d." % (xyz1.shape[1], xyz2.shape[1]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef float[:] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_atom_major(n_atoms, n_atoms,  &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i], 0, NULL)
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_atom_major(n_atoms, n_atoms, &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i], 0, NULL)
            distances[i] = sqrtf(msd)

    return np.array(distances, copy=False)


@cython.boundscheck(False)
@cython.wraparound(False)
def superpose_atom_major(float[:, :, ::1] xyz_align_target not None,
                         float[:, :, ::1] xyz_align_mobile not None,
                         float[::1] g_target not None, float[::1] g_mobile not None,
                         float[:, :, ::1] xyz_displace_mobile not None, int target_frame,
                         bool parallel=True):
    """superpose_atom_major(xyz_target, xyz_mobile, g_target, g_mobile, xyz_mobile_displace)

    Superpose each frame in xyz2 upon a frame in xyz1

    Parameters
    ----------
    xyz_align_target : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        Coordinates of reference frame to align to.
    xyz_align_mobile: np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        Coordinates of the mobile trajectory to align to target the target.
    g_target : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz_target
    g_mobile : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz_mobile
    xyz_displace_mobile : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        The coordinates of the mobile trajectory to displace
    target_frame : int
        The particular frame in xyz_align_target / g_target to align the mobile
        trajectory  to.
    parallel : bool, default=True
        Run the calculation using multiple cores simultaneously.
    """
    cdef int i
    cdef int n_frames = xyz_align_mobile.shape[0]
    cdef int n_atoms_align = xyz_align_target.shape[1]
    cdef int n_atoms_displace = xyz_displace_mobile.shape[1]
    if not (xyz_align_target.shape[1]  == xyz_align_mobile.shape[1]):
        raise ValueError("Input arrays must have same number of atoms. "
                         "found %d and %d." % (
                         xyz_align_target.shape[1],
                         xyz_align_mobile.shape[1]))
    if not xyz_align_mobile.shape[0] == xyz_displace_mobile.shape[0]:
        raise ValueError("xyz_align_mobile and xyz_displace_mobile must contain the same number of frames")

    # we could get away with using less memory here: we only need 1 rotation
    # matrix in serial mode and 1 rotation matrix per thread in parallel mode.
    cdef float[:, :, ::1] rot = np.zeros((n_frames, 3, 3), dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd_atom_major(n_atoms_align, n_atoms_align, &xyz_align_mobile[i, 0, 0],
                           &xyz_align_target[target_frame, 0, 0],
                           g_target[target_frame], g_mobile[i], 1, &rot[i, 0, 0])
            rot_atom_major(n_atoms_displace, &xyz_displace_mobile[i, 0, 0], &rot[i, 0, 0])
    else:
        for i in range(n_frames):
            msd_atom_major(n_atoms_align, n_atoms_align, &xyz_align_mobile[i, 0, 0],
                           &xyz_align_target[target_frame, 0, 0],
                           g_target[target_frame], g_mobile[i], 1, &rot[i, 0, 0])
            rot_atom_major(n_atoms_displace, &xyz_displace_mobile[i, 0, 0], &rot[i, 0, 0])


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleAlignDisplaceRMSDs_atom_major(float[:, :, ::1] xyz_align1 not None,
    float[:, :, ::1] xyz_align2 not None, float[::1] g_align1 not None,
    float[::1] g_align2 not None, float[:, :, ::1] xyz_displ1 not None,
    float[:, :, ::1] xyz_displ2 not None, int n_atoms_align,
    int n_atoms_displ, int frame, bool parallel=True):
    """getMultipleAlignDisplaceRMSDs_atom_major(xyz_align1, xyz_align2, g_align1, g_align2, xyz_displ1, xyz_displ2, n_atoms_align, n_atoms_displ, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in atom-major orders, where alignments are computed
    by determining the optimal rotation matrix to align `xyz_align1` against
    `xyz_align2`, and then using this rotation matrix displace `xyz_displ1`,
    and computing its mean-squared-displacement with respect to `xyz_displ2`.

    Parameters
    ----------
    xyz_align1 : np.ndarray, shape=(n_frames, n_atoms_align_padded, 3), dtype=float32
        Alignment coordinates of reference frame.
    xyz_align2 : np.ndarray, shape=(n_frames, n_atoms_align_padded, 3), dtype=float32
        Alignment coordinates of frames to iterate over.
    g_align1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g_align2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    xyz_displ1 : np.ndarray, shape=(n_frames, n_atoms_displ_padded, 3), dtype=float32
        Displacement coordinates of reference frame.
    xyz_displ2 : np.ndarray, shape=(n_frames, n_atoms_displ_padded, 3), dtype=float32
        Displacement coordinates of frames to iterate over.
    n_atoms_align : int
        The number of true atoms in the alignment coordinates
    n_atoms_displ : int
        The number of true atoms in the displacement coordinates
    frame : int
        Index of the desired reference frame in xyz_align1 / xyz_displ1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDs between xyz_displ1[frame] and all of xyz_displ2, after the alignment that maximizes
        the overlap between xyz_align1[frame] and all of xyz_align2
    rotations: np.ndarray, shape=(len(xyz2), 3, 3)
        The rotation matrix for each of the alignments
    """
    assert (xyz_align1.shape[2] == 3) and (xyz_align2.shape[2] == 3)
    assert (xyz_displ1.shape[2] == 3) and (xyz_displ2.shape[2] == 3)
    if not ((xyz_align1.shape[1] % 4 == 0) & (xyz_align2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz_align1.shape[1], xyz_align2.shape[1]))
    if not ((xyz_displ1.shape[1] % 4 == 0) & (xyz_displ2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz_displ1.shape[1], xyz_displ2.shape[1]))
    if not xyz_align1.shape[0] == xyz_displ1.shape[0]:
        raise ValueError("xyz_align1 and xyz_displ1 must contain the same number of frames")
    if not xyz_align2.shape[0] == xyz_displ2.shape[0]:
        raise ValueError("xyz_align2 and xyz_displ2 must contain the same number of frames")
    if frame >= xyz_align1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz_align1.shape[0]))

    cdef float msd
    cdef Py_ssize_t i
    cdef int n_frames = xyz_align2.shape[0]
    cdef int n_align_atoms_padded = xyz_align1.shape[1]
    cdef int n_displ_atoms_padded = xyz_displ1.shape[1]

    cdef float[:, :, ::1] rot = np.zeros((n_frames, 3, 3), dtype=np.float32)
    cdef float[:] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd_atom_major(n_atoms_align, n_align_atoms_padded, &xyz_align1[frame, 0, 0], &xyz_align2[i, 0, 0], g_align1[frame], g_align2[i], 1, &rot[i, 0, 0])
            msd = rot_msd_atom_major(n_atoms_displ, n_displ_atoms_padded, &xyz_displ1[frame, 0, 0], &xyz_displ2[i, 0, 0], &rot[i, 0, 0])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd_atom_major(n_atoms_align, n_align_atoms_padded, &xyz_align1[frame, 0, 0], &xyz_align2[i, 0, 0], g_align1[frame], g_align2[i], 1, &rot[i, 0, 0])
            msd = rot_msd_atom_major(n_atoms_displ, n_displ_atoms_padded, &xyz_displ1[frame, 0, 0], &xyz_displ2[i, 0, 0], &rot[i, 0, 0])
            distances[i] = sqrtf(msd)

    return np.array(distances, copy=False), np.array(rot, copy=False)

