# Copyright 2012-present mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

# This code was written by Imran S. Haque, with some contributes by
# Robert McGibbon


##############################################################################
# Imports
##############################################################################

import numpy as np
from mdtraj import _rmsd
from mdtraj.utils.six.moves import xrange

##############################################################################
# Globals
##############################################################################

__all__ = ['rmsd_cache', 'RMSDCache', 'align_array']

##############################################################################
# Code
##############################################################################


def rmsd_cache(trajectory, major='axis'):
    """Create a specialized copy of a trajectory's cartesian coordinates for fast RMSD calculations.

    Parameters
    ----------
    trajectory : md.Trajectory
        A Trajectory object, containing the cartesian coordinates of the system.
    major : {'atom' or 'axis'}
        Resulting memory layout for the coordinates array. Axis-major ordering
        performs better for typical structure sizes. Axis-major ordering has
        all the x-coordinates, followed by all the y-coordinates, followed by
        all the z-coordinates. Atom-major ordering has each atom's coordinates
        in atom order.

    Returns
    -------
    cache : RMSDCache
        The returned data structured contains the coordinate data from the
        trajectory, prepared for the RMSD calculation. This includes alignment
        to the appropriate byte boundaries, and other low-level stuff.

    Notes
    -----
    This operation will make a copy of the cartesian coordinates. This can be
    problematic if you have a lot of data (and limited memory). In that case,
    you can construct the RMSDCache directory from its constructor in a
    copy-free manner.

    Examples
    --------
    >>> import mdtraj as md
    >>> t = md.load('trajectory.h5')                          # doctest: +SKIP
    >>> c = md.rmsd_cache(t)                                  # doctest: +SKIP
    >>> c.rmsds_to(r, 0)                                      # doctest: +SKIP
    array([0.0001953,  0.1906953, 0.37336711, ...,  0.29543064,
           0.68428138, 0.35189939], dtype=float32)

    See Also
    --------
    mdtraj.rmsd.RMSDCache : The constructed RMSDCache object. Note that you can manually construct an RMSDCache without a trajectory, using that classes's constructor.
    """

    if major == 'atom':
        aligned = align_array(trajectory.xyz, major)
    elif major == 'axis':
        aligned = _allocate_aligned_array((trajectory.n_frames, 3, trajectory.n_atoms), major)
        aligned[:, :, 0:trajectory.n_atoms] = np.swapaxes(trajectory.xyz, 1, 2)
    else:
        raise ValueError("Must specify 'atom' or 'axis' major ordering")

    return RMSDCache(aligned, major, trajectory.n_atoms)


def _allocate_aligned_array(shape, major):
    '''Allocate an array with 16-byte alignment and padding for IRMSD routines.

    Parameters
    ----------
    shape : tuple
        shape of the array to allocate. The returned array may be larger than
        requested if the number of atoms (as specified by the shape and
        majority) is not a multiple of four; in this case, the number of
        elements along the atoms dimension of the array will be rounded up to
        the nearest multiple of four.
    major : {'axis', 'atoms'}
        See RMSDCache.__init__` for the definition of axis and atom
        majority.

    Returns
    -------
    array : np.ndarray
        An np.ndarray of three dimensions and float32 dtype, such that the
        first element of each conformation is aligned to a 16-byte boundary,
        the number of elements along the atom axis is the number of atoms (as
        specified in ``shape`` rounded up to the nearest multiple of four),
        and all padding atoms (atom elements not corresponding to an actual
        atom) are set to zero.
    '''
    dtype = np.dtype(np.float32)
    alignment = 16 / dtype.itemsize
    if major not in ('atom', 'axis'):
        raise ValueError("Must specify atom or axis major coordinates")

    if major == 'axis':
        n_confs, n_dims, n_atoms = shape
    else:
        n_confs, n_atoms, n_dims = shape

    n_padded_atoms = (n_atoms + alignment - 1) // alignment * alignment
    n_floats = n_dims * n_padded_atoms * n_confs
    nbytes = n_floats * dtype.itemsize
    buf = np.empty(nbytes + alignment * dtype.itemsize, dtype=np.uint8)
    start_index = -buf.ctypes.data % (alignment * dtype.itemsize)
    aligned = buf[start_index:start_index + nbytes].view(dtype)

    # Reshape the array to the correct majority and zero out the padding atoms
    if major == 'axis':
        aligned = aligned.reshape(n_confs, n_dims, n_padded_atoms)
        aligned[:, :, n_atoms:] = 0
    else:
        aligned = aligned.reshape(n_confs, n_padded_atoms, n_dims)
        aligned[:, n_atoms:, :] = 0

    return aligned


def align_array(coordinates, major):
    '''Given an ndarray of coords, return a copy satisfying IRMSD requirements.

    Parameters
    ----------
    coordinates: np.ndarray, ndim=3
        a 3-D ndarray of coordinates for one or more structures
    major: {'axis', 'atom'}
        Specifies the meaning of each dimension of the `coordinates` argument.
        See `Conformations.__init__` for the definition of atom and axis
        majority.

    Returns
    -------
    coord_copy : np.ndarray
        A copy of the given coordinates such that the copy satisfies the
        alignment, padding, and data-type requirements for IRMSD.
    '''
    aligned = _allocate_aligned_array(coordinates.shape, major)
    if major == 'axis':
        aligned[:, :, :coordinates.shape[2]] = coordinates
    else:
        aligned[:, :coordinates.shape[1], :] = coordinates
    return aligned


class RMSDCache(object):
    """Structure to store coordinates and compute RMSDs between conformations.

    RMSDCache wraps a 3-dimensional numpy array of coordinates, transparently
    handling structure centering and G (matrix trace) computation required to
    use the IRMSD fast-Theobald RMSD routines.

    Note that `RMSDCache` will modify the array of coordinates it
    is given, when those structures are centered!

    Parameters
    ----------
    coordinates : np.ndarray, shape=(M, N, P)
        an M x N x P numpy ndarray of type float32. See `major` for
        definition of dimensions.
    major: {'atom' or 'axis'}
        Specifies the storage format of M x N x P ndarray ``coordinates``.
        If ``major=='axis'``, then coordinates should be of shape
        ``(n_conformations, 3, n_padded_atoms)``. If ``major=='atom'``, then
        coordinates should be of shape ``(n_conformations, n_padded_atoms, 3)``
    n_atoms : int
        the number of actual, not padding, atoms in each structure in the
        array

    Notes
    -----
    There are special restrictions on ``coordinates`` to use the IRMSD
    routines. First, if the number of atoms ``n_atoms`` is not a multiple of
    4, then n_padded_atoms should be the next multiple of 4 that is larger
    than ``n_atoms``, and the corresponding 'padding atoms' in the coordinates
    array must be all-zero. Second, the coordinates array must be aligned
    to a 16-byte boundary.

    If your input array does not satisfy these requirements, you may use
    the ``mdtraj.rmsd.align_array`` function to create a copy meeting them.
    ``mdtraj.rmsd.align_array`` is NOT automatically called in this
    constructor, to avoid silently allocating/copying large memory structures.

    Attributes
    ----------
    n_atoms : int
        The number of real atoms in the conformations
    n_padded_atoms : int
        The number of atoms, including padding atoms which are added
        to the conformations (with coordinates 0) so that the number of atoms
        is a multiple of four.
    n_frames : int
        The number of conformations in the RMSDCache
    n_dims : int
        The number of cartesian dimensions. This is always three.
    major : {'axis', 'atom'}
        Specifies the memory layout of the coordinates.

    See Also
    --------
    mdtraj.rmsd_cache : Convenience function to construct an RMSDCache from a Trajectory

    Examples
    --------
    >>> import mdtraj as md
    >>> from mdtraj.rmsd import RMSDCache
    >>> axis_coords = np.sin(np.arange(30).reshape(2, 3, 5))
    >>> print axis_coords                                     # doctest: +SKIP
     [[[ 0.          0.84147098  0.90929743  0.14112001 -0.7568025 ]
       [-0.95892427 -0.2794155   0.6569866   0.98935825  0.41211849]
       [-0.54402111 -0.99999021 -0.53657292  0.42016704  0.99060736]]
     <BLANKLINE>
      [[ 0.65028784 -0.28790332 -0.96139749 -0.75098725  0.14987721]
       [ 0.91294525  0.83665564 -0.00885131 -0.8462204  -0.90557836]
       [-0.13235175  0.76255845  0.95637593  0.27090579 -0.66363388]]]
    >>> c = RMSDCache(md.rmsd.align_array(axis_coords, 'axis'), 'axis', n_atoms=5)
    >>> c.rmsds_to(c, 0)
    array([ 0.        ,  0.17355414], dtype=float32)
    """

    def __init__(self, coordinates, major, n_atoms):
        """Initialize a `RMSDCache` object.
        """
        if major not in ('atom', 'axis'):
            raise ValueError("Must specify atom or axis major coordinates")
        if coordinates.dtype != np.float32:
            raise ValueError("IRMSD can only handle single-precision float")
        if len(coordinates.shape) != 3:
            raise ValueError("Coordinates must have three dimensions "
                             "(conformations, atoms, axes)")
        self.n_atoms = n_atoms
        self.n_frames = coordinates.shape[0]

        if major == 'axis':
            self.n_dims = coordinates.shape[1]
            self.n_padded_atoms = coordinates.shape[2]
            self.axis_major = True
            self.atom_major = False
        else:
            self.n_dims = coordinates.shape[2]
            self.n_padded_atoms = coordinates.shape[1]
            self.axis_major = False
            self.atom_major = True

        self.major = major
        self.cords = coordinates
        self._g = None
        self._centered = False
        assert self.axis_major ^ self.atom_major

        if self.n_dims != 3:
            raise ValueError("IRMSD only supports operation in 3 dimensions")
        if self.n_padded_atoms % 4 != 0:
            raise ValueError("(Padded) number of atoms must be a multiple "
                             "of 4 to use IRMSD routines. You supplied "
                             "self.n_padded_atoms=%s" % self.n_padded_atoms)
        if self.cords.ctypes.data % 16 != 0:
            raise ValueError("Coordinate array must be aligned to a 16-byte "
                             "boundary to use IRMSD routines")

        return

    def __len__(self):
        return self.n_frames

    def _center(self):
        """Transform conformations so that each is centered about 0.

        This function is automatically called if necessary for an alignment to
        proceed, but is exposed in case the user wants to center structures at
        a different time.

        Modifies data in-place since it might be very large.
        """
        for ci in xrange(self.n_frames):
            if self.atom_major:
                centroid = np.mean(self.cords[ci, :self.n_atoms, :], axis=0) \
                             .reshape(1, 3)
                repcent = np.tile(centroid, (self.n_atoms, 1))
                self.cords[ci, :self.n_atoms, :] -= repcent
            else:
                centroid = np.mean(self.cords[ci, :, :self.n_atoms], axis=1) \
                             .reshape(3, 1)
                repcent = np.tile(centroid, (1, self.n_atoms))
                self.cords[ci, :, :self.n_atoms] -= repcent
        self._centered = True
        return

    def _traces(self):
        """Conformation traces. These are also known as the "G values"

        For a structure S made of column vectors Sx, Sy, Sz representing the x,
        y, and z coordinates of each atom in the structure, G(S) = tr(S'S) =
        dot(x,x) + dot(y,y) + dot(z,z). This quantity is related to the radius
        of gyration and is needed in the Theobald RMSD computation.
        """
        if self._g is None:
            self._compute_traces()
        return self._g

    def _compute_traces(self):
        if not self._centered:
            self._center()
        self._g = np.zeros((self.n_frames,), dtype=np.float32)
        for i in xrange(self.n_frames):
            for j in xrange(self.n_dims):
                if self.axis_major:
                    self._g[i] += np.dot(self.cords[i, j, :],
                                         self.cords[i, j, :])
                else:
                    self._g[i] += np.dot(self.cords[i, :, j],
                                         self.cords[i, :, j])
        return

    def rmsds_to(self, other_cache, other_index):
        """Compute RMSD of all conformations to a reference conformation.

        The underlying computation uses OpenMP and so will automatically
        parallelize across cores; it parallelizes over independent
        conformations in this Conformations object.

        To be compared, two Conformations objects must have the same atom/axis
        majority, the same number of atoms, and the same number of padded
        atoms. If the structures have not been centered or G values have not
        been calculated, they will be transparently computed here before the
        RMSDs are computed, and the results saved for future invocations.

        Parameters
        ----------
        other_cache : RMSDCache
            For each conformation in this RMSDCache object compute the RMSD
            to a particular 'reference' conformation in another RMSDCache
            object ``other_cache``, identified by index ``other_index``.
        other_index : int
            The index of the conformation in ``other_cache`` to measure
            distances to

        Notes
        -----
        This function uses OpenMP to parallelize the calculation across
        multiple cores. To control the number of threads launched by OpenMP,
        you can set the environment variable ``OMP_NUM_THREADS``.

        Returns
        -------
        rmsds : np.ndarray, shape=(n_frames,)
            A 1-D numpy array of the RMSDs.
        """
        if other_cache.major != self.major:
            raise ValueError("Cannot align two conformation sets of differing "
                             "atom/axis majority")
        if other_cache.n_atoms != self.n_atoms:
            raise ValueError("Cannot align two conformation sets of differing "
                             "number of atoms")
        if other_cache.n_padded_atoms != self.n_padded_atoms:
            raise ValueError("Cannot align two conformation sets of differing "
                             "number of padded atoms")

        if self.axis_major:
            return _rmsd.getMultipleRMSDs_axis_major(other_cache.cords, self.cords,
                        other_cache._traces(), self._traces(), self.n_atoms,
                        other_index, parallel=True)
        elif self.atom_major:
            return _rmsd.getMultipleRMSDs_atom_major(other_cache.cords, self.cords,
                        other_cache._traces(), self._traces(), self.n_atoms,
                        other_index, parallel=True)
        raise RuntimeError()
