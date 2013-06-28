import numpy as np


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
        See TheobaldConformations.__init__` for the definition of axis and atom
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

    n_padded_atoms = (n_atoms + alignment - 1) / alignment * alignment
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


class TheobaldConformations(object):
    """Structure to store coordinates and compute RMSDs between conformations.

    TheobaldConformations wraps a 3-dimensional numpy array of coordinates,
    transparently handling structure centering and G (matrix trace) computation
    required to use the IRMSD fast-Theobald RMSD routines.

    Note that `TheobaldConformations` will modify the array of coordinates it
    is given, when those structures are centered!

    Parameters
    ----------
    coordinates : np.ndarray, shape=(M, N, P)
        an M x N x P numpy ndarray of type float32. See `major` for
        definition of dimensions.
    major: {'atom' or 'axis'}
        Specifies the storage format of M x N x P ndarray ``coordinates``
        'axis': M = # conformations; N = # dimensions;
                P = # padded atoms
        'atom': M = # conformations; N = # padded atoms;
                P = # dimensions
        note that 'dimensions' must be 3 (points in 3D space)
    n_atoms : int
        the number of actual, not padding, atoms in each structure in the
        array

    Notes
    -----
    There are special restrictions on ``coordinates`` to use the IRMSD
    routines. First, if the number of atoms ``n_atoms` is not a multiple of
    4, then n_padded_atoms should be the next multiple of 4 that is larger
    than n_atoms, and the corresponding 'padding atoms' in the coordinates
    array must be all-zero. Second, the coordinates array must be aligned
    to a 16-byte boundary.

    If your input array does not satisfy these requirements, you may use
    the `align_array` function to create a copy meeting them. `align_array`
    is NOT automatically called in this constructor, to avoid silently
    allocating/copying large memory structures.

    Attributes
    ----------
    n_atoms : int
        The number of real atoms in the conformations
    n_padded_atoms : int
        The number of atoms, including padding atoms which are added
        to the conformations (with coordinates 0) so that the number of atoms
        is a multiple of four.
    n_frames : int
        The number of conformations in the TheobaldConformations
    n_dims : int
        The number of cartesian dimensions. This is always three.
    major : {'axis', 'atom'}
        Specifies the memory layout of the coordinates.
    G : np.ndarray
        Conformation traces
    """

    def __init__(self, coordinates, major, n_atoms):
        """Initialize a `TheobaldConformations` object.
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
        self._G = None
        self._centered = False
        assert self.axis_major ^ self.atom_major

        if self.n_dims != 3:
            raise ValueError("IRMSD only supports operation in 3 dimensions")
        if self.n_padded_atoms % 4 != 0:
            raise ValueError("(Padded) number of atoms must be a multiple "
                             "of 4 to use IRMSD routines")
        if self.cords.ctypes.data % 16 != 0:
            raise ValueError("Coordinate array must be aligned to a 16-byte "
                             "boundary to use IRMSD routines")

        return

    def __len__(self):
        return self.n_frames

    def center(self):
        """Transform conformations so that each is centered about 0.

        This function is automatically called if necessary for an alignment to
        proceed, but is exposed in case the user wants to center structures at
        a different time.

        Modifies data in-place since it might be very large.
        """
        for ci in xrange(self.nconfs):
            if self.atom_major:
                centroid = np.mean(self.cords[ci, :self.natoms, :], axis=0) \
                             .reshape(1, 3)
                repcent = np.tile(centroid, (self.natoms, 1))
                self.cords[ci, :self.natoms, :] -= repcent
            else:
                centroid = np.mean(self.cords[ci, :, :self.natoms], axis=1) \
                             .reshape(3, 1)
                repcent = np.tile(centroid, (1, self.natoms))
                self.cords[ci, :, :self.natoms] -= repcent
        self._centered = True
        return

    @property
    def G(self):
        """Conformation traces

        For a structure S made of column vectors Sx, Sy, Sz representing the x,
        y, and z coordinates of each atom in the structure, G(S) = tr(S'S) =
        dot(x,x) + dot(y,y) + dot(z,z). This quantity is related to the radius
        of gyration and is needed in the Theobald RMSD computation.
        """
        if self._G is None:
            self._compute_g()
        return self._G

    def _compute_g(self):
        if not self._centered:
            self.center()
        self._G = np.zeros((self.nconfs,), dtype=np.float32)
        for i in xrange(self.nconfs):
            for j in xrange(self.ndims):
                if self.axis_major:
                    self._G[i] += np.dot(self.cords[i, j, :],
                                         self.cords[i, j, :])
                else:
                    self._G[i] += np.dot(self.cords[i, :, j],
                                         self.cords[i, :, j])
        return

    def rmsds_to_reference(self, other_confs, ref_idx):
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
        other_confs : TheobaldConformations
            For each conformation in this TheobaldConformations object,
            compute the RMSD to a particular 'reference' conformation in
            another TheobaldConformations object ``other_confs``, identified
            by index ``ref_idx``.
        ref_idx : int
            The ind

        Returns
        -------
        rmsds : np.ndarray, shape=(n_frames,)
            A 1-D numpy array of the RMSDs.
        """
        if other_confs.major != self.major:
            raise ValueError("Cannot align two conformation sets of differing "
                             "atom/axis majority")
        if other_confs.natoms != self.n_atoms:
            raise ValueError("Cannot align two conformation sets of differing "
                             "number of atoms")
        if other_confs.n_padded_atoms != self.n_padded_atoms:
            raise ValueError("Cannot align two conformation sets of differing "
                             "number of padded atoms")

        ref_structure = other_confs.cords[ref_idx, :, :]
        # Getting G will implicitly center structures if needed
        G = self.G
        ref_G = other_confs.G[ref_idx]

        if self.axis_major:
            return rmsdcalc.getMultipleRMSDs_axis_major(
                self.n_atoms, self.n_padded_atoms, self.n_padded_atoms,
                self.cords, ref_structure,
                G, ref_G)
        else:
            return rmsdcalc.getMultipleRMSDs_atom_major(
                self.n_atoms, self.n_padded_atoms,
                self.cords, ref_structure,
                G, ref_G)
