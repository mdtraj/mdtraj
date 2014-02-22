##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A. Beauchamp, TJ Lane, Joshua Adelman, Lee-Ping Wang
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

from __future__ import print_function, division
import os
import warnings
import functools
from copy import deepcopy
import numpy as np

from mdtraj import (DCDTrajectoryFile, BINPOSTrajectoryFile, XTCTrajectoryFile,
                    TRRTrajectoryFile, HDF5TrajectoryFile, NetCDFTrajectoryFile,
                    LH5TrajectoryFile, PDBTrajectoryFile, MDCRDTrajectoryFile,
                    ArcTrajectoryFile, Topology)
from mdtraj.utils import unitcell, ensure_type, convert, cast_indices
from mdtraj.utils.six.moves import xrange
from mdtraj.utils.six import PY3
from mdtraj import _rmsd
from mdtraj import _FormatRegistry

##############################################################################
# Globals
##############################################################################

__all__ = ['open', 'load', 'iterload', 'load_frame', 'Trajectory']

##############################################################################
# Utilities
##############################################################################


def _assert_files_exist(filenames):
    """Throw an IO error if files don't exist

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings to check
    """
    if isinstance(filenames, str):
        filenames = [filenames]
    for fn in filenames:
        if not (os.path.exists(fn) and os.path.isfile(fn)):
            raise IOError('No such file: %s' % fn)


def _parse_topology(top):
    """Get the topology from a argument of indeterminate type
    If top is a string, we try loading a pdb, if its a trajectory
    we extract its topology.
    """

    if isinstance(top, str) and (os.path.splitext(top)[1] in ['.pdb', '.h5','.lh5']):
        topology = load_frame(top, 0).topology
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, Topology):
        topology = top
    else:
        raise TypeError('A topology is required. You supplied top=%s' % top)

    return topology



##############################################################################
# Utilities
##############################################################################


def open(filename, mode='r', force_overwrite=True, **kwargs):
    """Open a trajectory file-like object

    This factor function returns an instance of an open file-like
    object capable of reading/writing the trajectory (depending on
    'mode'). It does not actually load the trajectory from disk or
    write anything.

    Parameters
    ----------
    filename : str
        Path to the trajectory file on disk
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for
        write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?

    Other Parameters
    ----------------
    kwargs : dict
        Other keyword parameters are passed directly to the file object

    Returns
    -------
    fileobject : object
        Open trajectory file, whose type is determined by the filename
        extension

    See Also
    --------
    load, ArcTrajectoryFile, BINPOSTrajectoryFile, DCDTrajectoryFile,
    HDF5TrajectoryFile, LH5TrajectoryFile, MDCRDTrajectoryFile,
    NetCDFTrajectoryFile, PDBTrajectoryFile, TRRTrajectoryFile,
    XTCTrajectoryFile

    """
    extension = os.path.splitext(filename)[1]
    try:
        loader = _FormatRegistry.fileobjects[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files with extensions in %s'
                      % (filename, extension, _FormatRegistry.fileobjects.keys()))
    return loader(filename, mode=mode, force_overwrite=force_overwrite, **kwargs)


def load_frame(filename, index, top=None, atom_indices=None):
    """Load a single frame from a trajectory file

    Parameters
    ----------
    filename : str
        Path to the trajectory file on disk
    index : int
        Load the `index`-th frame from the specified file
    top : {str, Trajectory, Topology}
        Most trajectory formats do not contain topology information. Pass in
        either the path to a RCSB PDB file, a trajectory, or a topology to
        supply this information.
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. These indices are zero-based (not 1 based, as used by the PDB
        format).
        
    Examples
    --------
    >>> import mdtraj as md
    >>> first_frame = md.load_frame('traj.h5', 0)
    >>> print first_frame
    <mdtraj.Trajectory with 1 frames, 22 atoms>

    See Also
    --------
    load, load_frame

    Returns
    -------
    trajectory : md.Trajectory
        The resulting conformation, as an md.Trajectory object containing
        a single frame.
    """
    _assert_files_exist(filename)
    extension = os.path.splitext(filename)[1]
    try:
        loader = _FormatRegistry.loaders[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files with extensions in %s'
                      % (filename, extension, _FormatRegistry.loaders.keys()))

    kwargs = {'top': top, 'atom_indices': atom_indices}
    if loader.__name__ in ['load_hdf5', 'load_pdb']:
        kwargs.pop('top', None)

    return loader(filename, frame=index, **kwargs)


def load(filename_or_filenames, discard_overlapping_frames=False, **kwargs):
    """Load a trajectory from one or more files on disk.

    This function dispatches to one of the specialized trajectory loaders based
    on the extension on the filename. Because different trajectory formats save
    different information on disk, the specific keyword argument options supported
    depend on the specific loaded.

    Parameters
    ----------
    filename_or_filenames : {str, list of strings}
        Filename or list of filenames containing trajectory files of a single format.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them

    Other Parameters
    ----------------
    top : {str, Trajectory, Topology}
        Most trajectory formats do not contain topology information. Pass in
        either the path to a RCSB PDB file, a trajectory, or a topology to
        supply this information. This option is not required for the .h5, .lh5,
        and .pdb formats, which already contain topology information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it
        requires an extra copy, but will save memory.

    See Also
    --------
    load_frame, iterload

    Examples
    --------
    >>> import mdtraj as md
    >>> traj = md.load('output.xtc', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>
    
    >>> traj2 = md.load('output.xtc', stride=2, top='topology.pdb')
    >>> print traj2
    <mdtraj.Trajectory with 250 frames, 423 atoms at 0x11136e410>
    
    >>> traj3 = md.load_hdf5('output.xtc', atom_indices=[0,1] top='topology.pdb')
    >>> print traj3
    <mdtraj.Trajectory with 500 frames, 2 atoms at 0x18236e4a0>
    
    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """

    _assert_files_exist(filename_or_filenames)

    if "top" in kwargs:  # If applicable, pre-loads the topology from PDB for major performance boost.
        kwargs["top"] = _parse_topology(kwargs["top"])

    # grab the extension of the filename
    if isinstance(filename_or_filenames, str):  # If a single filename
        extension = os.path.splitext(filename_or_filenames)[1]
        filename = filename_or_filenames
    else:  # If multiple filenames, take the first one.
        extensions = [os.path.splitext(filename_i)[1] for filename_i in filename_or_filenames]
        if len(set(extensions)) != 1:
            raise(TypeError("All filenames must have same extension!"))
        else:
            t = [load(f, **kwargs) for f in filename_or_filenames]
            # we know the topology is equal because we sent the same topology kwarg
            # in, so there's no reason to spend extra time checking
            return t[0].join(t[1:], discard_overlapping_frames=discard_overlapping_frames,
                             check_topology=False)

    try:
        #loader = _LoaderRegistry[extension][0]
        loader = _FormatRegistry.loaders[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files '
                      'with extensions in %s' % (filename, extension, _FormatRegistry.loaders.keys()))

    if loader.__name__ in ['load_hdf5', 'load_pdb', 'load_lh5']:
        if 'top' in kwargs:
            warnings.warn('top= kwarg ignored since file contains topology information')
        # this is a little hack that makes calling load() more predicable. since
        # most of the loaders take a kwargs "top" except for load_hdf5, (since
        # it saves the topology inside the file), we often end up calling
        # load_hdf5 via this function with the top kwarg specified. but then
        # there would be a signature binding error. it's easier just to ignore
        # it.
        kwargs.pop('top', None)

    return loader(filename, **kwargs)


def iterload(filename, chunk=100, **kwargs):
    """An iterator over a trajectory from one or more files on disk, in fragments

    This may be more memory efficient than loading an entire trajectory at
    once

    Parameters
    ----------
    filename : str
        Path to the trajectory file on disk
    chunk : int
        Number of frames to load at once from disk per iteration.

    Other Parameters
    ----------------
    top : {str, Trajectory, Topology}
        Most trajectory formats do not contain topology information. Pass in
        either the path to a RCSB PDB file, a trajectory, or a topology to
        supply this information. This option is not required for the .h5, .lh5,
        and .pdb formats, which already contain topology information.
    stride : int, default=None
        Only read every stride-th frame.
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it
        requires an extra copy, but will save memory.

    See Also
    --------
    load, load_frame
        
    Examples
    --------
    >>> import mdtraj as md
    >>> for chunk in md.iterload('output.xtc', top='topology.pdb')
    ...    print chunk
    <mdtraj.Trajectory with 100 frames, 423 atoms at 0x110740a90>
    <mdtraj.Trajectory with 100 frames, 423 atoms at 0x110740a90>
    <mdtraj.Trajectory with 100 frames, 423 atoms at 0x110740a90>
    <mdtraj.Trajectory with 100 frames, 423 atoms at 0x110740a90>
    <mdtraj.Trajectory with 100 frames, 423 atoms at 0x110740a90>
    """
    stride = kwargs.get('stride', 1)
    atom_indices = cast_indices(kwargs.get('atom_indices', None))
    if chunk % stride != 0:
        raise ValueError('Stride must be a divisor of chunk. stride=%d does not go '
                         'evenly into chunk=%d' % (stride, chunk))

    if filename.endswith('.h5'):
        if 'top' in kwargs:
            warnings.warn('top= kwarg ignored since file contains topology information')
        with HDF5TrajectoryFile(filename) as f:
            if atom_indices is None:
                topology = f.topology
            else:
                topology = f.topology.subset(atom_indices)

            while True:
                data = f.read(chunk*stride, stride=stride, atom_indices=atom_indices)
                if data == []:
                    raise StopIteration()
                convert(data.coordinates, f.distance_unit, Trajectory._distance_unit, inplace=True)
                convert(data.cell_lengths, f.distance_unit, Trajectory._distance_unit, inplace=True)
                yield Trajectory(xyz=data.coordinates, topology=topology,
                                 time=data.time, unitcell_lengths=data.cell_lengths,
                                 unitcell_angles=data.cell_angles)

    if filename.endswith('.lh5'):
        if 'top' in kwargs:
            warnings.warn('top= kwarg ignored since file contains topology information')
        with LH5TrajectoryFile(filename) as f:
            if atom_indices is None:
                topology = f.topology
            else:
                topology = f.topology.subset(atom_indices)

            ptr = 0
            while True:
                xyz = f.read(chunk*stride, stride=stride, atom_indices=atom_indices)
                if len(xyz) == 0:
                    raise StopIteration()
                convert(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)
                time = np.arange(ptr, ptr+len(xyz)*stride, stride)
                ptr += len(xyz)*stride
                yield Trajectory(xyz=xyz, topology=topology, time=time)

    elif filename.endswith('.xtc'):
        topology = _parse_topology(kwargs.get('top', None))
        with XTCTrajectoryFile(filename) as f:
            while True:
                xyz, time, step, box = f.read(chunk*stride, stride=stride, atom_indices=atom_indices)
                if len(xyz) == 0:
                    raise StopIteration()
                convert(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)
                convert(box, f.distance_unit, Trajectory._distance_unit, inplace=True)
                trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
                trajectory.unitcell_vectors = box
                yield trajectory

    elif filename.endswith('.dcd'):
        topology = _parse_topology(kwargs.get('top', None))
        with DCDTrajectoryFile(filename) as f:
            ptr = 0
            while True:
                # for reasons that I have not investigated, dcdtrajectory file chunk and stride
                # together work like this method, but HDF5/XTC do not.
                xyz, box_length, box_angle = f.read(chunk, stride=stride, atom_indices=atom_indices)
                if len(xyz) == 0:
                    raise StopIteration()
                convert(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)
                convert(box_length, f.distance_unit, Trajectory._distance_unit, inplace=True)
                time = np.arange(ptr, ptr+len(xyz)*stride, stride)
                ptr += len(xyz)*stride
                yield Trajectory(xyz=xyz, topology=topology, time=time, unitcell_lengths=box_length,
                                 unitcell_angles=box_angle)

    else:
        t = load(filename, **kwargs)
        for i in range(0, len(t), chunk):
            yield t[i:i+chunk]


class Trajectory(object):
    """Container object for a molecular dynamics trajectory

    A Trajectory represents a collection of one or more molecular structures,
    generally (but not necessarily) from a molecular dynamics trajectory. The
    Trajectory stores a number of fields describing the system through time,
    including the cartesian coordinates of each atoms (``xyz``), the topology
    of the molecular system (``topology``), and information about the
    unitcell if appropriate (``unitcell_vectors``, ``unitcell_length``,
    ``unitcell_angles``).

    A Trajectory should generally be constructed by loading a file from disk.
    Trajectories can be loaded from (and saved to) the PDB, XTC, TRR, DCD,
    binpos, NetCDF or MDTraj HDF5 formats.

    Trajectory supports fancy indexing, so you can extract one or more frames
    from a Trajectory as a separate trajectory. For example, to form a
    trajectory with every other frame, you can slice with ``traj[::2]``.

    Trajectory uses the nanometer, degree & picosecond unit system.

    Examples
    --------
    >>> # loading a trajectory
    >>> import mdtraj as md
    >>> md.load('trajectory.xtc', top='native.pdb')
    <mdtraj.Trajectory with 1000 frames, 22 atoms at 0x1058a73d0>

    >>> # slicing a trajectory
    >>> t = md.load('trajectory.h5')
    >>> print(t)
    <mdtraj.Trajectory with 100 frames, 22 atoms>
    >>> print(t[::2])
    <mdtraj.Trajectory with 50 frames, 22 atoms>

    >>> # calculating the average distance between two atoms
    >>> import mdtraj as md
    >>> import numpy as np
    >>> t = md.load('trajectory.h5')
    >>> np.mean(np.sqrt(np.sum((t.xyz[:, 0, :] - t.xyz[:, 21, :])**2, axis=1)))

    See Also
    --------
    mdtraj.load : High-level function that loads files and returns an ``md.Trajectory``

    Attributes
    ----------
    n_frames : int
    n_atoms : int
    n_residues : int
    time : np.ndarray, shape=(n_frames,)
    timestep : float
    topology : md.Topology
    top : md.Topology
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
    unitcell_vectors : {np.ndarray, shape=(n_frames, 3, 3), None}
    unitcell_lengths : {np.ndarray, shape=(n_frames, 3), None}
    unitcell_angles : {np.ndarray, shape=(n_frames, 3), None}
    """

    # this is NOT configurable. if it's set to something else, things will break
    # (thus why I make it private)
    _distance_unit = 'nanometers'

    @property
    def topology(self):
        """Topology of the system, describing the organization of atoms into residues, bonds, etc

        Returns
        -------
        topology : md.Topology
            The topology object, describing the organization of atoms into
            residues, bonds, etc
        """

        return self._topology

    @topology.setter
    def topology(self, value):
        "Set the topology of the system, describing the organization of atoms into residues, bonds, etc"
        # todo: more typechecking
        self._topology = value

    @property
    def n_frames(self):
        """Number of frames in the trajectory

        Returns
        -------
        n_frames : int
            The number of frames in the trajectory
        """
        return self._xyz.shape[0]

    @property
    def n_atoms(self):
        """Number of atoms in the trajectory

        Returns
        -------
        n_atoms : int
            The number of atoms in the trajectory
        """
        return self._xyz.shape[1]

    @property
    def n_residues(self):
        """Number of residues (amino acids) in the trajectory

        Returns
        -------
        n_residues : int
            The number of residues in the trajectory's topology
        """
        if self.top is None:
            return 0
        return sum([1 for r in self.top.residues])

    @property
    def top(self):
        """Alias for self.topology, describing the organization of atoms into residues, bonds, etc

        Returns
        -------
        topology : md.Topology
            The topology object, describing the organization of atoms into
            residues, bonds, etc
        """
        return self._topology

    @property
    def timestep(self):
        """Timestep between frames, in picoseconds

        Returns
        -------
        timestep : float
            The timestep between frames, in picoseconds.
        """
        return self._time[1] - self._time[0]

    @property
    def time(self):
        """The simulation time corresponding to each frame, in picoseconds

        Returns
        -------
        time : np.ndarray, shape=(n_frames,)
            The simulation time corresponding to each frame, in picoseconds
        """
        return self._time

    @time.setter
    def time(self, value):
        "Set the simulation time corresponding to each frame, in picoseconds"
        if isinstance(value, list):
            value = np.array(value)

        if np.isscalar(value) and self.n_frames == 1:
            value = np.array([value])
        elif not value.shape == (self.n_frames,):
            raise ValueError('Wrong shape. Got %s, should be %s' % (value.shape,
                (self.n_frames)))

        self._time = value

    @property
    def unitcell_vectors(self):
        """The vectors that define the shape of the unit cell in each frame

        Returns
        -------
        vectors : np.ndarray, shape(n_frames, 3, 3)
            Vectors definiing the shape of the unit cell in each frame.
            The semantics of this array are that the shape of the unit cell
            in frame ``i`` are given by the three vectors, ``value[i, 0, :]``,
            ``value[i, 1, :]``, and ``value[i, 2, :]``.
        """
        if self._unitcell_lengths is None or self._unitcell_angles is None:
            return None

        v1, v2, v3 = unitcell.lengths_and_angles_to_box_vectors(
            self._unitcell_lengths[:, 0],  # a
            self._unitcell_lengths[:, 1],  # b
            self._unitcell_lengths[:, 2],  # c
            self._unitcell_angles[:, 0],   # alpha
            self._unitcell_angles[:, 1],   # beta
            self._unitcell_angles[:, 2],   # gamma
        )
        return np.swapaxes(np.dstack((v1, v2, v3)), 1, 2)

    @unitcell_vectors.setter
    def unitcell_vectors(self, vectors):
        """Set the three vectors that define the shape of the unit cell

        Parameters
        ----------
        vectors : tuple of three arrays, each of shape=(n_frames, 3)
            The semantics of this array are that the shape of the unit cell
            in frame ``i`` are given by the three vectors, ``value[i, 0, :]``,
            ``value[i, 1, :]``, and ``value[i, 2, :]``.
        """
        if vectors is None:
            self._unitcell_lengths = None
            self._unitcell_angles = None
            return

        if not len(vectors) == len(self):
            raise TypeError('unitcell_vectors must be the same length as '
                            'the trajectory. you provided %s' % vectors)

        v1 = vectors[:, 0, :]
        v2 = vectors[:, 1, :]
        v3 = vectors[:, 2, :]
        a, b, c, alpha, beta, gamma = unitcell.box_vectors_to_lengths_and_angles(v1, v2, v3)

        self._unitcell_lengths = np.vstack((a, b, c)).T
        self._unitcell_angles =  np.vstack((alpha, beta, gamma)).T

    @property
    def unitcell_lengths(self):
        """Lengths that define the shape of the unit cell in each frame.

        Returns
        -------
        lengths : {np.ndarray, shape=(n_frames, 3), None}
            Lengths of the unit cell in each frame, in nanometers, or None
            if the Trajectory contains no unitcell information.
        """
        return self._unitcell_lengths

    @property
    def unitcell_angles(self):
        """Angles that define the shape of the unit cell in each frame.

        Returns
        -------
        lengths : np.ndarray, shape=(n_frames, 3)
            The angles between the three unitcell vectors in each frame,
            ``alpha``, ``beta``, and ``gamma``. ``alpha' gives the angle
            between vectors ``b`` and ``c``, ``beta`` gives the angle between
            vectors ``c`` and ``a``, and ``gamma`` gives the angle between
            vectors ``a`` and ``b``. The angles are in degrees.
        """
        return self._unitcell_angles

    @unitcell_lengths.setter
    def unitcell_lengths(self, value):
        """Set the lengths that define the shape of the unit cell in each frame

        Parameters
        ----------
        value : np.ndarray, shape=(n_frames, 3)
            The distances ``a``, ``b``, and ``c`` that define the shape of the
            unit cell in each frame, or None
        """
        self._unitcell_lengths = ensure_type(value, np.float32, 2,
            'unitcell_lengths', can_be_none=True, shape=(len(self), 3),
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)

    @unitcell_angles.setter
    def unitcell_angles(self, value):
        """Set the lengths that define the shape of the unit cell in each frame

        Parameters
        ----------
        value : np.ndarray, shape=(n_frames, 3)
            The angles ``alpha``, ``beta`` and ``gamma`` that define the
            shape of the unit cell in each frame. The angles should be in
            degrees.
        """
        self._unitcell_angles = ensure_type(value, np.float32, 2,
            'unitcell_angles', can_be_none=True, shape=(len(self), 3),
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)

    @property
    def xyz(self):
        """Cartesian coordinates of each atom in each simulation frame

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            A three dimensional numpy array, with the cartesian coordinates
            of each atoms in each frame.
        """
        return self._xyz

    @xyz.setter
    def xyz(self, value):
        "Set the cartesian coordinates of each atom in each simulation frame"
        if self.top is not None:
            # if we have a topology and its not None
            shape = (None, self.topology._numAtoms, 3)
        else:
            shape = (None, None, 3)

        value = ensure_type(value, np.float32, 3, 'xyz', shape=shape,
                            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        self._xyz = value
        self._rmsd_traces = None

    def __len__(self):
        return self.n_frames

    def __add__(self, other):
        "Concatenate two trajectories"
        return self.join(other)

    def __str__(self):
        return "<mdtraj.Trajectory with %d frames, %d atoms>" % (self.n_frames, self.n_atoms)

    def __repr__(self):
        return "<mdtraj.Trajectory with %d frames, %d atoms at 0x%02x>" % (self.n_frames, self.n_atoms, id(self))

    def superpose(self, reference, frame=0, atom_indices=None, parallel=True):
        """Superpose each conformation in this trajectory upon a reference

        Parameters
        ----------
        reference : md.Trajectory
            For each conformation in this trajectory, aligned to a particular
            reference conformation in another trajectory object.
        frame : int
            The index of the conformation in `reference` to align to.
        atom_indices : array_like, or None
            The indices of the atoms to superpose. If not
            supplied, all atoms will be used.
        parallel : bool
            Use OpenMP to run the superposition in parallel over multiple cores

        Returns
        -------
        self
        """
        if atom_indices is None:
            atom_indices = slice(None)

        n_frames = self.xyz.shape[0]
        self_align_xyz = np.asarray(self.xyz[:, atom_indices, :], order='c')
        self_displace_xyz = np.asarray(self.xyz, order='c')
        ref_align_xyz = np.array(reference.xyz[frame, atom_indices, :], copy=True, order='c').reshape(1, -1, 3)

        offset = np.mean(self_align_xyz, axis=1, dtype=np.float64).reshape(n_frames, 1, 3)
        self_align_xyz -= offset
        if self_align_xyz.ctypes.data != self_displace_xyz.ctypes.data:
            # when atom_indices is None, these two arrays alias the same memory
            # so we only need to do the centering once
            self_displace_xyz -= offset

        ref_offset = ref_align_xyz[0].astype('float64').mean(0)
        ref_align_xyz[0] -= ref_offset

        self_g = np.einsum('ijk,ijk->i', self_align_xyz, self_align_xyz)
        ref_g = np.einsum('ijk,ijk->i', ref_align_xyz , ref_align_xyz)

        _rmsd.superpose_atom_major(
            ref_align_xyz, self_align_xyz, ref_g, self_g, self_displace_xyz,
            0, parallel=parallel)

        self.xyz = self_displace_xyz + ref_offset
        return self

    def join(self, other, check_topology=True, discard_overlapping_frames=False):
        """Join two trajectories together along the time/frame axis.

        This method joins trajectories along the time axis, giving a new trajectory
        of length equal to the sum of the lengths of `self` and `other`.
        It can also be called by using `self + other`

        Parameters
        ----------
        other : Trajectory or list of Trajectory
            One or more trajectories to join with this one. These trajectories
            are *appended* to the end of this trajectory.
        check_topology : bool
            Ensure that the topology of `self` and `other` are identical before
            joining them. If false, the resulting trajectory will have the
            topology of `self`.
        discard_overlapping_frames : bool, optional
            If True, compare coordinates at trajectory edges to discard overlapping
            frames.  Default: False.

        See Also
        --------
        stack : join two trajectories along the atom axis
        """

        if isinstance(other, Trajectory):
            other = [other]
        if isinstance(other, list):
            if not all(isinstance(o, Trajectory) for o in other):
                raise TypeError('You can only join Trajectory instances')
            if not all(self.n_atoms == o.n_atoms for o in other):
                raise  ValueError('Number of atoms in self (%d) is not equal '
                          'to number of atoms in other' % (self.n_atoms))
            if check_topology and not all(self.topology == o.topology for o in other):
                raise ValueError('The topologies of the Trajectories are not the same')
            if not all(self._have_unitcell == o._have_unitcell for o in other):
                raise ValueError('Mixing trajectories with and without unitcell')
        else:
            raise TypeError('`other` must be a list of Trajectory. You supplied %d' % type(other))


        # list containing all of the trajs to merge, including self
        trajectories = [self] + other
        if discard_overlapping_frames:
            for i in range(len(trajectories)-1):
                # last frame of trajectory i
                x0 = trajectories[i].xyz[-1]
                # first frame of trajectory i+1
                x1 = trajectories[i + 1].xyz[0]

                # check that all atoms are within 2e-3 nm
                # (this is kind of arbitrary)
                if np.all(np.abs(x1 - x0) < 2e-3):
                    trajectories[i] = trajectories[i][:-1]

        xyz = np.concatenate([t.xyz for t in trajectories])
        time = np.concatenate([t.time for t in trajectories])
        angles = lengths = None
        if self._have_unitcell:
            angles = np.concatenate([t.unitcell_angles for t in trajectories])
            lengths = np.concatenate([t.unitcell_lengths for t in trajectories])

        # use this syntax so that if you subclass Trajectory,
        # the subclass's join() will return an instance of the subclass
        return self.__class__(xyz, deepcopy(self._topology), time=time,
            unitcell_lengths=lengths, unitcell_angles=angles)

    def stack(self, other):
        """Stack two trajectories along the atom axis

        This method joins trajectories along the atom axis, giving a new trajectory
        with a number of atoms equal to the sum of the number of atoms in
        `self` and `other`.

        Notes
        -----
        The resulting trajectory will have the unitcell and time information
        the left operand.

        Examples
        --------
        >>> t1 = md.load('traj1.h5')
        >>> t2 = md.load('traj2.h5')
        >>> # even when t2 contains no unitcell information
        >>> t2.unitcell_vectors = None
        >>> stacked = t1.stack(t2)
        >>> # the stacked trajectory inherits the unitcell information
        >>> # from the first trajectory
        >>> np.all(stacked.unitcell_vectors == t1.unitcell_vectors)
        True

        Parameters
        ----------
        other : Trajectory
            The other trajectory to join

        See Also
        --------
        join : join two trajectories along the time/frame axis.
        """
        if not isinstance(other, Trajectory):
            raise TypeError('You can only stack two Trajectory instances')
        if self.n_frames != other.n_frames:
            raise ValueError('Number of frames in self (%d) is not equal '
                             'to number of frames in other (%d)' % (self.n_frames, other.n_frames))
        if self.topology is not None:
            topology = self.topology.join(other.topology)
        else:
            topology = None

        xyz = np.hstack((self.xyz, other.xyz))
        return self.__class__(xyz=xyz, topology=topology, unitcell_angles=self.unitcell_angles,
                              unitcell_lengths=self.unitcell_lengths, time=self.time)

    def __getitem__(self, key):
        "Get a slice of this trajectory"
        return self.slice(key)

    def slice(self, key, copy=True):
        """Slice trajectory, by extracting one or more frames into a separate object

        This method can also be called using index bracket notation, i.e
        `traj[1] == traj.slice(1)`

        Parameters
        ----------
        key : {int, np.ndarray, slice}
            The slice to take. Can be either an int, a list of ints, or a slice
            object.
        copy : bool, default=True
            Copy the arrays after slicing. If you set this to false, then if
            you modify a slice, you'll modify the original array since they
            point to the same data.
        """
        xyz = self.xyz[key]
        time = self.time[key]
        unitcell_lengths, unitcell_angles = None, None
        if self.unitcell_angles is not None:
            unitcell_angles = self.unitcell_angles[key]
        if self.unitcell_lengths is not None:
            unitcell_lengths = self.unitcell_lengths[key]

        if copy:
            xyz = xyz.copy()
            time = time.copy()
            topology = deepcopy(self._topology)

            if self.unitcell_angles is not None:
                unitcell_angles = unitcell_angles.copy()
            if self.unitcell_lengths is not None:
                unitcell_lengths = unitcell_lengths.copy()

        newtraj = self.__class__(xyz, topology, time, unitcell_lengths=unitcell_lengths,
                                 unitcell_angles=unitcell_angles)
        return newtraj

    def __init__(self, xyz, topology, time=None, unitcell_lengths=None, unitcell_angles=None):
        # install the topology into the object first, so that when setting
        # the xyz, we can check that it lines up (e.g. n_atoms), with the topology
        self.topology = topology
        self.xyz = xyz

        # _rmsd_traces are the inner product of each centered conformation,
        # which are required for computing RMSD. Normally these values are
        # calculated on the fly in the cython code (rmsd/_rmsd.pyx), but
        # optionally, we enable the use precomputed values which can speed
        # up the calculation (useful for clustering), but potentially be unsafe
        # if self._xyz is modified without a corresponding change to
        # self._rmsd_traces. This array is populated computed by
        # center_conformations, and no other methods should really touch it.
        self._rmsd_traces = None

        # box has no default, it'll just be none normally
        self.unitcell_lengths = unitcell_lengths
        self.unitcell_angles = unitcell_angles

        # time will take the default 1..N
        if time is None:
            time = np.arange(len(self.xyz))
        self.time = time

        if (topology is not None) and (topology._numAtoms != self.n_atoms):
             raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (self.n_atoms, topology._numAtoms))

    def openmm_positions(self, frame):
        """OpenMM-compatable positions of a single frame.

        Examples
        --------
        >>> t = md.load('trajectory.h5')
        >>> context.setPositions(t.openmm_positions(0))

        Parameters
        ----------
        frame : int
            The index of frame of the trajectory that you wish to extract

        Returns
        -------
        positions : list
            The cartesian coordinates of specific trajectory frame, formatted
            for input to OpenMM

        """
        from simtk.openmm import Vec3
        from simtk.unit import nanometer

        Pos = []
        for xyzi in self.xyz[frame]:
            Pos.append(Vec3(xyzi[0], xyzi[1], xyzi[2]))

        return Pos * nanometer

    def openmm_boxes(self, frame):
        """OpenMM-compatable box vectors of a single frame.

        Examples
        --------
        >>> t = md.load('trajectory.h5')
        >>> context.setPeriodicBoxVectors(t.openmm_positions(0))

        Parameters
        ----------
        frame : int
            Return box for this single frame.

        Returns
        -------
        box : tuple
            The periodic box vectors for this frame, formatted for input to
            OpenMM.
        """
        from simtk.openmm import Vec3
        from simtk.unit import nanometer

        vectors = self[frame].unitcell_vectors
        if vectors is None:
            raise ValueError("this trajectory does not contain box size information")

        v1, v2, v3 = vectors
        return (Vec3(*v1), Vec3(*v2), Vec3(*v3)) * nanometer

    @staticmethod
    # im not really sure if the load function should be just a function or a method on the class
    # so effectively, lets make it both?
    def load(filenames, **kwargs):
        """Load a trajectory from disk

        Parameters
        ----------
        filenames : {str, [str]}
            Either a string or list of strings

        Other Parameters
        ----------------
        As requested by the various load functions -- it depends on the extension
        """
        return load(filenames, **kwargs)

    def save(self, filename, **kwargs):
        """Save trajectory to disk, in a format determined by the filename extension

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory. The extension will
            be parsed and will control the format.

        Other Parameters
        ----------------
        lossy : bool
            For .h5 or .lh5, whether or not to use compression.
        no_models: bool
            For .pdb. TODO: Document this?
        force_overwrite : bool
            For .binpos, .xtc, .dcd. If `filename` already exists, overwrite it.
        """
        # grab the extension of the filename
        extension = os.path.splitext(filename)[1]

        savers = {'.xtc': self.save_xtc,
                  '.trr': self.save_trr,
                  '.pdb': self.save_pdb,
                  '.dcd': self.save_dcd,
                  '.h5': self.save_hdf5,
                  '.binpos': self.save_binpos,
                  '.nc': self.save_netcdf,
                  '.crd': self.save_mdcrd,
                  '.mdcrd': self.save_mdcrd,
                  '.ncdf': self.save_netcdf,
                  '.lh5': self.save_lh5,
                  }

        try:
            saver = savers[extension]
        except KeyError:
            raise IOError('Sorry, no saver for filename=%s (extension=%s) '
                          'was found. I can only save files '
                          'with extensions in %s' % (filename, extension, savers.keys()))

        # run the saver, and return whatever output it gives
        return saver(filename, **kwargs)

    def save_hdf5(self, filename, force_overwrite=True):
        """Save trajectory to MDTraj HDF5 format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with HDF5TrajectoryFile(filename, 'w', force_overwrite=True) as f:
            f.write(coordinates=self.xyz, time=self.time,
                    cell_angles=self.unitcell_angles,
                    cell_lengths=self.unitcell_lengths)
            f.topology = self.topology

    def save_pdb(self, filename, force_overwrite=True):
        """Save trajectory to RCSB PDB format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        self._check_valid_unitcell()

        with PDBTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            for i in xrange(self.n_frames):

                if self._have_unitcell:
                    f.write(convert(self._xyz[i], Trajectory._distance_unit, f.distance_unit),
                            self.topology,
                            modelIndex=i,
                            unitcell_lengths=convert(self.unitcell_lengths[i], Trajectory._distance_unit, f.distance_unit),
                            unitcell_angles=self.unitcell_angles[i])
                else:
                    f.write(convert(self._xyz[i], Trajectory._distance_unit, f.distance_unit),
                            self.topology,
                            modelIndex=i)

    def save_xtc(self, filename, force_overwrite=True):
        """Save trajectory to Gromacs XTC format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with XTCTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=self.xyz, time=self.time, box=self.unitcell_vectors)

    def save_trr(self, filename, force_overwrite=True):
        """Save trajectory to Gromacs TRR format

        Notes
        -----
        Only the xyz coordinates and the time are saved, the velocities
        and forces in the trr will be zeros

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with TRRTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=self.xyz, time=self.time, box=self.unitcell_vectors)

    def save_dcd(self, filename, force_overwrite=True):
        """Save trajectory to CHARMM/NAMD DCD format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filenames, if its already there
        """
        self._check_valid_unitcell()
        with DCDTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(convert(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    cell_lengths=convert(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)


    def save_binpos(self, filename, force_overwrite=True):
        """Save trajectory to AMBER BINPOS format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with BINPOSTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(convert(self.xyz, Trajectory._distance_unit, f.distance_unit))


    def save_mdcrd(self, filename, force_overwrite=True):
        """Save trajectory to AMBER mdcrd format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        self._check_valid_unitcell()
        if self._have_unitcell:
            if not np.all(self.unitcell_angles == 90):
                raise ValueError('Only rectilinear boxes can be saved to mdcrd files')

        with MDCRDTrajectoryFile(filename, mode='w', force_overwrite=force_overwrite) as f:
            f.write(convert(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    convert(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit))


    def save_netcdf(self, filename, force_overwrite=True):
        """Save trajectory in AMBER NetCDF format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        self._check_valid_unitcell()
        with NetCDFTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(coordinates=convert(self._xyz, Trajectory._distance_unit, NetCDFTrajectoryFile.distance_unit),
                    time=self.time,
                    cell_lengths=convert(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)

    def save_lh5(self, filename):
        """Save trajectory in deprecated MSMBuilder2 LH5 (lossy HDF5) format.

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        """
        with LH5TrajectoryFile(filename, 'w', force_overwrite=True) as f:
            f.write(coordinates=self.xyz)
            f.topology = self.topology

    def center_coordinates(self, mass_weighted=False):
        """Center each trajectory frame at the origin (0,0,0).

        This method acts inplace on the trajectory.  The centering can
        be either uniformly weighted (mass_weighted=False) or weighted by
        the mass of each atom (mass_weighted=True).

        Parameters
        ----------
        mass_weighted : bool, optional (default = False)
            If True, weight atoms by mass when removing COM.

        Returns
        -------
        self
        """
        if mass_weighted and self.top is not None:
            masses = np.array([a.element.mass for a in self.top.atoms])
            masses /= masses.sum()
            for x in self._xyz:
                x -= (x.astype('float64').T.dot(masses))
        else:
            self._rmsd_traces = _rmsd._center_inplace_atom_major(self._xyz)

        return self

    def restrict_atoms(self, atom_indices):
        """Retain only a subset of the atoms in a trajectory (inplace)

        Deletes atoms not in `atom_indices`, and re-indexes those that remain

        Parameters
        ----------
        atom_indices : list([int])
            List of atom indices to keep.

        Returns
        -------
        self
        """
        if self._topology is not None:
            self._topology = self._topology.subset(atom_indices)
        self._xyz = np.array(self.xyz[:,atom_indices], order='C')
        return self


    def _check_valid_unitcell(self):
        """Do some sanity checking on self.unitcell_lengths and self.unitcell_angles
        """
        if self.unitcell_lengths is not None and self.unitcell_angles is None:
            raise AttributeError('unitcell length data exists, but no angles')
        if self.unitcell_lengths is None and self.unitcell_angles is not None:
            raise AttributeError('unitcell angles data exists, but no lengths')

        if self.unitcell_lengths is not None and np.any(self.unitcell_lengths < 0):
            raise ValueError('unitcell length < 0')

        if self.unitcell_angles is not None and np.any(self.unitcell_angles < 0):
            raise ValueError('unitcell angle < 0')

    @property
    def _have_unitcell(self):
        return self._unitcell_lengths is not None and self._unitcell_angles is not None
