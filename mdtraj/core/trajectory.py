##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2022 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A. Beauchamp, TJ Lane, Joshua Adelman, Lee-Ping Wang, Jason Swails
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
from copy import deepcopy
from collections.abc import Iterable
import numpy as np
import functools

from mdtraj.formats import DCDTrajectoryFile
from mdtraj.formats import BINPOSTrajectoryFile
from mdtraj.formats import XTCTrajectoryFile
from mdtraj.formats import TRRTrajectoryFile
from mdtraj.formats import HDF5TrajectoryFile
from mdtraj.formats import NetCDFTrajectoryFile
from mdtraj.formats import LH5TrajectoryFile
from mdtraj.formats import PDBTrajectoryFile
from mdtraj.formats import MDCRDTrajectoryFile
from mdtraj.formats import DTRTrajectoryFile
from mdtraj.formats import LAMMPSTrajectoryFile
from mdtraj.formats import XYZTrajectoryFile
from mdtraj.formats import GroTrajectoryFile
from mdtraj.formats import TNGTrajectoryFile
from mdtraj.formats import AmberNetCDFRestartFile
from mdtraj.formats import AmberRestartFile

from mdtraj.formats.prmtop import load_prmtop
from mdtraj.formats.psf import load_psf
from mdtraj.formats.mol2 import load_mol2
from mdtraj.formats.gro import load_gro
from mdtraj.formats.arc import load_arc
from mdtraj.formats.hoomdxml import load_hoomdxml
from mdtraj.formats.gsd import write_gsd, load_gsd_topology
from mdtraj.core.topology import Topology
from mdtraj.core.residue_names import _SOLVENT_TYPES
from mdtraj.utils import (ensure_type, in_units_of, lengths_and_angles_to_box_vectors,
                          box_vectors_to_lengths_and_angles, cast_indices,
                          deprecated)
from mdtraj.utils.six.moves import xrange
from mdtraj.utils.six import PY3, string_types
from mdtraj import _rmsd
from mdtraj import FormatRegistry
from mdtraj.geometry import distance
from mdtraj.geometry import _geometry

##############################################################################
# Globals
##############################################################################

__all__ = ['open', 'load', 'iterload', 'load_frame', 'load_topology', 'join',
           'Trajectory']
# supported extensions for constructing topologies
_TOPOLOGY_EXTS = ['.pdb', '.pdb.gz', '.h5','.lh5', '.prmtop', '.parm7', '.prm7',
                  '.psf', '.mol2', '.hoomdxml', '.gro', '.arc', '.hdf5', '.gsd']


##############################################################################
# Utilities
##############################################################################


def _assert_files_exist(filenames):
    """Throw an IO error if files don't exist

    Parameters
    ----------
    filenames : {path-like, [path-like]}
        Path or list of paths to check
    """
    if isinstance(filenames, (string_types, os.PathLike)):
        filenames = [filenames]
    for fn in filenames:
        if not (os.path.exists(fn) and os.path.isfile(fn)):
            raise IOError('No such file: %s' % fn)


def _assert_files_or_dirs_exist(names):
    """Throw an IO error if files don't exist

    Parameters
    ----------
    filenames : {path-like, [path-like]}
        Path or list of paths to check
    """
    if isinstance(names, (string_types, os.PathLike)):
        names = [names]
    for fn in names:
        if not (os.path.exists(fn) and \
                        (os.path.isfile(fn) or os.path.isdir(fn))):
            raise IOError('No such file: %s' % fn)

if PY3:
    def _hash_numpy_array(x):
        hash_value = hash(x.shape)
        hash_value ^= hash(x.strides)
        hash_value ^= hash(x.data.tobytes())
        return hash_value
else:
    def _hash_numpy_array(x):
        writeable = x.flags.writeable
        try:
            x.flags.writeable = False
            hash_value = hash(x.shape)
            hash_value ^= hash(x.strides)
            hash_value ^= hash(x.data)
        finally:
            x.flags.writeable = writeable
        return hash_value


def load_topology(filename, **kwargs):
    """Load a topology

    Parameters
    ----------
    filename : path-like
        Path to a file containing a system topology. The following extensions
        are supported: '.pdb', '.pdb.gz', '.h5','.lh5', '.prmtop', '.parm7',
            '.prm7', '.psf', '.mol2', '.hoomdxml', '.gsd'

    Returns
    -------
    topology : md.Topology
    """
    return _parse_topology(filename, **kwargs)


def _parse_topology(top, **kwargs):
    """Get the topology from a argument of indeterminate type
    If top is a string, we try loading a pdb, if its a trajectory
    we extract its topology.

    Returns
    -------
    topology : md.Topology
    """

    if isinstance(top, (string_types, os.PathLike)):
        ext = _get_extension(top)
    else:
        ext = None  # might not be a string

    if isinstance(top, Topology):
        topology = top
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.pdb', '.pdb.gz', '.pdbx', '.cif', '.h5','.lh5']):
        _traj = load_frame(top, 0, **kwargs)
        topology = _traj.topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.prmtop', '.parm7', '.prm7']):
        topology = load_prmtop(top, **kwargs)
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.psf']):
        topology = load_psf(top, **kwargs)
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.mol2']):
        topology = load_mol2(top, **kwargs).topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.gro']):
        topology = load_gro(top, **kwargs).topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.arc']):
        topology = load_arc(top, **kwargs).topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.hoomdxml']):
        topology = load_hoomdxml(top, **kwargs).topology
    elif isinstance(top, (string_types, os.PathLike)) and (ext in ['.gsd']):
        topology = load_gsd_topology(top, **kwargs)
    elif isinstance(top, (string_types, os.PathLike)):
        raise IOError('The topology is loaded by filename extension, and the '
                      'detected "%s" format is not supported. Supported topology '
                      'formats include %s and "%s".' % (
                          ext, ', '.join(['"%s"' % e for e in _TOPOLOGY_EXTS[:-1]]),
                          _TOPOLOGY_EXTS[-1]))
    else:
        raise TypeError('A topology is required. You supplied top=%s' % str(top))

    return topology

def _get_extension(filename):
    (base, extension) = os.path.splitext(filename)
    if extension == '.gz':
        extension2 = os.path.splitext(base)[1]
        return extension2 + extension
    return extension

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
    filename : path-like
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
    XTCTrajectoryFile, TNGTrajectoryFile

    """
    extension = _get_extension(filename)
    try:
        loader = FormatRegistry.fileobjects[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files with extensions in %s'
                      % (filename, extension, FormatRegistry.fileobjects.keys()))
    return loader(filename, mode=mode, force_overwrite=force_overwrite, **kwargs)


def load_frame(filename, index, top=None, atom_indices=None, **kwargs):
    """Load a single frame from a trajectory file

    Parameters
    ----------
    filename : path-like
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

    extension = _get_extension(filename)
    try:
        loader = FormatRegistry.loaders[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files with extensions in %s'
                      % (filename, extension, FormatRegistry.loaders.keys()))

    kwargs['atom_indices'] = atom_indices
    if extension not in _TOPOLOGY_EXTS:
        kwargs['top'] = top

    if loader.__name__ not in ['load_dtr']:
        _assert_files_exist(filename)
    else:
        _assert_files_or_dirs_exist(filename)

    return loader(filename, frame=index, **kwargs)


def load(filename_or_filenames, discard_overlapping_frames=False, **kwargs):
    """Load a trajectory from one or more files on disk.

    This function dispatches to one of the specialized trajectory loaders based
    on the extension on the filename. Because different trajectory formats save
    different information on disk, the specific keyword argument options supported
    depend on the specific loaded.

    Parameters
    ----------
    filename_or_filenames : {path-like, list of path-like objects}
        Filename or list of filenames containing trajectory files of a single format.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them

    Other Parameters
    ----------------
    top : {path-like, Trajectory, Topology}
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

    >>> traj3 = md.load_hdf5('output.xtc', atom_indices=[0,1], top='topology.pdb')
    >>> print traj3
    <mdtraj.Trajectory with 500 frames, 2 atoms at 0x18236e4a0>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """

    # If a single filename make a list out of it in
    #order to have an easier function later on
    if isinstance(filename_or_filenames, (string_types, os.PathLike)):
        filename_or_filenames = [filename_or_filenames]


    extensions = [_get_extension(f) for f in filename_or_filenames]
    extension = extensions[0]
    #Make the needed checks
    if len(set(extensions)) == 0:
        raise ValueError('No trajectories specified. '
                            'filename_or_filenames was an empty list')
    elif len(set(extensions)) > 1:
        raise TypeError("Each filename must have the same extension. "
                        "Received: %s" % ', '.join(set(extensions)))

    #pre-loads the topology from PDB for major performance boost.
    topkwargs = kwargs.copy()
    #if top is not given try with one of the trajectory files
    topkwargs.pop("top", None)
    topkwargs.pop("atom_indices", None)
    topkwargs.pop("frame", None)
    topkwargs.pop("stride", None)
    topkwargs.pop("start", None)
    kwargs["top"] = _parse_topology(kwargs.get("top", filename_or_filenames[0]), **topkwargs)

    #get the right loader
    try:
        #loader = _LoaderRegistry[extension][0]
        loader = FormatRegistry.loaders[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                    'was found. I can only load files '
                    'with extensions in %s' % (
                        filename_or_filenames[0], extension, FormatRegistry.loaders.keys()))

    if loader.__name__ not in ['load_dtr']:
            _assert_files_exist(filename_or_filenames)
    else:
        _assert_files_or_dirs_exist(filename_or_filenames)


    if extension not in _TOPOLOGY_EXTS:
        # standard_names is a valid keyword argument only for files containing topologies
        kwargs.pop('standard_names', None)

    trajectories = []
    tmp_file = filename_or_filenames[0]
    filename_or_filenames = filename_or_filenames[1:]  # ignore first file
    try:
        # this is a little hack that makes calling load() more predictable. since
        # most of the loaders take a kwargs "top" except for load_hdf5, (since
        # it saves the topology inside the file), we often end up calling
        # load_hdf5 via this function with the top kwarg specified. but then
        # there would be a signature binding error. it's easier just to ignore
        # it.
        #TODO make all the loaders accept a pre parsed topology (top) in order to avoid
        #this part and have a more consistent interface and a faster load function
        t = loader(tmp_file, **kwargs)

    except TypeError as e:

        #Don't want to intercept legit
        #TypeErrors
        if "got an unexpected keyword argument 'top'" not in str(e):
            raise

        warnings.warn('top= kwargs ignored since this file parser does not support it')

        kwargs.pop('top', None)

        t = loader(tmp_file, **kwargs)

    except ValueError as e:

        if 'xyz must be shape' in str(e):

            raise ValueError('The topology and the trajectory files might not contain the same atoms\n'
            'The input topology must contain all atoms even if '
            'you want to select a subset of them with atom_indices'
            ) from e

        raise

    trajectories.append(t)

    # Only do this monkey patching if needed in order not to
    # modify the output topology
    if ('top' in kwargs) and (
        kwargs.get('atom_indices', None) is not None) and (
        len(filename_or_filenames) > 0):

        # In case only a part of the atoms were selected
        # I get the right topology that
        # kwargs['top'].subset shall return
        subset_topology = trajectories[0].topology

        # Little monkey-patch to prevent further subsetting Topologies
        # this modified version of the topology will never exit this function
        kwargs['top'].subset = lambda atom_indices : subset_topology



    # We know the topology is equal because we send the same topology
    # kwarg in. Therefore, we explictly throw away the topology on all
    # but the first trajectory by making them all point to None
    #  and use check_topology=False on the join.
    # Throwing the topology away explictly allows a large number of pdb
    # files to be read in without using ridiculous amounts of memory.
    for f in filename_or_filenames:
        t = loader(f, **kwargs)

        t.topology = None
        trajectories.append(t)


    if len(trajectories) == 1: #if only one file was given there is nothing to join
        return trajectories[0]

    return join(trajectories, check_topology=False,
                discard_overlapping_frames=discard_overlapping_frames)


def iterload(filename, chunk=100, **kwargs):
    """An iterator over a trajectory from one or more files on disk, in fragments

    This may be more memory efficient than loading an entire trajectory at
    once

    Parameters
    ----------
    filename : path-like
        Path to the trajectory file on disk
    chunk : int
        Number of frames to load at once from disk per iteration.  If 0, load all.

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
    skip : int, default=0
        Skip first n frames.

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
    stride = kwargs.pop('stride', 1)
    atom_indices = cast_indices(kwargs.pop('atom_indices', None))
    top = kwargs.pop('top', None)
    skip = kwargs.pop('skip', 0)

    extension = _get_extension(filename)
    if extension not in _TOPOLOGY_EXTS:
        topology = _parse_topology(top)

    if chunk == 0:
        # If chunk was 0 then we want to avoid filetype-specific code
        # in case of undefined behavior in various file parsers.
        # TODO: this will first apply stride, then skip!
        if extension not in _TOPOLOGY_EXTS:
            kwargs['top'] = top
        yield load(filename, **kwargs)[skip:]
    elif extension in ('.pdb', '.pdb.gz'):
        # the PDBTrajectortFile class doesn't follow the standard API. Fixing it
        # to support iterload could be worthwhile, but requires a deep refactor.
        t = load(filename, stride=stride, atom_indices=atom_indices)
        for i in range(0, len(t), chunk):
            yield t[i:i+chunk]
    elif extension in ('.gsd'):
        i = 0
        while True:
            traj = load(filename, stride=stride, atom_indices=atom_indices,
                    start=i, n_frames=chunk)
            if len(traj) ==0 :
                return
            i += chunk
            yield traj
    else:
        with (lambda x: open(x, n_atoms=topology.n_atoms)
              if extension in ('.crd', '.mdcrd')
              else open(filename))(filename) as f:
            if skip > 0:
                f.seek(skip)
            while True:
                if extension not in _TOPOLOGY_EXTS:
                    traj = f.read_as_traj(topology, n_frames=chunk, stride=stride, atom_indices=atom_indices, **kwargs)
                else:
                    traj = f.read_as_traj(n_frames=chunk, stride=stride, atom_indices=atom_indices, **kwargs)

                if len(traj) == 0:
                    return

                yield traj

def join(trajs, check_topology=True, discard_overlapping_frames=False):
    """Concatenate multiple trajectories into one long trajectory

    Parameters
    ----------
    trajs : iterable of trajectories
        Combine these into one trajectory
    check_topology : bool
        Make sure topologies match before joining
    discard_overlapping_frames : bool
        Check for overlapping frames and discard
    """
    return functools.reduce(
        lambda x, y:
        x.join(y, check_topology=check_topology,
               discard_overlapping_frames=discard_overlapping_frames),
        trajs
    )

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
    def n_chains(self):
        """Number of chains in the trajectory

        Returns
        -------
        n_chains : int
            The number of chains in the trajectory's topology
        """
        if self.top is None:
            return 0
        return sum([1 for c in self.top.chains])

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

    @top.setter
    def top(self, value):
        "Set the topology of the system, describing the organization of atoms into residues, bonds, etc"
        # todo: more typechecking
        self._topology = value

    @property
    def timestep(self):
        """Timestep between frames, in picoseconds

        Returns
        -------
        timestep : float
            The timestep between frames, in picoseconds.
        """
        if self.n_frames <= 1:
            raise(ValueError("Cannot calculate timestep if trajectory has one frame."))
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
            Vectors defining the shape of the unit cell in each frame.
            The semantics of this array are that the shape of the unit cell
            in frame ``i`` are given by the three vectors, ``value[i, 0, :]``,
            ``value[i, 1, :]``, and ``value[i, 2, :]``.
        """
        if self._unitcell_lengths is None or self._unitcell_angles is None:
            return None

        v1, v2, v3 = lengths_and_angles_to_box_vectors(
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
        if vectors is None or np.all(np.abs(vectors) < 1e-15):
            self._unitcell_lengths = None
            self._unitcell_angles = None
            return

        if not len(vectors) == len(self):
            raise TypeError('unitcell_vectors must be the same length as '
                            'the trajectory. you provided %s' % str(vectors))

        v1 = vectors[:, 0, :]
        v2 = vectors[:, 1, :]
        v3 = vectors[:, 2, :]
        a, b, c, alpha, beta, gamma = box_vectors_to_lengths_and_angles(v1, v2, v3)

        self._unitcell_lengths = np.vstack((a, b, c)).T
        self._unitcell_angles =  np.vstack((alpha, beta, gamma)).T


    @property
    def unitcell_volumes(self):
        """Volumes of unit cell for each frame.

        Returns
        -------
        volumes : {np.ndarray, shape=(n_frames), None}
            Volumes of the unit cell in each frame, in nanometers^3, or None
            if the Trajectory contains no unitcell information.
        """
        if self.unitcell_lengths is not None:
            return np.array(list(map(np.linalg.det, self.unitcell_vectors)))
        else:
            return None

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

    def _string_summary_basic(self):
        """Basic summary of traj in string form."""
        unitcell_str = 'and unitcells' if self._have_unitcell else 'without unitcells'
        value = "mdtraj.Trajectory with %d frames, %d atoms, %d residues, %s" % (
                    self.n_frames, self.n_atoms, self.n_residues, unitcell_str)
        return value

    def __len__(self):
        return self.n_frames

    def __add__(self, other):
        "Concatenate two trajectories"
        return self.join(other)

    def __str__(self):
        return "<%s>" % (self._string_summary_basic())

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __hash__(self):
        hash_value = hash(self.top)
        # combine with hashes of arrays
        hash_value ^= _hash_numpy_array(self._xyz)
        hash_value ^= _hash_numpy_array(self.time)
        if self._unitcell_lengths is not None:
            hash_value ^= _hash_numpy_array(self._unitcell_lengths)
        if self._unitcell_angles is not None:
            hash_value ^= _hash_numpy_array(self._unitcell_angles)
        return hash_value

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    # def describe(self):
    #     """Diagnostic summary statistics on the trajectory"""
    #     # What information do we want to display?
    #     # Goals: easy to figure out if a trajectory is blowing up or contains
    #     # bad data, easy to diagonose other problems. Generally give a
    #     # high-level description of the data in the trajectory.
    #     # Possibly show std. dev. of differnt coordinates in the trajectory
    #     # or maybe its RMSD drift or something?
    #     # Also, check for any NaNs or Infs in the data. Or other common issues
    #     # like that?
    #     # Note that pandas.DataFrame has a describe() method, which gives
    #     # min/max/mean/std.dev./percentiles of each column in a DataFrame.
    #     raise NotImplementedError()

    def superpose(self, reference, frame=0, atom_indices=None,
                  ref_atom_indices=None, parallel=True):
        """Superpose each conformation in this trajectory upon a reference

        Parameters
        ----------
        reference : md.Trajectory
            Align self to a particular frame in `reference`
        frame : int
            The index of the conformation in `reference` to align to.
        atom_indices : array_like, or None
            The indices of the atoms to superpose. If not
            supplied, all atoms will be used.
        ref_atom_indices : array_like, or None
            Use these atoms on the reference structure. If not supplied,
            the same atom indices will be used for this trajectory and the
            reference one.
        parallel : bool
            Use OpenMP to run the superposition in parallel over multiple cores

        Returns
        -------
        self
        """

        if atom_indices is None:
            atom_indices = slice(None)

        if ref_atom_indices is None:
            ref_atom_indices = atom_indices

        if not isinstance(ref_atom_indices, slice) and (
            len(ref_atom_indices) != len(atom_indices)):
            raise ValueError("Number of atoms must be consistent!")

        n_frames = self.xyz.shape[0]
        self_align_xyz = np.asarray(self.xyz[:, atom_indices, :], order='c')
        self_displace_xyz = np.asarray(self.xyz, order='c')
        ref_align_xyz = np.array(reference.xyz[frame, ref_atom_indices, :],
                                 copy=True, order='c').reshape(1, -1, 3)

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

        self_displace_xyz += ref_offset
        self.xyz = self_displace_xyz
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

    def stack(self, other, keep_resSeq=True):
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
        keep_resSeq : bool, optional, default=True
            see ```mdtraj.core.topology.Topology.join``` method documentation

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
            topology = self.topology.join(other.topology, keep_resSeq=keep_resSeq)
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
        rmsd_traces = None
        if self.unitcell_angles is not None:
            unitcell_angles = self.unitcell_angles[key]
        if self.unitcell_lengths is not None:
            unitcell_lengths = self.unitcell_lengths[key]
        if self._rmsd_traces is not None:
            rmsd_traces = self._rmsd_traces

        if copy:
            xyz = xyz.copy()
            time = time.copy()
            topology = deepcopy(self._topology)

            if self.unitcell_angles is not None:
                unitcell_angles = unitcell_angles.copy()
            if self.unitcell_lengths is not None:
                unitcell_lengths = unitcell_lengths.copy()
            if rmsd_traces is not None :
                rmsd_traces = rmsd_traces.copy()
        else:
            topology = self._topology

        newtraj = self.__class__(
            xyz, topology, time, unitcell_lengths=unitcell_lengths,
            unitcell_angles=unitcell_angles)

        if rmsd_traces is not None:
            newtraj._rmsd_traces = rmsd_traces

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
        self._time_default_to_arange = (time is None)
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
        from openmm import Vec3
        from openmm.unit import nanometer

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
        from openmm import Vec3
        from openmm.unit import nanometer

        vectors = self.unitcell_vectors[frame]
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
        filenames : {path-like, [path-like]}
            Either a path or list of paths

        Other Parameters
        ----------------
        As requested by the various load functions -- it depends on the extension
        """
        return load(filenames, **kwargs)

    def _savers(self):
        """Return a dictionary mapping extensions to the appropriate format-specific save function"""
        return {'.xtc': self.save_xtc,
                '.trr': self.save_trr,
                '.pdb': self.save_pdb,
                '.pdb.gz': self.save_pdb,
                '.dcd': self.save_dcd,
                '.h5': self.save_hdf5,
                '.binpos': self.save_binpos,
                '.nc': self.save_netcdf,
                '.netcdf': self.save_netcdf,
                '.ncrst' : self.save_netcdfrst,
                '.crd': self.save_mdcrd,
                '.mdcrd': self.save_mdcrd,
                '.ncdf': self.save_netcdf,
                '.lh5': self.save_lh5,
                '.lammpstrj': self.save_lammpstrj,
                '.xyz': self.save_xyz,
                '.xyz.gz': self.save_xyz,
                '.gro': self.save_gro,
                '.rst7' : self.save_amberrst7,
                '.tng' : self.save_tng,
                '.dtr': self.save_dtr,
                '.gsd': self.save_gsd,
            }

    def save(self, filename, **kwargs):
        """Save trajectory to disk, in a format determined by the filename extension

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory. The extension will
            be parsed and will control the format.

        Other Parameters
        ----------------
        lossy : bool
            For .h5 or .lh5, whether or not to use compression.
        no_models: bool
            For .pdb. TODO: Document this?
        force_overwrite : bool
            If `filename` already exists, overwrite it.
        """
        # grab the extension of the filename
        extension = _get_extension(filename)
        savers = self._savers()

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
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with HDF5TrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(coordinates=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    time=self.time,
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)
            f.topology = self.topology

    def save_lammpstrj(self, filename, force_overwrite=True):
        """Save trajectory to LAMMPS custom dump format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with LAMMPSTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)

    def save_xyz(self, filename, force_overwrite=True):
        """Save trajectory to .xyz format.

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with XYZTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    types=[a.name for a in self.top.atoms])

    def save_pdb(self, filename, force_overwrite=True, bfactors=None):
        """Save trajectory to RCSB PDB format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        bfactors : array_like, default=None, shape=(n_frames, n_atoms) or (n_atoms,)
            Save bfactors with pdb file. If the array is two dimensional it should
            contain a bfactor for each atom in each frame of the trajectory.
            Otherwise, the same bfactor will be saved in each frame.
        """
        self._check_valid_unitcell()

        if not bfactors is None:
            if len(np.array(bfactors).shape) == 1:
                if len(bfactors) != self.n_atoms:
                    raise ValueError("bfactors %s should be shaped as (n_frames, n_atoms) or (n_atoms,)" % str(np.array(bfactors).shape))

                bfactors = [bfactors] * self.n_frames

            else:
                if np.array(bfactors).shape != (self.n_frames, self.n_atoms):
                    raise ValueError("bfactors %s should be shaped as (n_frames, n_atoms) or (n_atoms,)" % str(np.array(bfactors).shape))

        else:
            bfactors = [None] * self.n_frames


        with PDBTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            for i in xrange(self.n_frames):

                if self._have_unitcell:
                    f.write(in_units_of(self._xyz[i], Trajectory._distance_unit, f.distance_unit),
                            self.topology,
                            modelIndex=i,
                            bfactors=bfactors[i],
                            unitcell_lengths=in_units_of(self.unitcell_lengths[i], Trajectory._distance_unit, f.distance_unit),
                            unitcell_angles=self.unitcell_angles[i])
                else:
                    f.write(in_units_of(self._xyz[i], Trajectory._distance_unit, f.distance_unit),
                            self.topology,
                            modelIndex=i,
                            bfactors=bfactors[i])

    def save_xtc(self, filename, force_overwrite=True):
        """Save trajectory to Gromacs XTC format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with XTCTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    time=self.time,
                    box=in_units_of(self.unitcell_vectors, Trajectory._distance_unit, f.distance_unit))

    def save_trr(self, filename, force_overwrite=True):
        """Save trajectory to Gromacs TRR format

        Notes
        -----
        Only the xyz coordinates and the time are saved, the velocities
        and forces in the trr will be zeros

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with TRRTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    time=self.time,
                    box=in_units_of(self.unitcell_vectors, Trajectory._distance_unit, f.distance_unit))

    def save_dcd(self, filename, force_overwrite=True):
        """Save trajectory to CHARMM/NAMD DCD format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filenames, if its already there
        """
        self._check_valid_unitcell()
        with DCDTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)

    def save_dtr(self, filename, force_overwrite=True):
        """Save trajectory to DESMOND DTR format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filenames, if its already there
        """
        self._check_valid_unitcell()
        with DTRTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles,
                    times=self.time)

    def save_binpos(self, filename, force_overwrite=True):
        """Save trajectory to AMBER BINPOS format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        with BINPOSTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit))


    def save_mdcrd(self, filename, force_overwrite=True):
        """Save trajectory to AMBER mdcrd format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        self._check_valid_unitcell()
        if self._have_unitcell:
            if not np.all(self.unitcell_angles == 90):
                raise ValueError('Only rectilinear boxes can be saved to mdcrd files. '
                                 'Your angles are {}'.format(self.unitcell_angles))

        with MDCRDTrajectoryFile(filename, mode='w', force_overwrite=force_overwrite) as f:
            f.write(xyz=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit))


    def save_netcdf(self, filename, force_overwrite=True):
        """Save trajectory in AMBER NetCDF format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if it's already there
        """
        self._check_valid_unitcell()
        with NetCDFTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(coordinates=in_units_of(self._xyz, Trajectory._distance_unit, NetCDFTrajectoryFile.distance_unit),
                    time=self.time,
                    cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                    cell_angles=self.unitcell_angles)

    def save_netcdfrst(self, filename, force_overwrite=True):
        """Save trajectory in AMBER NetCDF restart format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the restart
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if it's already there

        Notes
        -----
        NetCDF restart files can only store a single frame. If only one frame
        exists, "filename" will be written.  Otherwise, "filename.#" will be
        written, where # is a zero-padded number from 1 to the total number of
        frames in the trajectory
        """
        self._check_valid_unitcell()
        if self.n_frames == 1:
            with AmberNetCDFRestartFile(filename, 'w', force_overwrite=force_overwrite) as f:
                coordinates = in_units_of(self._xyz, Trajectory._distance_unit,
                                          AmberNetCDFRestartFile.distance_unit)
                lengths = in_units_of(self.unitcell_lengths, Trajectory._distance_unit,
                                      AmberNetCDFRestartFile.distance_unit)
                f.write(coordinates=coordinates, time=self.time[0],
                        cell_lengths=lengths, cell_angles=self.unitcell_angles)
        else:
            fmt = '%s.%%0%dd' % (filename, len(str(self.n_frames)))
            for i in xrange(self.n_frames):
                with AmberNetCDFRestartFile(fmt % (i+1), 'w', force_overwrite=force_overwrite) as f:
                    coordinates = in_units_of(self._xyz, Trajectory._distance_unit,
                                              AmberNetCDFRestartFile.distance_unit)
                    lengths = in_units_of(self.unitcell_lengths, Trajectory._distance_unit,
                                          AmberNetCDFRestartFile.distance_unit)
                    f.write(coordinates=coordinates[i], time=self.time[i],
                            cell_lengths=lengths[i], cell_angles=self.unitcell_angles[i])

    def save_amberrst7(self, filename, force_overwrite=True):
        """Save trajectory in AMBER ASCII restart format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the restart
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if it's already there

        Notes
        -----
        Amber restart files can only store a single frame. If only one frame
        exists, "filename" will be written.  Otherwise, "filename.#" will be
        written, where # is a zero-padded number from 1 to the total number of
        frames in the trajectory
        """
        self._check_valid_unitcell()
        if self.n_frames == 1:
            with AmberRestartFile(filename, 'w', force_overwrite=force_overwrite) as f:
                coordinates = in_units_of(self._xyz, Trajectory._distance_unit,
                                          AmberRestartFile.distance_unit)
                lengths = in_units_of(self.unitcell_lengths, Trajectory._distance_unit,
                                      AmberRestartFile.distance_unit)
                f.write(coordinates=coordinates, time=self.time[0],
                        cell_lengths=lengths, cell_angles=self.unitcell_angles)
        else:
            fmt = '%s.%%0%dd' % (filename, len(str(self.n_frames)))
            for i in xrange(self.n_frames):
                with AmberRestartFile(fmt % (i+1), 'w', force_overwrite=force_overwrite) as f:
                    coordinates = in_units_of(self._xyz, Trajectory._distance_unit,
                                              AmberRestartFile.distance_unit)
                    lengths = in_units_of(self.unitcell_lengths, Trajectory._distance_unit,
                                          AmberRestartFile.distance_unit)
                    f.write(coordinates=coordinates[i], time=self.time[0],
                            cell_lengths=lengths[i], cell_angles=self.unitcell_angles[i])

    def save_lh5(self, filename, force_overwrite=True):
        """Save trajectory in deprecated MSMBuilder2 LH5 (lossy HDF5) format.

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if it's already there
        """
        with LH5TrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(coordinates=self.xyz)
            f.topology = self.topology

    def save_gro(self, filename, force_overwrite=True, precision=3):
        """Save trajectory in Gromacs .gro format

        Parameters
        ----------
        filename : path-like
            Path to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at that filename if it exists
        precision : int, default=3
            The number of decimal places to use for coordinates in GRO file
        """
        self._check_valid_unitcell()
        with GroTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(self.xyz, self.topology, self.time, self.unitcell_vectors,
                    precision=precision)

    def save_tng(self, filename, force_overwrite=True):
        """Save trajectory to Gromacs TNG format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        self._check_valid_unitcell()
        with TNGTrajectoryFile(os.fspath(filename), 'w', force_overwrite=force_overwrite) as f:
            f.write(self.xyz, time=self.time, box=self.unitcell_vectors)

    def save_gsd(self, filename, force_overwrite=True):
        """Save trajectory to HOOMD GSD format

        Parameters
        ----------
        filename : path-like
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filenames, if its already there
        """
        if os.path.exists(filename) and not force_overwrite:
            raise IOError('"%s" already exists' % filename)

        self._check_valid_unitcell()
        write_gsd(filename, self.xyz, self.topology,
                cell_lengths=self.unitcell_lengths,
                cell_angles=self.unitcell_angles)

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
            self.xyz -= distance.compute_center_of_mass(self)[:, np.newaxis, :]
        else:
            self._rmsd_traces = _rmsd._center_inplace_atom_major(self._xyz)

        return self

    @deprecated('restrict_atoms was replaced by atom_slice and will be removed in 2.0')
    def restrict_atoms(self, atom_indices, inplace=True):
        """Retain only a subset of the atoms in a trajectory

        Deletes atoms not in `atom_indices`, and re-indexes those that remain

        Parameters
        ----------
        atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of atom indices to keep.
        inplace : bool, default=True
            If ``True``, the operation is done inplace, modifying ``self``.
            Otherwise, a copy is returned with the restricted atoms, and
            ``self`` is not modified.

        Returns
        -------
        traj : md.Trajectory
            The return value is either ``self``, or the new trajectory,
            depending on the value of ``inplace``.
        """
        return self.atom_slice(atom_indices, inplace=inplace)

    def atom_slice(self, atom_indices, inplace=False):
        """Create a new trajectory from a subset of atoms

        Parameters
        ----------
        atom_indices : array-like, dtype=int, shape=(n_atoms)
            List of indices of atoms to retain in the new trajectory.
        inplace : bool, default=False
            If ``True``, the operation is done inplace, modifying ``self``.
            Otherwise, a copy is returned with the sliced atoms, and
            ``self`` is not modified.

        Returns
        -------
        traj : md.Trajectory
            The return value is either ``self``, or the new trajectory,
            depending on the value of ``inplace``.

        See Also
        --------
        stack : stack multiple trajectories along the atom axis
        """
        xyz = np.array(self.xyz[:, atom_indices], order='C')
        topology = None
        if self._topology is not None:
            topology = self._topology.subset(atom_indices)

        if inplace:
            if self._topology is not None:
                self._topology = topology
            self._xyz = xyz

            return self

        unitcell_lengths = unitcell_angles = None
        if self._have_unitcell:
            unitcell_lengths = self._unitcell_lengths.copy()
            unitcell_angles = self._unitcell_angles.copy()
        time = self._time.copy()

        return Trajectory(xyz=xyz, topology=topology, time=time,
                          unitcell_lengths=unitcell_lengths,
                          unitcell_angles=unitcell_angles)

    def remove_solvent(self, exclude=None, inplace=False):
        """
        Create a new trajectory without solvent atoms

        Parameters
        ----------
        exclude : array-like, dtype=str, shape=(n_solvent_types)
            List of solvent residue names to retain in the new trajectory.
        inplace : bool, default=False
            The return value is either ``self``, or the new trajectory,
            depending on the value of ``inplace``.

        Returns
        -------
        traj : md.Trajectory
            The return value is either ``self``, or the new trajectory,
            depending on the value of ``inplace``.
        """
        solvent_types = list(_SOLVENT_TYPES)

        if exclude is not None:

            if isinstance(exclude, str):
                raise TypeError('exclude must be array-like')
            if not isinstance(exclude, Iterable):
                raise TypeError('exclude is not iterable')

            for type in exclude:
                if type not in solvent_types:
                    raise ValueError(type + 'is not a valid solvent type')
                solvent_types.remove(type)

        atom_indices = [atom.index for atom in self.topology.atoms if
                atom.residue.name not in solvent_types]

        return self.atom_slice(atom_indices, inplace = inplace)

    def smooth(self, width, order=3, atom_indices=None, inplace=False):
        """Smoothen a trajectory using a zero-delay Buttersworth filter. Please
        note that for optimal results the trajectory should be properly aligned
        prior to smoothing (see `md.Trajectory.superpose`).

        Parameters
        ----------
        width : int
            This acts very similar to the window size in a moving average
            smoother. In this implementation, the frequency of the low-pass
            filter is taken to be two over this width, so it's like
            "half the period" of the sinusiod where the filter starts
            to kick in. Must be an integer greater than one.
        order : int, optional, default=3
            The order of the filter. A small odd number is recommended. Higher
            order filters cutoff more quickly, but have worse numerical
            properties.
        atom_indices : array-like, dtype=int, shape=(n_atoms), default=None
            List of indices of atoms to retain in the new trajectory.
            Default is set to `None`, which applies smoothing to all atoms.
        inplace : bool, default=False
            The return value is either ``self``, or the new trajectory,
            depending on the value of ``inplace``.

        Returns
        -------
        traj : md.Trajectory
            The return value is either ``self``, or the new smoothed trajectory,
            depending on the value of ``inplace``.

        References
        ----------
        .. [1] "FiltFilt". Scipy Cookbook. SciPy. <http://www.scipy.org/Cookbook/FiltFilt>.
        """
        from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

        if width < 2.0 or not isinstance(width, int):
            raise ValueError('width must be an integer greater than 1.')
        if not atom_indices:
            atom_indices = range(self.n_atoms)

        # find nearest odd integer
        pad = int(np.ceil((width + 1)/2)*2 - 1)

        # Use lfilter_zi to choose the initial condition of the filter.
        b, a = butter(order, 2.0 / width)
        zi = lfilter_zi(b, a)

        xyz = self.xyz.copy()

        for i in atom_indices:
            for j in range(3):

                signal = xyz[:, i, j]
                padded = np.r_[signal[pad - 1: 0: -1], signal, signal[-1: -pad: -1]]

                # Apply the filter to the width.
                z, _ = lfilter(b, a, padded, zi=zi*padded[0])

                # Apply the filter again, to have a result filtered at an order
                # the same as filtfilt.
                z2, _ = lfilter(b, a, z, zi=zi*z[0])

                # Use filtfilt to apply the filter.
                output = filtfilt(b, a, padded)

                xyz[:, i, j] = output[(pad-1): -(pad-1)]


        if not inplace:
            return Trajectory(xyz=xyz, topology=self.topology,
                              time=self.time,
                              unitcell_lengths=self.unitcell_lengths,
                              unitcell_angles=self.unitcell_angles)

        self.xyz = xyz

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

    def make_molecules_whole(self, inplace=False, sorted_bonds=None):
        """Only make molecules whole

        Parameters
        ----------
        inplace : bool
            If False, a new Trajectory is created and returned.
            If True, this Trajectory is modified directly.
        sorted_bonds : array of shape (n_bonds, 2)
            Pairs of atom indices that define bonds, in sorted order.
            If not specified, these will be determined from the trajectory's
            topology.

        See Also
        --------
        image_molecules
        """
        unitcell_vectors = self.unitcell_vectors
        if unitcell_vectors is None:
            raise ValueError('This Trajectory does not define a periodic unit cell')

        if inplace:
            result = self
        else:
            result = Trajectory(xyz=self.xyz, topology=self.topology,
                                time=self.time,
                                unitcell_lengths=self.unitcell_lengths,
                                unitcell_angles=self.unitcell_angles)

        if sorted_bonds is None:
            sorted_bonds = sorted(self._topology.bonds, key=lambda bond: bond[0].index)
            sorted_bonds = np.asarray([[b0.index, b1.index] for b0, b1 in sorted_bonds], dtype=np.int32)

        box = np.asarray(result.unitcell_vectors, order='c')
        _geometry.whole_molecules(result.xyz, box, sorted_bonds)
        if not inplace:
            return result
        return self

    def image_molecules(self, inplace=False, anchor_molecules=None, other_molecules=None, sorted_bonds=None, make_whole=True):
        """Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory.

        This method is useful for visualizing a trajectory in which molecules were not wrapped
        to the periodic unit cell, or in which the macromolecules are not centered with respect
        to the solvent.  It tries to be intelligent in deciding what molecules to center, so you
        can simply call it and trust that it will "do the right thing".

        Parameters
        ----------
        inplace : bool, default=False
            If False, a new Trajectory is created and returned.  If True, this Trajectory
            is modified directly.
        anchor_molecules : list of atom sets, optional, default=None
            Molecule that should be treated as an "anchor".
            These molecules will be centered in the box and put near each other.
            If not specified, anchor molecules are guessed using a heuristic.
        other_molecules : list of atom sets, optional, default=None
            Molecules that are not anchors. If not specified,
            these will be molecules other than the anchor molecules
        sorted_bonds : array of shape (n_bonds, 2)
            Pairs of atom indices that define bonds, in sorted order.
            If not specified, these will be determined from the trajectory's
            topology. Only relevant if ``make_whole`` is True.
        make_whole : bool
            Whether to make molecules whole.

        Returns
        -------
        traj : md.Trajectory
            The return value is either ``self`` or the new trajectory,
            depending on the value of ``inplace``.

        See Also
        --------
        Topology.guess_anchor_molecules
        """
        unitcell_vectors = self.unitcell_vectors
        if unitcell_vectors is None:
            raise ValueError('This Trajectory does not define a periodic unit cell')

        if anchor_molecules is None:
            anchor_molecules = self.topology.guess_anchor_molecules()

        if other_molecules is None:
            # Determine other molecules by which molecules are not anchor molecules
            molecules = self._topology.find_molecules()
            other_molecules = [mol for mol in molecules if mol not in anchor_molecules]

        # Expand molecules into atom indices
        anchor_molecules_atom_indices = [np.fromiter((a.index for a in mol), dtype=np.int32) for mol in anchor_molecules]
        other_molecules_atom_indices  = [np.fromiter((a.index for a in mol), dtype=np.int32) for mol in other_molecules]

        if inplace:
            result = self
        else:
            result = Trajectory(xyz=self.xyz, topology=self.topology, time=self.time,
                unitcell_lengths=self.unitcell_lengths, unitcell_angles=self.unitcell_angles)

        if make_whole and sorted_bonds is None:
            sorted_bonds = sorted(self._topology.bonds, key=lambda bond: bond[0].index)
            sorted_bonds = np.asarray([[b0.index, b1.index] for b0, b1 in sorted_bonds], dtype=np.int32)
        elif not make_whole:
            sorted_bonds = None

        box = np.asarray(result.unitcell_vectors, order='c')
        _geometry.image_molecules(result.xyz, box, anchor_molecules_atom_indices, other_molecules_atom_indices, sorted_bonds)
        if not inplace:
            return result
        return self
