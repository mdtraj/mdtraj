# Copyright 2012 mdtraj developers
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

##############################################################################
# Imports
##############################################################################

import os
import warnings
import logging
import functools
from itertools import izip
from copy import deepcopy
import numpy as np

from mdtraj import (DCDTrajectoryFile, BINPOSTrajectoryFile, XTCTrajectoryFile,
                    TRRTrajectoryFile, HDF5TrajectoryFile, NetCDFTrajectoryFile,
                    PDBTrajectoryFile, Topology)
from mdtraj.utils import unitcell, ensure_type

try:
    from simtk.openmm import Vec3
    from simtk.unit import nanometer
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


__all__ = ['Trajectory', 'load', 'load_pdb', 'load_xtc', 'load_trr', 'load_binpos',
           'load_dcd', 'load_netcdf', 'load_hdf5', 'load_netcdf', 'load_xml']

##############################################################################
# Globals
##############################################################################

logger = logging.getLogger(__name__)

# note, there's another global named "_LoaderRegistry" that's declared at
# the bottom of the file

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
    if isinstance(filenames, basestring):
        filenames = [filenames]
    for fn in filenames:
        if not os.path.exists(fn):
            raise IOError("I'm sorry, the file you requested does not seem to "
            "exist: %s" % fn)


def _parse_topology(top):
    """Get the topology from a argument of indeterminate type
    If top is a string, we try loading a pdb, if its a trajectory
    we extract its topology.
    """
    if isinstance(top, basestring):
        topology = PDBTrajectoryFile(top).topology
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, Topology):
        topology = top
    else:
        raise TypeError('Could not interpreted top=%s' % top)

    return topology


def load(filename_or_filenames, discard_overlapping_frames=False, **kwargs):
    """Load a trajectory from one or more files on disk.
    
    This function dispatches to one of the specialized trajectory loaders based
    on the extension on the filename. Because different trajectory formats save
    different information on disk, the specific keyword argument options supported
    depend on the specific loaded.
    
    Parameters
    ----------
    filename_or_filenames : {str, list of strings}
        filename or list of filenames containing trajectory files of a single format.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them


    Other Parameters
    ----------------
    DEPENDS ON THE LOADER. NOT FILLED IN YET. NEED TO ADD STRIDING AND SUCH

    See Also
    --------
    load_pdb, load_xtc, load_trr, load_hdf5, load_netcdf, load_dcd, load_binpos

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """

    _assert_files_exist(filename_or_filenames)

    # grab the extension of the filename
    if isinstance(filename_or_filenames, basestring):  # If a single filename
        extension = os.path.splitext(filename_or_filenames)[1]
        filename = filename_or_filenames
    else:  # If multiple filenames, take the first one.
        extensions = [os.path.splitext(filename_i)[1] for filename_i in filename_or_filenames]
        if len(set(extensions)) != 1:
            raise(TypeError("All filenames must have same extension!"))
        else:
            return functools.reduce(lambda a, b: a.join(b, discard_overlapping_frames=discard_overlapping_frames), (load(f,**kwargs) for f in filename_or_filenames))

    try:
        loader = _LoaderRegistry[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files '
                      'with extensions in %s' % (filename, extension, _LoaderRegistry.keys()))

    return loader(filename, **kwargs)


def load_pdb(filename):
    """Load a RCSB Protein Data Bank file from disk.

    Parameters
    ----------
    filename : str
        Path to the PDB file on disk.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
        
    See Also
    --------
    mdtraj.PDBTrajectoryFile : Low level interface to PDB files
    """
    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_pdb. '
            'you supplied %s' % type(filename))

    filename = str(filename)
    f = PDBTrajectoryFile(filename)

    # convert from angstroms to nm
    coords = f.positions / 10.0

    assert coords.ndim == 3, 'internal shape error'
    n_frames = len(coords)

    trajectory = Trajectory(xyz=coords, topology=f.topology)

    if f.unitcell_vectors is not None:
        a, b, c = f.unitcell_vectors
        # we need to convert the distances from angstroms to nanometers
        unitcell_lengths = np.array([[a / 10.0, b / 10.0, c / 10.0]])
        unitcell_angles = np.array([[90.0, 90.0, 90.0]])

        # we need to project these unitcell parameters
        # into each frame, for a multiframe PDB

        trajectory.unitcell_lengths = np.repeat(unitcell_lengths, n_frames, axis=0)
        trajectory.unitcell_angles = np.repeat(unitcell_angles, n_frames, axis=0)

    return trajectory


def load_xml(filename, top=None):
    """Load a single conformation from an XML file, such as those
    produced by OpenMM

    Note: The OpenMM serialized state XML format contains information that
    is not read by this method, including forces, energies, and velocities.
    Here, we just read the positions and the box vectors.

    Parameters
    ----------
    filename : string
        The path on disk to the XML file
    top : {str, Trajectory, Topology}
        The XML format does not contain topology information. Pass in either the
        path to a pdb file, a trajectory, or a topology to supply this information.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """
    topology = _parse_topology(top)

    import xml.etree.cElementTree as etree
    tree = etree.parse(filename)

    # get all of the positions from the XML into a list of tuples
    # then convert to a numpy array
    positions = []
    for position in tree.getroot().find('Positions'):
        positions.append((float(position.attrib['x']),
                          float(position.attrib['y']),
                          float(position.attrib['z'])))

    box = []
    vectors = tree.getroot().find('PeriodicBoxVectors')
    for name in ['A', 'B', 'C']:
        box.append((float(vectors.find(name).attrib['x']),
                    float(vectors.find(name).attrib['y']),
                    float(vectors.find(name).attrib['z'])))

    traj = Trajectory(xyz=np.array(positions), topology=topology)
    traj.unitcell_vectors = np.array(box).reshape(1,3,3)

    return traj


def load_xtc(filename, top=None, chunk=None):
    """Load an Gromacs XTC file from disk.
    
    Since the Gromacs XTC format doesn't contain information to specify the
    topolgy, you need to supply the topology yourself.

    Parameters
    ----------
    filename : str
        Filename (string) of xtc trajectory.
    top : {str, Trajectory, Topology}
        The XTC format does not contain topology information. Pass in either the
        path to a RCSB PDB file, a trajectory, or a topology to supply this
        information.
    chunk : None
        This option is depricated.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.XTCTrajectoryFile :  Low level interface to XTC files
    """
    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_xtc')

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_xtc. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)

    with XTCTrajectoryFile(filename, 'r') as f:
        xyz, time, step, box = f.read()

    # note we're not doing anything with the box vectors
    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    trajectory.unitcell_vectors = box

    return trajectory


def load_trr(filename, top=None, chunk=None):
    """Load a trr file. Since the trr doesn't contain information
    to specify the topolgy, you need to supply the topology yourself

    Parameters
    ----------
    filename : str
        Filename of TRR trajectory file.
    top : {str, Trajectory, Topology}
        The TRR format does not contain topology information. Pass in either the
        path to a pdb file, a trajectory, or a topology to supply this information.
    chunk : None
        This option is depricated.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.TRRTrajectoryFile :  Low level interface to TRR files
    """
    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_trr')

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_trr. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)

    with TRRTrajectoryFile(filename) as f:
        xyz, time, step, box, lambd = f.read()

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    trajectory.unitcell_vectors = box
    return trajectory


def load_dcd(filename, top=None):
    """Load an xtc file. Since the dcd format doesn't contain information
    to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filename : str
        String filename of DCD file.
    top : {str, Trajectoray, Topology}
        DCD XTC format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.DCDTrajectoryFile :  Low level interface to DCD files
    """
    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_dcd')

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_trr. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    with DCDTrajectoryFile(filename) as f:
        xyz, box_length, box_angle = f.read()

    xyz /= 10.  # convert from anstroms to nanometer
    box_length /= 10.  # convert from anstroms to nanometer

    trajectory = Trajectory(xyz=xyz, topology=topology)
    trajectory.unitcell_lengths = box_length
    trajectory.unitcell_angles = box_angle

    return trajectory


def load_hdf5(filename, stride=None, frame=None):
    """Load an MDTraj hdf5 trajectory file from disk.

    Parameters
    ----------
    filename : str
        String filename of HDF Trajectory file.
    stride : int, default=None
        Only read every stride-th frame
    frame : {None, int}, default=None
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.HDF5TrajectoryFile :  Low level interface to HDF5 files
    """
    tf = HDF5TrajectoryFile(filename)
    if frame is None:
        data = tf.read(stride=stride)
    else:
        tf._frame_index += int(frame)
        data = tf.read(n_frames=1)

    trajectory = Trajectory(xyz=data.coordinates, topology=tf.topology, time=data.time)
    trajectory.unitcell_lengths = data.cell_lengths
    trajectory.unitcell_angles = data.cell_angles

    return trajectory


def load_binpos(filename, top=None, chunk=500):
    """Load an AMBER BINPOS file.

    Parameters
    ----------
    filename : str
        String filename of AMBER binpos file.
    top : {str, Trajectory, Topology}
        The BINPOS format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.BINPOSTrajectoryFile :  Low level interface to BINPOS files
    """
    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_binpos')

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_binpos. '
            'you supplied %s' % type(filename))


    topology = _parse_topology(top)

    with BINPOSTrajectoryFile(filename) as f:
        xyz = f.read()
    xyz /= 10.0  # convert from anstroms to nanometer

    return Trajectory(xyz=xyz, topology=topology)


def load_netcdf(filename, top=None, stride=None):
    """Load an AMBER NetCDF file. Since the NetCDF format doesn't contain
    information to specify the topolgy, you need to supply a topology

    Parameters
    ----------
    filename : str
        filename of AMBER NetCDF file.
    top : {str, Trajectory, Topology}
        The NetCDF format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    stride : int, default=None
        Only read every stride-th frame

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.NetCDFTrajectoryFile :  Low level interface to NetCDF files
    """
    topology = _parse_topology(top)
    with NetCDFTrajectoryFile(filename) as f:
        xyz, time, cell_lengths, cell_angles = f.read(stride=stride)
        xyz /= 10.0           # convert from angstroms to nanometer
        cell_lengths /= 10.0  # convert from angstroms to nanometer

    if isinstance(time, np.ma.masked_array) and np.all(time.mask):
        # if time is a masked array and all the entries are masked
        # then we just tread it as if we never found it
        time = None

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    trajectory.unitcell_lengths = cell_lengths
    trajectory.unitcell_angles = cell_angles
    return trajectory


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

    Trajectory uses the nanometers, angstroms & picoseconds unit system.

    Examples
    --------
    >>> # loading a trajectory
    >>> import mdtraj as md
    >>> md.load('trajectory.xtc', top='native.pdb')           # doctest: +SKIP
    <mdtraj.Trajectory with 1000 frames, 22 atoms at 0x1058a73d0>

    >>> # slicing a trajectory
    >>> t = md.load('trajectory.h5')                          # doctest: +SKIP
    >>> print t                                               # doctest: +SKIP
    <mdtraj.Trajectory with 100 frames, 22 atoms>
    >>> print t[::2]                                          # doctest: +SKIP
    <mdtraj.Trajectory with 50 frames, 22 atoms>

    >>> # calculating the average distance between two atoms
    >>> import mdtraj as md
    >>> import numpy as np
    >>> t = md.load('trajectory.h5')                                              # doctest: +SKIP
    >>> np.mean(np.sqrt(np.sum((t.xyz[:, 0, :] - t.xyz[:, 21, :])**2, axis=1)))   # doctest: +SKIP

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
        if hasattr(self, 'topology'):
            shape = (None, self.topology._numAtoms, 3)
        else:
            shape = (None, None, 3)

        value = ensure_type(value, np.float32, 3, 'xyz', shape=shape,
                            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        self._xyz = value

    def __len__(self):
        return self.n_frames

    def __add__(self, other):
        "Concatenate two trajectories"
        return self.join(other)

    def __str__(self):
        return "<mdtraj.Trajectory with %d frames, %d atoms>" % (self.n_frames, self.n_atoms)

    def __repr__(self):
        return "<mdtraj.Trajectory with %d frames, %d atoms at 0x%02x>" % (self.n_frames, self.n_atoms, id(self))


    def join(self, other, check_topology=True, discard_overlapping_frames=False):
        """Join two trajectories together

        This method can also be called by using `self + other`

        Parameters
        ----------
        other : Trajectory
            The other trajectory to join
        check_topology : bool
            Ensure that the topology of `self` and `other` are identical before
            joining them. If false, the resulting trajectory will have the
            topology of `self`.
        discard_overlapping_frames : bool, optional
            If True, compare coordinates at trajectory edges to discard overlapping
            frames.  Default: False.
        """
        if not isinstance(other, Trajectory):
            raise TypeError('You can only add two Trajectory instances')

        if self.n_atoms != other.n_atoms:
            raise ValueError('Number of atoms in self (%d) is not equal '
                'to number of atoms in other (%d)' % (self.n_atoms, other.n_atoms))

        if check_topology:
            if self.topology != other.topology:
                raise ValueError('The topologies of the two Trajectories are not the same')

        lengths2 = None
        angles2 = None
        other_has_unitcell = (other.unitcell_lengths is not None and other.unitcell_angles is not None)

        if discard_overlapping_frames:
            x0 = self.xyz[-1]
            x1 = other.xyz[-1]
            if np.linalg.norm(x1 - x0) < 1e-8:
                xyz = other.xyz[1:]
                time = other.time[1:]
                if other_has_unitcell:
                    lengths2 = other.unitcell_lengths[1:]
                    angles2 = other.unitcell_angles[1:]
        else:
            xyz = other.xyz
            time = other.time
            if other_has_unitcell:
                lengths2 = other.unitcell_lengths
                angles2 = other.unitcell_angles

        xyz = np.concatenate((self.xyz, xyz))
        time = np.concatenate((self.time, time))

        if self.unitcell_lengths is None and self.unitcell_angles is None and not other_has_unitcell:
            lengths, angles = None, None
        elif self.unitcell_lengths is not None and self.unitcell_angles is not None and other_has_unitcell:
            lengths = np.concatenate((self.unitcell_lengths, lengths2))
            angles = np.concatenate((self.unitcell_angles, angles2))
        else:
            raise ValueError("One trajectory has box size, other doesn't. "
                             "I don't know what to do")

        # use this syntax so that if you subclass Trajectory,
        # the subclass's join() will return an instance of the subclass
        return self.__class__(xyz, deepcopy(self._topology), time=time,
            unitcell_lengths=lengths, unitcell_angles=angles)

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

        # box has no default, it'll just be none normally
        self.unitcell_lengths = unitcell_lengths
        self.unitcell_angles = unitcell_angles

        # time will take the default 1..N
        if time is None:
            time = np.arange(len(self.xyz))
        self.time = time

        if not topology._numAtoms == self.n_atoms:
            raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (self.n_atoms, topology._numAtoms))

    def openmm_positions(self, frame):
        """OpenMM-compatable positions of a single frame.

        Examples
        --------
        >>> t = md.load('trajectory.h5')                      # doctest: +SKIP
        >>> context.setPositions(t.openmm_positions(0))       # doctest: +SKIP

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
        # copied from Lee-Ping Wang's Molecule.py
        if not HAVE_OPENMM:
            raise ImportError('OpenMM was not imported')

        Pos = []
        for xyzi in self.xyz[frame]:
            Pos.append(Vec3(xyzi[0], xyzi[1], xyzi[2]))

        return Pos * nanometer

    def openmm_boxes(self, frame):
        """OpenMM-compatable box vectors of a single frame.

        Examples
        --------
        >>> t = md.load('trajectory.h5')                          # doctest: +SKIP
        >>> context.setPeriodicBoxVectors(t.openmm_positions(0))  # doctest: +SKIP

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
        # copied from Lee-Ping Wang's Molecule.py
        if not HAVE_OPENMM:
            raise ImportError('OpenMM was not imported')

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
          '.ncdf': self.save_netcdf}

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
        topology = self.topology

        # convert to angstroms
        if self.unitcell_lengths is not None and self.unitcell_angles is not None:
            a, b, c = 10*self.unitcell_lengths[0]

            if (np.abs(self.unitcell_angles[0, 0] - 90.0) > 1e-8 or
                    np.abs(self.unitcell_angles[0, 1] - 90.0) > 1e-8 or
                    np.abs(self.unitcell_angles[0, 2] - 90.0) > 1e-8):
                warnings.warn('Unit cell information not saved correctly to PDB')

            topology.setUnitCellDimensions((a, b, c))

        with PDBTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            for i in xrange(self.n_frames):
                # need to convert internal nm to angstroms for output
                f.write(self._xyz[i] * 10, topology, modelIndex=i)


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
        # convert from internal nm representation to angstroms for output
        xyz = self.xyz * 10
        lengths = self.unitcell_lengths
        if lengths is not None:
            lengths = lengths * 10

        with DCDTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz, lengths, self.unitcell_angles)


    def save_binpos(self, filename, force_overwrite=True):
        """Save trajectory to AMBER BINPOS format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        # convert from internal nm representation to angstroms for output
        xyz = self.xyz * 10
        with BINPOSTrajectoryFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(xyz)


    def save_netcdf(self, filename, force_overwrite=True):
        """Save trajectory in AMBER NetCDF format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        from mdtraj.netcdf import NetCDFFile

        xyz = self.xyz * 10
        with NetCDFFile(filename, 'w', force_overwrite=force_overwrite) as f:
            f.write(coordinates=xyz, time=self.time)

    def center_coordinates(self):
        """Remove the center of mass from each frame in trajectory.

        This method acts inplace on the trajectory
        """
        for x in self._xyz:
            x -= (x.astype('float64').mean(0))

    def restrict_atoms(self, atom_indices):
        """Retain only a subset of the atoms in a trajectory (inplace)

        Deletes atoms not in `atom_indices`, and re-indexes those that remain

        Parameters
        ----------
        atom_indices : list([int])
            List of atom indices to keep.
        """
        self.top.restrict_atoms(atom_indices)
        self._xyz = self.xyz[:,atom_indices]


##############################################################################
# Global hack for registration of loaders.
##############################################################################

_LoaderRegistry = {
    '.xtc': load_xtc,
    '.xml': load_xml,
    '.trr': load_trr,
    '.pdb': load_pdb,
    '.dcd': load_dcd,
    '.h5': load_hdf5,
    #'.lh5': _load_legacy_hdf,
    '.binpos': load_binpos,
    '.ncdf': load_netcdf,
    '.nc': load_netcdf
}