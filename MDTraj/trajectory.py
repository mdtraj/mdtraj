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

import os
import numpy as np
import tables
from mdtraj import dcd, xtc, binpos, trr
from mdtraj.pdb import pdbfile
from mdtraj import io
from mdtraj.utils import unitcell
import mdtraj.topology
import functools
from itertools import izip
from copy import deepcopy
try:
    from simtk.openmm import Vec3
    from simtk.unit import nanometer
    import simtk.openmm.app.topology
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

import logging
logger = logging.getLogger(__name__)

MAXINT16 = np.iinfo(np.int16).max
MAXINT32 = np.iinfo(np.int32).max
DEFAULT_PRECISION = 1000


def _convert_to_lossy_integers(X, precision=DEFAULT_PRECISION):
    """Implementation of the lossy compression used in Gromacs XTC using
    the pytables library.  Convert 32 bit floats into 16 bit integers.
    These conversion functions have been optimized for memory use.
    Further memory reduction would require an in-place astype() operation,
    which one could create using ctypes."""

    if np.max(X) * float(precision) < MAXINT16 and np.min(X) * float(precision) > -MAXINT16:
        X *= float(precision)
        Rounded = X.astype("int16")
        X /= float(precision)
    else:
        X *= float(precision)
        Rounded = X.astype("int32")
        X /= float(precision)
        logger.error("Data range too large for int16: try removing center of mass motion, check for 'blowing up, or just use .h5 or .xtc format.'")
    return(Rounded)


def _convert_from_lossy_integers(X, precision=DEFAULT_PRECISION):
    """Implementation of the lossy compression used in Gromacs XTC using
    the pytables library.  Convert 16 bit integers into 32 bit floats."""
    X2 = X.astype("float32")
    X2 /= float(precision)
    return(X2)


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
        topology = pdbfile.PDBFile(top).topology
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, mdtraj.topology.Topology):
        topology = top
    elif HAVE_OPENMM and isinstance(top, simtk.openmm.app.topology.Topology):
        topology = top
    else:
        raise TypeError('Could not interpreted top=%s' % top)

    return topology


def load(filename_or_filenames, discard_overlapping_frames=False, **kwargs):
    """Load a trajectory from a file or list of files

    Parameters
    ----------
    filename_or_filenames : {str, list of strings}
        filename or list of filenames containing trajectory files of a single format.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them


    Other Parameters
    -------------------------
    DEPENDS ON THE LOADER. NOT FILLED IN YET. NEED TO ADD STRIDING AND SUCH

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
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

    # We have only a single trajectory now.
    loaders = {'.xtc':    load_xtc,
               '.trr':    load_trr,
               '.pdb':    load_pdb,
               '.dcd':    load_dcd,
               '.h5':     load_hdf,
               '.lh5':    load_hdf,
               '.binpos': load_binpos}

    try:
        loader = loaders[extension]
    except KeyError:
        raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                      'was found. I can only load files '
                      'with extensions in %s' % (filename, extension, loaders.keys()))

    return loader(filename, **kwargs)

def load_pdb(filename):
    """Load a pdb file.

    Parameters
    ----------
    filename : str
        Filename of PDB

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """
    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_pdb. '
            'you supplied %s' % type(filename))

    filename = str(filename)
    f = pdbfile.PDBFile(filename)

    # convert from angstroms to nm
    coords = f.positions / 10

    trajectory = Trajectory(xyz=coords, topology=f.topology)

    if f.topology.getUnitCellDimensions() is not None:
        a, b, c = f.topology.getUnitCellDimensions()
        unitcell_params = np.recarray((1,),
                dtype=[('a', '<f8'), ('b', '<f8'), ('c', '<f8'),
                       ('alpha', '<f8'), ('beta', '<f8'), ('gamma', '<f8')])
        unitcell_params[0] = (a, b, c, 90.0, 90.0, 90.0)

        trajectory.unitcell_parameters = unitcell_params

    return trajectory

def load_xtc(filename, top=None, chunk=500):
    """Load an xtc file. Since the xtc doesn't contain information
    to specify the topolgy, you need to supply the topology yourself

    Parameters
    ----------
    filename : str
        Filename (string) of xtc trajectory.
    top : {str, Trajectory, Topology}
        The XTC format does not contain topology information. Pass in either the
        path to a pdb file, a trajectory, or a topology to supply this information.
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
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

    xyz, time, step, box, prec = xtc.read(filename, chunk)
    # note we're not doing anything with the box vectors
    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    trajectory.unitcell_vectors = box

    return trajectory


def load_trr(filename, top=None, chunk=500):
    """Load a trr file. Since the trr doesn't contain information
    to specify the topolgy, you need to supply the topology yourself

    Parameters
    ----------
    filename : str
        Filename of TRR trajectory file.
    top : {str, Trajectory, Topology}
        The TRR format does not contain topology information. Pass in either the
        path to a pdb file, a trajectory, or a topology to supply this information.
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
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

    xyz, time, step, box, lambd = trr.read(filename, chunk)

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
    trajectory : Trajectory
        A trajectory file!
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
    xyz, box_length, box_angle = dcd.read(filename)

    xyz /= 10.  # convert from anstroms to nanometer
    box_length /= 10.  # convert from anstroms to nanometer

    trajectory = Trajectory(xyz=xyz, topology=topology)
    unitcell = np.recarray((len(xyz),), dtype=[('a', '<f8'), ('b', '<f8'), ('c', '<f8'),
                                               ('alpha', '<f8'), ('beta', '<f8'), ('gamma', '<f8')])
    unitcell['a'] = box_length[:, 0]
    unitcell['b'] = box_length[:, 1]
    unitcell['c'] = box_length[:, 2]
    unitcell['alpha'] = box_angle[:, 0]
    unitcell['beta'] = box_angle[:, 1]
    unitcell['gamma'] = box_angle[:, 2]
    trajectory.unitcell_parameters = unitcell

    return trajectory


def load_hdf(filename, top=None, stride=None, frame=None, chunk=50000,
             upconvert_int16=True):
    """
    Load an hdf file.

    Parameters
    ----------
    filename : str
        String filename of HDF Trajectory file.
    top : topology
        Replace the topology in the file with this topology
    stride : int, default=None
        Only read every stride-th frame
    frame : {None, int}, default=None
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
    chunk : int, default=50000
       Size of the chunk to use for loading the file.
    upconvert_int16 : bool, default=True
        If the data is encoded in the file as int16, an automatic upconversion
        to float32 based on the gromacs lossy encoding scheme is done.

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_hdf. '
            'you supplied %s' % type(filename))

    F = tables.File(filename, 'r')
    try:
        # get the topology from the file, or input argument list
        if top is None:
            if hasattr(F.root, 'topology'):
                topology = mdtraj.topology.from_bytearray(F.root.topology[:])
            else:
                raise ValueError("I can't read the topology from your old HDF "
                    "file, so how about you give me a new topology with the top "
                    "optional argument")
        else:
            topology = _parse_topology(top)

        # find the xyz coordinates in the file
        if hasattr(F.root, 'xyz'):
            xyz_node = F.root.xyz
        elif hasattr(F.root, 'XYZList'):
            xyz_node = F.root.XYZList
        else:
            raise ValueError('XYZ data not found in %s' % filename)

        # find the time in the file
        time_node = None
        if hasattr(F.root, 'time'):
            time_node = F.root.time

        # find box in the file
        box_node = None
        if hasattr(F.root, 'box'):
            box_node = F.root.box

        if chunk < stride:
            raise ValueError('Chunk must be greater than or equal to the stride')

        # if they want only a single frame, quit early without
        # doing the chunking
        if frame is not None:
            time = None
            if time_node is not None:
                time = time_node[int(frame)]

            box = None
            if box_node is not None:
                box = box_node[int(frame)]

            xyz = xyz_node[int(frame)]
            if upconvert_int16 and xyz.dtype == np.int16:
                xyz = _convert_from_lossy_integers(xyz)

            trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
            trajectory.unitcell_vectors = box
            return trajectory

        # adjust the stride/chunk
        if stride is not None:
            # Need to do this in order to make sure we stride correctly.
            # since we read in chunks, and then we need the strides
            # to line up
            while chunk % stride != 0:
                chunk -= 1
        else:
            stride = 1


        shape = xyz_node.shape
        begin_range_list = np.arange(0, shape[0], chunk)
        end_range_list = np.concatenate((begin_range_list[1:], [shape[0]]))
        def enum_chunks():
            for r0, r1 in zip(begin_range_list, end_range_list):
                xyz = np.array(xyz_node[r0: r1: stride])

                time = None
                if time_node is not None:
                    time = np.array(time_node[r0: r1: stride])

                box = None
                if box_node is not None:
                    box = np.array(box_node[r0: r1: stride])

                yield xyz, time, box

        zipper = tuple(izip(*enum_chunks()))
        xyz = np.concatenate(zipper[0])

        time = None
        if time_node is not None:
            time = np.concatenate(zipper[1])

        box = None
        if box_node is not None:
            box = np.concatenate(zipper[2])

        if upconvert_int16 and xyz.dtype == np.int16:
            xyz = _convert_from_lossy_integers(xyz)
    except:
        raise
    finally:
        F.close()

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    trajectory.unitcell_vectors = box
    return trajectory


def load_binpos(filename, top=None, chunk=500):
    """Load an AMBER binpos file. Since the dcd format doesn't contain
    information to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filename : str
        String filename of Amber binpos file.
    top : {str, Trajectory, Topology}
        DCD XTC format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
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

    xyz = binpos.read(filename)
    xyz /= 10.0  # convert from anstroms to nanometer

    return Trajectory(xyz=xyz, topology=topology)


class Trajectory(object):
    """
    Attributes
    ----------
    n_frames : int
    n_atoms : int
    time : np.ndarray dims=1, shape=(n_frames)
    timestep
    topology : pdb.topology
    top : pdb.topology
       Alias for topology
    """

    @property
    def n_frames(self):
        return self._xyz.shape[0]

    @property
    def n_atoms(self):
        return self._xyz.shape[1]

    @property
    def top(self):
        return self.topology

    @property
    def timestep(self):
        return self._time[1] - self._time[0]

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
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
        """Get the vectors that define the shape of the unit cell
        in each frame

        Returns
        -------
        vectors : np.ndarray, shape(n_frames, 3, 3)
            Vectors definiing the shape of the unit cell in each frame.
            The semantics of this array are that the shape of the unit cell
            in frame `i` are given by the three vectors, value[i, 0, :],
            `value[i, 1, :]`, and `value[i, 2, :]`.
        """
        if self._unitcell_parameters is None:
            return None

        v1, v2, v3 = unitcell.lengths_and_angles_to_box_vectors(
            self._unitcell_parameters['a'],
            self._unitcell_parameters['b'],
            self._unitcell_parameters['c'],
            self._unitcell_parameters['alpha'] * np.pi / 180,
            self._unitcell_parameters['beta']  * np.pi / 180,
            self._unitcell_parameters['gamma']  * np.pi / 180,
        )
        return np.swapaxes(np.dstack((v1, v2, v3)), 1, 2)

    @unitcell_vectors.setter
    def unitcell_vectors(self, vectors):
        """Set the three vectors that define the shape of the unit cell

        Parameters
        ----------
        vectors : tuple of three arrays, each of shape=(n_frames, 3)

            The semantics of this array are that the shape of the unit cell
            in frame `i` are given by the three vectors, value[i, 0, :],
            `value[i, 1, :]`, and `value[i, 2, :]`.
        """
        if vectors is None:
            self._unitcell_parameters = None
            return

        if not len(vectors) == len(self):
            raise TypeError('unitcell_vectors must be the same length as '
                            'the trajectory. you provided %s' % vectors)

        v1 = vectors[:, 0, :]
        v2 = vectors[:, 1, :]
        v3 = vectors[:, 2, :]
        a, b, c, alpha, beta, gamma = unitcell.box_vectors_to_lengths_and_angles(v1, v2, v3)

        self._unitcell_parameters = np.recarray((len(self),),
            dtype=[('a', '<f8'), ('b', '<f8'), ('c', '<f8'),
                   ('alpha', '<f8'), ('beta', '<f8'), ('gamma', '<f8')])

        self._unitcell_parameters['a'] = a
        self._unitcell_parameters['b'] = b
        self._unitcell_parameters['c'] = c
        self._unitcell_parameters['alpha'] = alpha * 180 / np.pi
        self._unitcell_parameters['beta'] = beta * 180 / np.pi
        self._unitcell_parameters['gamma'] = gamma * 180 / np.pi

    @property
    def unitcell_parameters(self):
        """Get the six parameters, (a, b, c, \alpha, \beta, \gamma) that define
        the shape of the unit cell in each frame

        Returns
        -------
        parameters : np.recarray, dype=[('a', float), ('b', float), ('c', float), ('alpha', float), ('beta', float), ('gamma', float)]
            Parameters describing the shape of the unit cell in each frame. The
            return value is a six field  record array. The records 'a', 'b', and 'c'
            give the length of the three unit cell vectors, 'alpha' gives the
            angle between vectors **b** and **c**, beta gives the angle between
            vectors **c** and **a**, and gamma gives the angle between vectors
            **a** and **b**. The angles are in degrees.
        """
        return self._unitcell_parameters

    @unitcell_parameters.setter
    def unitcell_parameters(self, value):
        """Set the six parameters, (a, b, c, \alpha, \beta, \gamma) that define
        the shape of the unit cell in each frame

        Parameters
        ----------
        value : np.recarray, dype=[('a', float), ('b', float), ('c', float), ('alpha', float), ('beta', float), ('gamma', float)]
        """
        if value is None:
            self._unitcell_parameters = None
            return


        correct_dtype = np.dtype([('a', '<f8'), ('b', '<f8'), ('c', '<f8'),
                                  ('alpha', '<f8'), ('beta', '<f8'), ('gamma', '<f8')])
        if not value.dtype == correct_dtype:
            raise TypeError('unitcell_parameters must be a six field recarray '
                'with dtype %s' % correct_dtype)

        if isinstance(value, np.core.records.record) and len(self) == 1:
            # if this is a length 1 array, then you can provide a single
            # record instead of an array of records. this helps make slicing
            # work sensibly
            value = np.array([value], dtype=value.dtype)

        if not len(value) == len(self):
            raise TypeError('unitcell_parameters must be the same length as '
                            'the trajectory. you provided %s' % value)

        if np.all(value['alpha'] < 2*np.pi) and np.all(value['beta'] < 2*np.pi) and np.all(value['gamma'] < 2*np.pi):
            warnings.warn('All of the angles you provided are less than 2*pi. '
                'Did you give me radians? I want degrees.')

        self._unitcell_parameters = value

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, value):
        #TODO control the precision of xyz
        if value.ndim == 2:
            n_frames = 1
            n_atoms, n_dims = value.shape
            assert n_dims == 3

            xyz = value.reshape((n_frames, n_atoms, n_dims))
        elif value.ndim == 3:
            xyz = value
        else:
            raise ValueError('xyz is wrong shape')

        if hasattr(self, 'topology') and self.topology._numAtoms != xyz.shape[1]:
            raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (xyz.shape[1], self.topology._numAtoms))

        self._xyz = xyz

    def __len__(self):
        return self.n_frames

    def __add__(self, other):
        "Concatenate two trajectories"
        return self.join(other)

    def join(self, other, check_topology=True, discard_overlapping_frames=False):
        """
        Join two trajectories together

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
            if not mdtraj.topology.equal(self.topology, other.topology):
                raise ValueError('The topologies are not the same')

        unitcell_parameters2 = None

        if discard_overlapping_frames:
            x0 = self.xyz[-1]
            x1 = other.xyz[-1]
            if np.linalg.norm(x1 - x0) < 1e-8:
                xyz = other.xyz[1:]
                time = other.time[1:]
                if other.unitcell_parameters is not None:
                    unitcell_parameters2 = other.unitcell_parameters[1:]
        else:
            xyz = other.xyz
            time = other.time
            if other.unitcell_parameters is not None:
                unitcell_parameters2 = other.unitcell_parameters

        xyz = np.concatenate((self.xyz, xyz))
        time = np.concatenate((self.time, time))

        if self.unitcell_parameters is None and unitcell_parameters2 is None:
            unitcell_parameters = None
        elif self.unitcell_parameters is not None and unitcell_parameters2 is not None:
            unitcell_parameters = np.concatenate((self.unitcell_parameters,
                                                  unitcell_parameters2))
        else:
            raise ValueError("One trajectory has box size, other doesn't. "
                             "I don't know what to do")

        # use this syntax so that if you subclass Trajectory,
        # the subclass's join() will return an instance of the subclass
        return self.__class__(xyz, deepcopy(self.topology), time, unitcell_parameters)

    def __getitem__(self, key):
        "Get a slice of this trajectory"
        return self.slice(key)

    def slice(self, key, copy=True):
        """
        Slice a trajectory, by extracting one or more frames a separate traj

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
        if self.unitcell_parameters is not None:
            unitcell_parameters = self.unitcell_parameters[key]
        else:
            unitcell_parameters = None

        if copy:
            xyz = xyz.copy()
            time = time.copy()
            topology = deepcopy(self.topology)

            if self.unitcell_parameters is not None:
                unitcell_parameters = unitcell_parameters.copy()

        newtraj = self.__class__(xyz, topology, time, unitcell_parameters)
        return newtraj


    def __init__(self, xyz, topology, time=None, unitcell_parameters=None):
        self.xyz = xyz
        self.topology = topology

        # box has no default, it'll just be none normally
        self._unitcell_parameters = None
        if unitcell_parameters is not None:
            self.unitcell_parameters = unitcell_parameters

        # time will take the default 1..N
        if time is None:
            time = np.arange(len(self.xyz))
        self.time = time

        if not topology._numAtoms == self.n_atoms:
            raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (self.n_atoms, topology._numAtoms))

    def openmm_positions(self, frame):
        """
        Return OpenMM compatable positions of a single frame.

        Parameters
        ----------
        frame : int
            Which trajectory frame to return.

        Returns
        -------
        positions : list of XYZ coordinates of specific trajectory frame, formatted
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
        """Return OpenMM compatable box vectors of a single frame.

        Parameters
        ----------
        frame : int
            Return box for this single frame.

        Returns
        -------
        box : list of XYZ coordinates of periodic box vectors, formatted
            for input to OpenMM
        """
        # copied from Lee-Ping Wang's Molecule.py
        if not HAVE_OPENMM:
            raise ImportError('OpenMM was not imported')

        if self.unitcell_parameters is None:
            raise ValueError("this trajectory does not contain box size information")

        params = self.unitcell_parameters[frame]
        v1, v2, v3 = lengths_and_angles_to_box_vectors(params['a'], params['b'], params['c'],
            params['alpha']*np.pi/180, params['beta']*np.pi/180, params['gamma']*np.pi/180)

        return (Vec3(*v1), Vec3(*v2), Vec3(*v3)) * nanometer


    @staticmethod
    # im not really sure if the load function should be just a function or a method on the class
    # so effectively, lets make it both?
    def load(filenames, **kwargs):
        """
        Load a trajectory

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
        """
        Save trajectory to a file

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory. The extension will be parsed and will
            control the format.

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
                  '.h5': self.save_hdf,
                  '.lh5': self.save_hdf,
                  '.binpos': self.save_binpos}

        try:
            saver = savers[extension]
        except KeyError:
            raise IOError('Sorry, no saver for filename=%s (extension=%s) '
                          'was found. I can only save files '
                          'with extensions in %s' % (filename, extension, savers.keys()))

        # run the saver, and return whatever output it gives
        return saver(filename)

    def save_hdf(self, filename, lossy=True):
        """
        Save trajectory to an hdf5

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        lossy : bool
            Use gromacs style compression
        """
        if lossy:
            xyz=_convert_to_lossy_integers(self.xyz)

        kwargs = {
            'xyz': xyz,
            'time': self.time,
            'topology': mdtraj.topology.to_bytearray(self.topology)
        }
        if self.box is not None:
            kwargs['box'] = self.box

        return io.saveh(filename, **kwargs)

    def save_pdb(self, filename, no_models=False):
        """
        Save a trajectory to pdb

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        no_models : bool
            TODO: Document this feature. What does it do?
        """

        f = open(filename, 'w')

        topology = self.topology

        # convert to angstroms
        if self.topology.getUnitCellDimensions() is not None:
            topology.setUnitCellDimensions(tuple([10*i for i in self.topology.getUnitCellDimensions()]))

        pdbfile.PDBFile.writeHeader(topology, file=f)

        for i in xrange(self.n_frames):
            if no_models:
                mind = None
            else:
                mind = i

            # need to convert internal nm to angstroms for output
            positions = [ list(self._xyz[i,j,:].flatten()*10) for j in xrange(self.n_atoms) ]
            pdbfile.PDBFile.writeModel(topology,
                                      positions,
                                      file=f,
                                      modelIndex=mind)

        pdbfile.PDBFile.writeFooter(topology, file=f)
        f.close()

        return

    def save_xtc(self, filename, force_overwrite=True):
        """
        Save a trajectory to gromacs XTC format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        return xtc.write(filename, self.xyz, time=self.time,
            force_overwrite=force_overwrite)

    def save_trr(self, filename, force_overwrite=True):
        """
        Save a trajectory to gromacs TRR format

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
        return trr.write(filename, self.xyz, time=self.time,
            force_overwrite=force_overwrite)

    def save_dcd(self, filename, force_overwrite=True):
        """
        Save a trajectory to CHARMM dcd format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        # convert from internal nm representation to angstroms for output
        xyz = self.xyz * 10
        return dcd.write(filename, xyz, force_overwrite=force_overwrite)

    def save_binpos(self, filename, force_overwrite=True):
        """
        Save a trajectory to amber BINPOS format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        """
        # convert from internal nm representation to angstroms for output
        xyz = self.xyz * 10
        return binpos.write(filename, xyz, force_overwrite=force_overwrite)
