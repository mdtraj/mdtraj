import os
import numpy as np
import tables
from tables import NoSuchNodeError
from mdtraj import dcd
from mdtraj import xtc
from mdtraj import binpos
from mdtraj.pdb import pdbfile
from mdtraj import io
from topology import Topology
from itertools import izip
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
    elif isinstance(top, Topology):
        topology = top
    else:
        raise TypeError('Could not interpreted top=%s' % top)

    return topology


def load(filename, **kwargs):
    """Load a trajectory from a file or list of files

    Parameters
    ----------
    filename : {str, list of strings}
        filesystem path from which to load the trajectory

    Other Parameters
    -------------------------
    DEPENDS ON THE LOADER. NOT FILLED IN YET. NEED TO ADD STRIDING AND SUCH

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """

    _assert_files_exist(filename)
    # grab the extension of the filename
    extension = os.path.splitext(filename)[1]

    loaders = {'.xtc':    load_xtc,
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
        filesystem path from which to load the trajectory

    Other Parameters
    ----------------
    None

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
    return Trajectory(xyz=f.positions, topology=f.topology)


def load_xtc(filenames, top, discard_overlapping_frames=False, chunk=500):
    """Load an xtc file. Since the xtc doesn't contain information
    to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings giving one or multiple filenames that together form an xtx
    top : {str, Trajectory, Topology}
        The XTC format does not contain topology information. Pass in either the path to a pdb file,
        a trajectory, or a topology to supply this information.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.


    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """
    topology = _parse_topology(top)
    if isinstance(filenames, basestring):
        filenames = [filenames]

    coords = []
    times = []
    for i, filename in enumerate(filenames):
        xyz, box, time, prec, step = xtc.read(filename, chunk)
        if i > 0 and discard_overlapping_frames:
            if np.sum(np.square(xyz[0] - coords[-1][-1])) < 1e-8:
                xyz = xyz[0:-1]
                time = time[0:-1]
        coords.append(xyz)
        times.append(time)

    times = np.concatenate(time)
    coords = np.vstack(coords)

    return Trajectory(xyz=coords, topology=topology, time=times)


def load_dcd(filenames, top, discard_overlapping_frames=False):
    """Load an xtc file. Since the dcd format doesn't contain information
    to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings giving one or multiple filenames that together
        form an dcd
    top : {str, Trajectoray, Topology}
        DCD XTC format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """

    topology = _parse_topology(top)
    if isinstance(filenames, basestring):
        filenames = [filenames]

    coords = []
    for i, filename in enumerate(filenames):
        xyz = dcd.read_xyz(filename)
        if i > 0 and discard_overlapping_frames:
            if np.sum(np.square(xyz[0] - coords[-1][-1])) < 1e-8:
                xyz = xyz[0:-1]
        coords.append(xyz)

    coords = np.vstack(coords)

    return Trajectory(xyz=coords, topology=topology)


def load_hdf(filename, top=None, stride=None, chunk=50000, upconvert_int16=True):
    """
    Load an hdf file.

    Parameters
    ----------
    filename : str
        filesystem path from which to load the trajectory
    top : topology
        Replace the topology in the file with this topology
    stride : int, default=None
        Only read every stride-th frame
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

    F = tables.File(filename, 'r')
    try:
        # get the topology from the file, or input argument list
        if top is None:
            if hasattr(F.root, 'topology'):
                topology = Topology.from_bytearray(F.root.topology)
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
        if hasattr(F.root, 'time'):
            time_node = F.root.time
        else:
            time_node = None

        if chunk < stride:
            raise ValueError('Chunk must be greater than or equal to the stride')

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
                if time_node is not None:
                    time = np.array(time_node[r0: r1: stride])
                else:
                    time = None

                yield xyz, time

        zipper = tuple(izip(*enum_chunks()))
        xyz = np.concatenate(zipper[0])
        if time_node is not None:
            time = np.concatenate(zipper[1])
        else:
            time = None

        if upconvert_int16 and xyz.dtype == np.int16:
            xyz = _convert_from_lossy_integers(xyz)
    except:
        raise
    finally:
        F.close()

    return Trajectory(xyz=xyz, topology=topology, time=time)


def load_binpos(filenames, top, chunk=500, discard_overlapping_frames=False):
    """Load an AMBER binpos file. Since the dcd format doesn't contain
    information to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings giving one or multiple filenames that together
        form a binpos
    top : {str, Trajectory, Topology}
        DCD XTC format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    chunk : int, default=500
        Size of the chunk to use for loading the xtc file. Memory is allocated
        in units of the chunk size, so larger chunk can be more time-efficient.
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """

    topology = _parse_topology(top)
    if isinstance(filenames, basestring):
        filenames = [filenames]

    coords = []
    for i, filename in enumerate(filenames):
        xyz = binpos.read_xyz(filename)
        if i > 0 and discard_overlapping_frames:
            if np.sum(np.square(xyz[0] - coords[-1][-1])) < 1e-8:
                xyz = xyz[0:-1]
        coords.append(xyz)

    coords = np.vstack(coords)

    return Trajectory(xyz=coords, topology=topology)


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

            self._xyz = value.reshape((n_frames, n_atoms, n_dims))

        elif value.ndim == 3:
            self._xyz = value

        else:
            raise ValueError('xyz is wrong shape')

    def __len__(self):
        return len(self.xyz)

    def __add__(self, other):
        "Concatenate two trajectories"
        raise NotImplementedError

    def __iadd__(self, other):
        "Concatenate two trajectories, override the += operator"
        raise NotImplementedError

    def __getitem__(self, key):
        "Get a slice of this trajectory"
        raise NotImplementedError

    def __init__(self, xyz, topology, time=None):
        self.xyz = xyz
        self.topology = topology

        if time is None:
            time = np.arange(len(xyz))
        self._time = time

        if not topology._numAtoms == self.n_atoms:
            raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (self.n_atoms, topology._numAtoms))


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
        """
        # grab the extension of the filename
        extension = os.path.splitext(filename)[1]

        savers = {'.xtc': self.save_xtc,
                  '.pdb': self.save_pdb,
                  '.dcd': self.save_dcd,
                  '.h5': self.save_hdf,
                  '.binpos': self.save_binpos}

        try:
            saver = savers[extension]
        except KeyError:
            raise IOError('Sorry, no saver for filename=%s (extension=%s) '
                          'was found. I can only load files '
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

        return io.saveh(filename,
            xyz=xyz,
            time=self.time,
            topology=self.topology.to_bytearray())

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
        pdbfile.PDBFile.writeHeader(self.topology, file=f)

        for i in range(self._xyz.shape[0]):
            if no_models:
                mind = None
            else:
                mind = i
            positions = [ list(self._xyz[i,j,:].flatten()) for j in range(self._xyz.shape[1]) ]
            pdbfile.PDBFile.writeModel(self.topology,
                                      positions,
                                      file=f,
                                      modelIndex=mind)

        pdbfile.PDBFile.writeFooter(self.topology, file=f)
        f.close()

        return

    def save_xtc(self, filename):
        """
        Save a trajectory to gromacs XTC format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        """
        raise NotImplementedError

    def save_dcd(self, filename):
        """
        Save a trajectory to CHARMM dcd format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        """
        raise NotImplementedError

    def save_binpos(self, filename):
        """
        Save a trajectory to amber BINPOS format

        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        """

        raise NotImplementedError
