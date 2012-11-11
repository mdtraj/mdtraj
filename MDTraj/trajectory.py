import os
import numpy as np
import xtc
import dcd
from pdb import pdbfile
from topology import Topology
import io
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

def load(filename, **kwargs):
    """Load a trajectory from a file or list of files

    Parameters
    ---------
    filename : {str, list of strings}
       The file or list of files to load from

    Optional Keyword Aruments
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

    loaders = {'.xtc': load_xtc,
               '.pdb': load_pdb,
               '.dcd': load_dcd,
               '.h5': load_hdf}

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
        Path to the pdb file on disk

    Optional Arguments
    ------------------
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


def load_xtc(filenames, top, **kwargs):
    """Load an xtc file. Since the xtc doesn't contain information
    to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings giving one or multiple filenames that together form an xtx
    top : {str, Trajectory, Topology}
        The XTC format does not contain topology information. Pass in either the path to a pdb file,
        a trajectory, or a topology to supply this information.

    Optional Arguments
    ------------------
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

    discard_overlapping_frames = kwargs.pop('discard_overlapping_frames', False)
    chunk = kwargs.pop('chunk', 500)

    if isinstance(top, basestring):
        topology = pdbfile.PDBFile(top).topology
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, Topology):
        topology = top
    else:
        raise TypeError('Could not interpreted top=%s' % top)

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


def load_dcd(filenames, top, **kwargs):
    """Load an xtc file. Since the dcd format doesn't contain information
    to specify the topolgy, you need to supply a pdb_filename

    Parameters
    ----------
    filenames : {str, [str]}
        String or list of strings giving one or multiple filenames that together
        form an xtx
    top : {str, Trajectory, Topology}
        DCD XTC format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.

    Optional Arguments
    ------------------
    discard_overlapping_frames : bool, default=False
        Look for overlapping frames between the last frame of one filename and
        the first frame of a subsequent filename and discard them

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """

    discard_overlapping_frames = kwargs.pop('discard_overlapping_frames', False)

    if isinstance(top, basestring):
        topology = pdbfile.PDBFile(top).topology
    elif isinstance(top, Trajectory):
        topology = top.topology
    elif isinstance(top, Topology):
        topology = top
    else:
        raise TypeError('Could not interpreted top=%s' % top)

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


def load_hdf(filename, top=None):
    """Load an hdf file.

    Parameters
    ----------
    filename : str
        Path to the file

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """
    container = io.loadh(filename, deferred=False)
    if 'xyz' in container:
        xyz = container['xyz']
    elif 'XYZList' in container:
        xyz = container['XYZList']
    else:
        raise IOError('xyz data not found in %s' % filename)

    if 'topology' in container:
        topology = Topology.from_bytearray(container['topology'])
    else:
        if top is None:
            raise ValueError("I can't read the topology from your old HDF "
                "file, so how about you give me a new topology with the top "
                "optional argument")
        topology = top


    if 'time' in container:
        time = container['time']
    else:
        time = None

    if xyz.dtype == np.int16:
        xyz = _convert_from_lossy_integers(xyz)

    return Trajectory(xyz=xyz, topology=topology, time=time)


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
    def load(filename, **kwargs):
        return load(filename, **kwargs)

    def save(self, filename):
        "Save trajectory to a file"
        # grab the extension of the filename
        extension = os.path.splitext(filename)[1]

        savers = {#'.xtc': self.save_xtc,
                  #'.pdb': self.save_pdb,
                  #'.dcd': self.save_dcd,
                  '.h5': self.save_hdf}

        try:
            saver = savers[extension]
        except KeyError:
            raise IOError('Sorry, no loader for filename=%s (extension=%s) '
                          'was found. I can only load files '
                          'with extensions in %s' % (filename, extension, savers.keys()))

        # run the saver, and return whatever output it gives
        return saver(filename)

    def save_hdf(self, filename):
        return io.saveh(filename,
            xyz=_convert_to_lossy_integers(self.xyz),
            time=self.time,
            topology=self.topology.to_bytearray())

    def join(self, other_trajectory):
        "Join two trajectories together"
        pass
