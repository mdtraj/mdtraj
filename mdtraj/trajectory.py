import os
import numpy as np
from numbers import Number

from xtc import xtc
from dcd import dcd
from pdb import pdbfile
from topology import Topology
import io

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
        Look for overlapping frames and discard them
    
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

    n_redundant = 0
    xyz = []
    for i, c in enumerate(xtc.XTCReader(filenames)):
        if discard_overlapping_frames:
            if i > 0:
                if np.sum(np.abs(c.coords - xyz[-1])) < 10. ** -8:
                    n_redundant += 1
                    xyz.append(np.array(c.coords).copy())
                    
        else:
            xyz.append(np.array(c.coords).copy())
        
    xyz = np.array(xyz)

    # TODO: Get the timestep
    timestep = None
    
    return Trajectory(xyz=xyz, topology=topology, timestep=timestep)
    


def load_pdb(filename, **kwargs):
    """Load a pdb file.

    Parameters
    ----------
    filename : str
        Path to the pdb file on disk

    Optional Arguments
    ------------------

    Returns
    -------
    trajectory : Trajectory
        A trajectory file!
    """
    
    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_pdb. you supplied %s' % type(filename))
    filename = str(filename)
    f = pdbfile.PDBFile(filename)
    return Trajectory(xyz=f.positions, topology=f.topology)


def load_dcd(filename):
    pass


def load_hdf(filename):
    pass


class Trajectory(object):
    """
    Attributes
    ----------
    n_frames : int
    n_atoms : int
    timestep : int
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
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        if (isinstance(value, Number) and value > 0) or value is None:
            self._timestep = value
        else:
            raise TypeError('Timestep must be a positive number or None')

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

    def __init__(self, xyz, topology, timestep=None):
        self.xyz = xyz
        self.topology = topology
        self.timestep = timestep
        
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
        # timestep needs to be an array for io to work
        if self.timestep is None:
            ts_array = np.array([-1])
        else:
            ts_array = np.array(self.timestep)

        #TODO: control the precision of xyz

        return io.saveh(filename, xyz=self.xyz, timestep=ts_array, topology=self.topology.to_bytearray())

    def join(self, other_trajectory):
        "Join two trajectories together"
        pass

if __name__ == '__main__':
    #load('/home/robert/msmbuilder/Tutorial/native.pdbd')
    traj = load('/home/robert/msmbuilder/Tutorial/native.pdb')
    print load('/home/robert/msmbuilder/Tutorial/XTC/RUN00/frame0.xtc', top=traj).top

    traj.save('file.h5')
