import os
import numpy as np

from xtc import xtc
from dcd import dcd
from pdb import pdbfile
from topology import Topology

#__all__ == ['load', 'load_xtc', 'load_pdb', 'load_dcd', 'Trajectory']

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
    
    return Trajectory(xyz=xyz, topology=topology)
    


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


class Trajectory():
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
        return self.xyz.shape[0]

    @property
    def n_atoms(self):
        return self.xyz.shape[1]

    @property
    def top(self):
        return self.topology

    def __init__(self, xyz, topology, timestep=-1):
        if xyz.ndim == 2:
            n_frames = 1
            n_atoms, n_dims = xyz.shape
            assert n_dims == 3

            self.xyz = xyz.reshape((n_frames, n_atoms, n_dims))
        elif xyz.ndim == 3:
            self.xyz = xyz
        else:
            print xyz
            raise ValueError()
            
        self.topology = topology
        self.timestep = timestep

        if not topology._numAtoms == self.n_atoms:
            raise ValueError("Number of atoms in xyz (%s) and "
                "in topology (%s) don't match" % (self.n_atoms, topology._numAtoms))

    def join(self, other_trajectory):
        "Join two trajectories together"

        pass


if __name__ == '__main__':
    #load('/home/robert/msmbuilder/Tutorial/native.pdbd')
    traj = load('/home/robert/msmbuilder/Tutorial/native.pdb')
    print load('/home/robert/msmbuilder/Tutorial/XTC/RUN00/frame0.xtc', top=traj).top
