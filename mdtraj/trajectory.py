import os
import IPython as ip
from mdtraj.xtc import xtc
from mdtraj.dcd import dcd
from mdtraj.pdb import pdbfile

def load(filename, **kwargs):
    """Load a trajectory from a file
    
    Parameters
    ---------
    filename : str
    
    Returns
    -------
    trajectory : Trajectory
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
    
    loader(filename, **kwargs)


def load_xtc(filename):
    pass


def load_pdb(filename, **kwargs):
    f = pdbfile.PDBFile(filename)
    print f.positions



def load_dcd(filename):
    pass


def load_hdf(filename):
    pass


class Trajectory():
    pass
        

if __name__ == '__main__':
    #load('/home/robert/msmbuilder/Tutorial/native.pdbd')
    load('/home/robert/msmbuilder/Tutorial/native.pdb')
