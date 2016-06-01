import mdtraj as md
from mdtraj.formats import TNGTrajectoryFile
from mdtraj.testing.testing import *

test_fn = get_fn('tng_example.tng')

def test_load_trajectory():
    """Compare a TNG file to the PDB file it was created from."""
    
    pdbtraj = md.load_pdb(get_fn('frame0.pdb'))
    tngtraj = md.load_tng(get_fn('frame0.tng'), top=pdbtraj.topology)
    eq(pdbtraj.n_frames, tngtraj.n_frames)
    eq(pdbtraj.unitcell_vectors, tngtraj.unitcell_vectors)
    eq(pdbtraj.xyz, tngtraj.xyz)
