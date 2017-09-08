from mdtraj import load, rmsd
import numpy as np

def test_make_whole(get_fn):
    # Check make_whole for issue #1274
    t = load(get_fn('t783_orig.pdb'))
    t2 = load(get_fn('t783_whole.pdb'))
    t.make_molecules_whole()
    r = rmsd(t, t2)[0]
    np.testing.assert_almost_equal(r, 0.0, decimal=3) 
    
