from mdtraj.trajectory import load
from mdtraj.testing import get_fn, eq

def test_load_multiframe():
    t = load(get_fn('multiframe.pdb'))
    yield lambda: eq(t.n_frames, 2)
    yield lambda: eq(t.n_atoms, 22)