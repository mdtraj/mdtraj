from mdtraj import load
from mdtraj.testing import eq

def test_load_single(get_fn):
    # Just check for any raised errors coming from loading a single file.
    load(get_fn('frame0.pdb'))


def test_load_single_list(get_fn):
    # See if a single-element list of files is successfully loaded.
    load([get_fn('frame0.pdb')])


def test_load_many_list(get_fn):
    # See if a multi-element list of files is successfully loaded.
    single = load(get_fn('frame0.pdb'))
    double = load(2 * [get_fn('frame0.pdb')], discard_overlapping_frames=False)
    assert 2 * single.n_frames == double.n_frames

def test_load_atom_indices_multiple_files(get_fn):
    ref_t = load(get_fn('native.pdb'))
    t = load([get_fn('native.pdb')]*2, atom_indices=[0])

    eq(t.topology, ref_t.topology.subset([0]))
