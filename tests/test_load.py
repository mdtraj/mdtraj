from mdtraj import load

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
