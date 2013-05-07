from mdtraj.testing import get_fn, eq
from mdtraj.trajectory import load


def test_0():
    t1 = load(get_fn('native2.xml'), top=get_fn('native2.pdb'))
    t2 = load(get_fn('native2.pdb'))

    t1.center_coordinates()
    t2.center_coordinates()

    yield lambda: eq(t1.xyz, t2.xyz)
    yield lambda: eq(t1.unitcell_vectors, t2.unitcell_vectors)
