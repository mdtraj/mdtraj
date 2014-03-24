import mdtraj as md
from mdtraj.testing import get_fn
from mdtraj import nmr

def test_1():
    t = md.load(get_fn('2EQQ.pdb'))
    result = nmr.chemical_shifts_spartaplus(t)

    print result
