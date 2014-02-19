from mdtraj.testing import get_fn, eq
from mdtraj.pdb import element
import mdtraj as md

def test_element_0():
    t = md.load(get_fn('bpti.pdb'))

    a = t.top.atom(15)
    H = element.Element.getBySymbol('H')

    eq(a.element, element.hydrogen)

