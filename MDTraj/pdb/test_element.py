from mdtraj.testing import get_fn, eq
from mdtraj.pdb import element
import mdtraj as md

def test_it():
    t = md.load(get_fn('bpti.pdb'))

    a = t.top.atom(15)

    print a.element.symbol, a.element.name, a.element.mass
    print element.hydrogen.symbol, element.hydrogen.name, element.hydrogen.mass
    H = element.Element.getBySymbol('H')
    print H.symbol, H.name, H.mass

    print a.element == element.hydrogen

    print element.Element.getBySymbol('N').mass

    eq(a.element, element.hydrogen)

