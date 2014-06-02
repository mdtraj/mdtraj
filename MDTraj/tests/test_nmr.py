from __future__ import print_function
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif
from mdtraj import nmr

def test_1():
    t = md.load(get_fn('2EQQ.pdb'))
    result = nmr.chemical_shifts_spartaplus(t)

    print(result)

def test_2_scalar_couplings():
    t = md.load(get_fn('frame0.h5'))  # This is Alanine dipeptide
    for model in ["Ruterjans1999", "Bax2007", "Bax1997"]:
        indices, J = md.compute_J3_HN_HA(t)
        eq(indices.shape, (1, 4))
        eq(J.shape, (501, 1))
        J = J.mean()
        assert abs(J - 6.06) <= 2.0, "Value is far from experimental value."
        # 6.06 [Hz] is the value from Baldwin PNAS 2006 Table 1.  
        # We expect the models to give something comparable to this
        # If it doesn't, something is fishy.
        # Typical ranges are between 1 and 9.
        # Obviously this isn't a perfect test, but it's still a useful sanity check.


def test_3_scalar_couplings():
    t = md.load(get_fn('1bpi.pdb'))
    for model in ["Ruterjans1999", "Bax2007", "Bax1997"]:
        indices, J = md.compute_J3_HN_HA(t)
