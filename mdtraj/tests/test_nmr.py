from __future__ import print_function
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif
from mdtraj.nmr.shift_wrappers import find_executable, SPARTA_PLUS, PPM, SHIFTX2
from mdtraj.utils import six

@skipif(not find_executable(SPARTA_PLUS), 'SPARTA+ binary not found')
def test_spartaplus():
    t = md.load(get_fn('2EQQ.pdb'))
    result = md.chemical_shifts_spartaplus(t)

    print(result)
    eq(result.shape[1], 20)  # 2EQQ is NMR structure with 20 frames
    eq(float(result.ix[(1, "HA")][0]), 4.378, decimal=4)
    # 4.378 taken from first entry in pred.tab, which looks like the following:
    #  1    E   HA     0.000     4.378     4.350     0.047     0.000     0.291


@skipif(not find_executable(PPM), 'PPM binary not found')
def test_ppm():
    t = md.load(get_fn('2EQQ.pdb'))
    result = md.chemical_shifts_ppm(t)

    print(result)
    
    eq(result.shape[1], 20)  # 2EQQ is NMR structure with 20 frames
    eq(float(result.ix[(2, "CA")][0]), 53.004, decimal=4)
    # taken from first entry in bb_details.dat, which looks like the following:
    #       2     ASN   CA    999.000  53.004  51.168  51.802  53.081  54.098  52.820  52.379  51.856  53.034  52.754  54.134  54.222  51.210  52.207  50.824  54.459  53.605  54.211  53.688  52.344  53.004  51.168  51.802  53.081  54.098  52.820  52.379  51.856  53.034  52.754  54.134  54.222  51.210  52.207  50.824  54.459  53.605  54.211  53.688  52.344


@skipif(not (find_executable(SHIFTX2) and six.PY2), 'SHIFTX2 binary not found or python2.X not found.')
def test_shiftx2():
    t = md.load(get_fn('2EQQ.pdb'))
    result = md.chemical_shifts_shiftx2(t)

    print(result)
    
    eq(result.shape[1], 20)  # 2EQQ is NMR structure with 20 frames
    eq(float(result.ix[(1, "C")][0]), 175.2570, decimal=4)
    # taken from first entry in trj.pdb.cs, which looks like the following:
    #  NUM,RES,ATOMNAME,SHIFT
    #  1,E,C,175.2570


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
