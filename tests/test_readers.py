from mdtraj import xtc
from mdtraj import dcd
import numpy.testing as npt

def test1():
    filename_dcd = '/Users/rmcgibbo/local/msmbuilder/Tutorial/DCD/RUN00/frame0.dcd'
    filename_xtc = '/Users/rmcgibbo/local/msmbuilder/Tutorial/XTC/RUN00/frame0.xtc'

    c = dcd.read_xyz(filename_dcd)
    a = xtc.read(filename_xtc, chunk=1)[0]
    b = xtc.read(filename_xtc, chunk=2)[0]

    
    npt.assert_array_equal(a,b)
    npt.assert_array_almost_equal(a, 0.1*c)
