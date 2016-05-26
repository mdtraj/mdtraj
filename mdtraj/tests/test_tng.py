'''
Created on 05.12.2015

@author: marscher
'''
from mdtraj.formats.tng import TNGTrajectoryFile
from mdtraj.testing.testing import get_fn

test_fn = get_fn('tng_example.tng')

def test_len_with_ctx():
    with TNGTrajectoryFile(test_fn) as fh:
        res = len(fh)
        print (res)
        assert res == 10
        
def test_len():
    fh = TNGTrajectoryFile(test_fn)
    res = len(fh)
    print (res)
    assert res == 10
    fh.close()