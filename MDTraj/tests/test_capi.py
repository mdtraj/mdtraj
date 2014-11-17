import os
import mdtraj


def test_1():
    capi = mdtraj.capi()
    assert os.path.isdir(capi['lib_path'])
    assert os.path.isdir(capi['include_path'])
