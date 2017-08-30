import os
import sys
import mdtraj


def test_1():
    capi = mdtraj.capi()
    assert os.path.isdir(capi['lib_dir'])
    assert os.path.isdir(capi['include_dir'])
    if sys.platform == 'win32':
        assert os.path.isfile(os.path.join(capi['lib_dir'], 'libtheobald.lib'))
    else:
        assert os.path.isfile(os.path.join(capi['lib_dir'], 'libtheobald.a'))
