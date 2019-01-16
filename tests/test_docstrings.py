import re
import mdtraj
from mdtraj.testing.docstrings import docstring_verifiers, import_all_modules
from mdtraj.testing.docstrings import import_all_modules
import pytest

SKIP_MODULES = [
    r'mdtraj\.utils\.external',
    r'mdtraj\.utils\.six',
    r'mdtraj\.utils\.unit\.unit_math',
    r'mdtraj\.utils\.unit\.baseunit',
    r'mdtraj\.utils\.unit\.prefix',
    r'mdtraj\.utils\.unit\.unit',
    r'mdtraj\.utils\.unit\.quantity',
    r'mdtraj\.utils\.unit\.mymatrix',
    r'mdtraj\.formats\.lh5',
    r'mdtraj\.formats\.hdf5',
    r'mdtraj\.formats\.pdb\.pdbstructure',
    r'mdtraj\.scripts',
    r'mdtraj\.testing\.docscrape',
    r'mdtraj\.io',
]

MODULES = import_all_modules(mdtraj)


@pytest.fixture(params=MODULES, ids=lambda x: x.__name__)
def module(request):
    mod = request.param
    if any(re.search(s, mod.__name__) for s in SKIP_MODULES):
        pytest.skip("Skip checking certain modules' docstrings")

    return mod


def test_docstrings(module):
    for test in docstring_verifiers(module):
        test()
