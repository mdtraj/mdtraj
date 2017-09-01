import re
import mdtraj
from mdtraj.testing.docstrings import docstring_verifiers, import_all_modules
from mdtraj.testing.docstrings import import_all_modules
import pytest

SKIP_MODULES = [
    'mdtraj\.utils\.external',
    'mdtraj\.utils\.six',
    'mdtraj\.utils\.unit\.unit_math',
    'mdtraj\.utils\.unit\.baseunit',
    'mdtraj\.utils\.unit\.prefix',
    'mdtraj\.utils\.unit\.unit',
    'mdtraj\.utils\.unit\.quantity',
    'mdtraj\.utils\.unit\.mymatrix',
    'mdtraj\.formats\.lh5',
    'mdtraj\.formats\.hdf5',
    'mdtraj\.formats\.pdb\.pdbstructure',
    'mdtraj\.scripts',
    'mdtraj\.testing\.docscrape',
    'mdtraj\.io',
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
