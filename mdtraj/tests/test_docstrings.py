import re
import mdtraj
from mdtraj.testing.docstrings import docstring_verifiers
from mdtraj.testing.docstrings import import_all_modules

SKIP_MODULES = [
    'mdtraj\.utils\.external',
    'mdtraj\.utils\.six',
    'mdtraj\.utils\.unit\.unit_math',
    'mdtraj\.utils\.unit\.baseunit',
    'mdtraj\.utils\.unit\.prefix',
]

SKIP_FUNCTIONS = [
]


def test_all_docstrings():
    for module in import_all_modules(mdtraj):
        if any(re.search(s, module.__name__) for s in SKIP_MODULES):
            continue

        for test in docstring_verifiers(module):
            if test.fname in SKIP_FUNCTIONS:
                continue
            yield test
