"""
Execute each notebook as a test, reporting an error if any cell throws an exception.
Adapted from https://gist.github.com/minrk/2620876.
"""

import os

import nbformat
import pytest
from jupyter_client import KernelManager

test_dir = os.path.dirname(os.path.abspath(__file__))
examples = [fn for fn in os.listdir(test_dir) if fn.endswith('.ipynb')]


@pytest.fixture(params=examples)
def example_fn(request):
    if 'openmm' in request.param:
        try:
            from simtk.openmm import app
        except ImportError:
            pytest.skip("Openmm required for example notebook `{}`".format(request.param))

    cwd = os.path.abspath('.')
    os.chdir(test_dir)
    yield request.param
    os.chdir(cwd)


def test_example_notebook(example_fn):
    with open(example_fn) as f:
        nb = nbformat.reads(f.read(), nbformat.NO_CONVERT)
    run_notebook(nb)


def run_notebook(nb):
    km = KernelManager()
    km.start_kernel(stderr=open(os.devnull, 'w'))
    kc = km.client()
    kc.start_channels()
    shell = kc.shell_channel
    # simple ping:
    kc.execute("pass")
    shell.get_msg()

    failures = 0
    for cell in nb.cells:
        if cell.cell_type != 'code':
            continue
        kc.execute(cell.source)
        # wait for finish, maximum 20s
        reply = shell.get_msg(timeout=60)['content']
        if reply['status'] == 'error':
            failures += 1
            print("\nFAILURE:")
            print('\n'.join(reply['traceback']))
            print()

    kc.stop_channels()
    km.shutdown_kernel()
    del km
    if failures > 0:
        raise Exception()
