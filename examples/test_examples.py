"""
Execute each notebook as a test, reporting an error if any cell throws an exception.
Adapted from https://gist.github.com/minrk/2620876.
"""
from __future__ import print_function
import os
import sys

import nbformat
import pytest
from jupyter_client import KernelManager
from six.moves.queue import Empty

FLAKEY_LIST = ['centroids.ipynb', 'native-contact.ipynb', 'hbonds.ipynb']
TIMEOUT = 60  # seconds

test_dir = os.path.dirname(os.path.abspath(__file__))
examples = [pytest.param(fn, marks=pytest.mark.flaky) if fn in FLAKEY_LIST else fn
            for fn in os.listdir(test_dir) if fn.endswith('.ipynb')]


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
        try:
            # wait for finish, w/ timeout
            reply = shell.get_msg(timeout=TIMEOUT)['content']
        except Empty:
            raise Exception(
                'Timeout (%.1f) when executing the following %s cell: "%s"' %
                (TIMEOUT, cell.cell_type, cell.source.strip()))
        if reply['status'] == 'error':
            failures += 1
            print("\nFAILURE:", file=sys.stderr)
            print('\n'.join(reply['traceback']), file=sys.stderr)
            print(file=sys.stderr)

    kc.stop_channels()
    km.shutdown_kernel()
    del km
    if failures > 0:
        raise Exception()
