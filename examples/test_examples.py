"""
Execute each notebook as a test, reporting an error if any cell throws an exception.
Adapted from https://gist.github.com/minrk/2620876.
"""
import os
import sys

import nbformat
from jupyter_client import KernelManager

def test_examples():
    for f in os.listdir('.'):
        if f.endswith('.ipynb'):
            yield check_one_notebook, f

def check_one_notebook(filename):
    with open(filename) as f:
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
        reply = shell.get_msg(timeout=20)['content']
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