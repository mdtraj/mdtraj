"""
Execute each notebook as a test, reporting an error if any cell throws an exception.
Adapted from https://gist.github.com/minrk/2620876.
"""

import os
import shutil
import socket
import sys
from queue import Empty

import pytest

try:
    import nbformat
    from jupyter_client import KernelManager
except Exception:
    pytest.skip("Skipping no nbformat/jupyter", allow_module_level=True)

FLAKEY_LIST = ["centroids.ipynb", "native-contact.ipynb", "hbonds.ipynb"]
SPARTA_PLUS = ["sparta+", "SPARTA+", "SPARTA+.linux"]
TIMEOUT = 60  # seconds

test_dir = os.path.dirname(os.path.abspath(__file__))
examples = [
    pytest.param(fn, marks=pytest.mark.flaky) if fn in FLAKEY_LIST else fn
    for fn in os.listdir(test_dir)
    if fn.endswith(".ipynb")
]


def is_network_connected():
    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection(("1.1.1.1", 53))
        return True
    except OSError:
        pass
    return False


def find_executable(names):
    for possible in names:
        result = shutil.which(possible)
        if result is not None:
            return result
    return None


@pytest.fixture(params=examples)
def example_fn(request):
    if "openmm" in request.param:
        try:
            from openmm import app  # noqa: F401
        except ImportError:
            pytest.skip(
                f"Openmm required for example notebook `{request.param}`",
            )

    if "nmr" in request.param:
        if find_executable(SPARTA_PLUS) is None:
            pytest.skip(
                f"Sparta+ not found for example notebook `{request.param}`",
            )

    if not is_network_connected():
        if any(x in request.param for x in ("native-contact", "hbonds")):
            pytest.skip("Network access required")

    cwd = os.path.abspath(".")
    os.chdir(test_dir)
    yield request.param
    os.chdir(cwd)


def test_example_notebook(example_fn):
    with open(example_fn) as f:
        nb = nbformat.reads(f.read(), nbformat.NO_CONVERT)
    run_notebook(nb)


def run_notebook(nb):
    km = KernelManager()
    km.start_kernel(stderr=open(os.devnull, "w"))
    kc = km.client()
    kc.start_channels()
    # simple ping:
    kc.execute("pass")
    kc.get_shell_msg()

    failures = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        kc.execute(cell.source)
        try:
            # wait for finish, w/ timeout
            reply = kc.get_shell_msg(timeout=TIMEOUT)["content"]
        except Empty:
            raise Exception(
                f'Timeout ({TIMEOUT:.1f}) when executing the following {cell.cell_type} cell: "{cell.source.strip()}"',
            )
        if reply["status"] == "error":
            failures += 1
            print("\nFAILURE:", file=sys.stderr)
            print("\n".join(reply["traceback"]), file=sys.stderr)
            print(file=sys.stderr)

    kc.stop_channels()
    km.shutdown_kernel()
    del km
    if failures > 0:
        raise Exception()
