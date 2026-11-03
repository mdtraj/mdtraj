import functools
import os
import tempfile
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.error import HTTPError

import pytest

from mdtraj import load, load_frame
from mdtraj.testing import eq


def test_load_single(get_fn):
    # Just check for any raised errors coming from loading a single file.
    load(get_fn("frame0.pdb"))


def test_load_single_top_none(get_fn):
    # Just check for any raised errors coming from loading a single file with top=None
    load(get_fn("frame0.pdb"), top=None)


def test_load_single_list(get_fn):
    # See if a single-element list of files is successfully loaded.
    load([get_fn("frame0.pdb")])


def test_load_many_list(get_fn):
    # See if a multi-element list of files is successfully loaded.
    single = load(get_fn("frame0.pdb"))
    double = load(2 * [get_fn("frame0.pdb")], discard_overlapping_frames=False)
    assert 2 * single.n_frames == double.n_frames


def test_load_atom_indices_multiple_files(get_fn):
    ref_t = load(get_fn("native.pdb"))
    t = load([get_fn("native.pdb")] * 2, atom_indices=[0])

    eq(t.topology, ref_t.topology.subset([0]))


@pytest.fixture(scope="module")
def data_url(get_fn):
    """Serve the test data directory over HTTP from a local thread.

    Loading from a URL is format-agnostic, so it is tested against a local
    server rather than a recorded RCSB cassette: no network, and any fixture
    file can be requested by name.
    """
    directory = os.path.dirname(get_fn("frame0.pdb"))
    handler = functools.partial(SimpleHTTPRequestHandler, directory=directory)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


@pytest.mark.parametrize("filename", ["4py5.cif", "frame0.gro", "frame0.h5", "frame0.dcd"])
def test_load_non_pdb_from_url(data_url, get_fn, filename):
    # Only the PDB reader can open a URL itself; every other loader takes a
    # path, so load() has to fetch the URL to a local file first (#2118).
    kwargs = {"top": get_fn("native.pdb")} if filename.endswith(".dcd") else {}
    from_url = load(f"{data_url}/{filename}", **kwargs)
    from_disk = load(get_fn(filename), **kwargs)
    eq(from_url.n_frames, from_disk.n_frames)
    eq(from_url.n_atoms, from_disk.n_atoms)
    eq(from_url.xyz, from_disk.xyz)


def test_load_frame_non_pdb_from_url(data_url, get_fn):
    from_url = load_frame(f"{data_url}/frame0.xtc", 3, top=get_fn("native.pdb"))
    from_disk = load_frame(get_fn("frame0.xtc"), 3, top=get_fn("native.pdb"))
    eq(from_url.xyz, from_disk.xyz)


def test_load_topology_from_url(data_url, get_fn):
    # A URL works for the topology argument as well as the trajectory.
    from_url = load(f"{data_url}/frame0.xtc", top=f"{data_url}/native.pdb")
    from_disk = load(get_fn("frame0.xtc"), top=get_fn("native.pdb"))
    eq(from_url.xyz, from_disk.xyz)
    eq(from_url.topology, from_disk.topology)


def test_load_url_list_mixed_with_paths(data_url, get_fn):
    double = load([f"{data_url}/frame0.pdb", get_fn("frame0.pdb")], discard_overlapping_frames=False)
    single = load(get_fn("frame0.pdb"))
    assert 2 * single.n_frames == double.n_frames


def test_load_from_url_cleans_up_temporary_files(data_url, tmp_path, monkeypatch):
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    load(f"{data_url}/frame0.gro")
    assert list(tmp_path.iterdir()) == []


def test_load_from_url_not_found(data_url):
    with pytest.raises(HTTPError):
        load(f"{data_url}/does-not-exist.gro")
