import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif

try:
    from simtk.openmm import app
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False

@skipif(not HAVE_OPENMM)
def test_topology_openmm():
    topology = md.load(get_fn('1bpi.pdb')).topology

    mm = topology.to_openmm()
    assert isinstance(mm, app.Topology)

    topology2 = md.Topology.from_openmm(mm)

    eq(topology, topology2)


@skipif(not HAVE_PANDAS)
def test_topology_pandas():
    topology = md.load(get_fn('native.pdb')).topology
    
    df = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(df)

    eq(topology, topology2)


test_topology_pandas()
