import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif

try:
    from simtk.openmm import app
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

@skipif(not HAVE_OPENMM)
def test_topology_to_openmm():
    topology = md.load(get_fn('1bpi.pdb')).topology

    mm = topology.to_openmm()
    assert isinstance(mm, app.Topology)

    topology2 = md.Topology.from_openmm(mm)

    eq(topology, topology2)
