import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq
from mdtraj.rmsd.irmsd import rmsd_cache, RMSDCache, align_array
from mdtraj.geometry.alignment import rmsd_qcp
import matplotlib.pyplot as pp

def test_axis_major():
    t = md.load(get_fn('traj.h5'))
    pt = rmsd_cache(t, major='axis')
    calculated = pt.rmsds_to_reference(pt, 10)


    reference = np.zeros(t.n_frames)
    for i in range(t.n_frames):
        reference[i] = rmsd_qcp(t.xyz[10], t.xyz[i])

    pp.clf()
    pp.scatter(reference, calculated)
    pp.xlabel('reference')
    pp.ylabel('calculated')

    pp.savefig('rmsd_axis.png')
    eq(calculated, reference, decimal=4)


def test_atom_major():
    t = md.load(get_fn('traj.h5'))
    pt = rmsd_cache(t, major='atom')
    calculated = pt.rmsds_to_reference(pt, 10)

    reference = np.zeros(t.n_frames)
    for i in range(t.n_frames):
        reference[i] = rmsd_qcp(t.xyz[10], t.xyz[i])

    pp.clf()
    pp.scatter(reference, calculated)
    pp.xlabel('reference')
    pp.ylabel('calculated')

    pp.savefig('rmsd_axis.png')
    eq(calculated, reference, decimal=4)


def test_rmsd_to_self():
    n_atoms = 16
    conf_atom_1 = np.array(np.random.randn(1, n_atoms, 3), dtype=np.float32)
    conf_atom_2 = np.copy(conf_atom_1)

    r_atom_1 = RMSDCache(align_array(conf_atom_1, 'atom'), major='atom', n_atoms=n_atoms)
    r_atom_2 = RMSDCache(align_array(conf_atom_2, 'atom'), major='atom', n_atoms=n_atoms)

    yield lambda: eq(float(r_atom_1.rmsds_to_reference(r_atom_2, 0)[0]), 0.0, decimal=4)

    conf_axis_1 = np.copy(conf_atom_1[0].T.reshape(1, 3, n_atoms))
    conf_axis_2 = np.copy(conf_axis_1)

    r_axis_1 = RMSDCache(align_array(conf_axis_1, 'axis'), major='axis', n_atoms=n_atoms)
    r_axis_2 = RMSDCache(align_array(conf_axis_2, 'axis'), major='axis', n_atoms=n_atoms)

    yield lambda: eq(float(r_axis_1.rmsds_to_reference(r_axis_2, 0)[0]), 0.0, decimal=4)
