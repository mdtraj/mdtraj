import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn
from mdtraj.rmsd.irmsd import rmsd_cache
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
    
    pp.savefig('rmsd_atom.png')

