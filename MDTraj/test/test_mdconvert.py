import os
import tempfile
import shutil
import numpy as np

import mdtraj as md
from mdtraj.testing import skipif, get_fn, eq

try:
    from scripttest import TestFileEnvironment
    HAVE_SCRIPTTEST = True
except ImportError:
    HAVE_SCRIPTTEST = False

test_dir = './mdconvert-test-output'
staging_dir = tempfile.mkdtemp()
def teardown_module(module):
    shutil.rmtree(staging_dir)
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)

def setup_module():
    global TRAJ
    
    xyz = np.arange(150, dtype=np.float32).reshape(10,5,3)
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue('ALA', chain)
    topology.add_atom('CA', md.pdb.element.carbon, residue)
    topology.add_atom('HG1', md.pdb.element.hydrogen, residue)
    topology.add_atom('SG', md.pdb.element.sulfur, residue)
    topology.add_atom('OD1', md.pdb.element.oxygen, residue)
    topology.add_atom('NE', md.pdb.element.nitrogen, residue)

    time = np.arange(10)**2
    unitcell_lengths = np.array([[1.1,1.2,1.3]] * 10)
    unitcell_angles = np.array([[80, 90, 110]] * 10)
    
    TRAJ = md.Trajectory(xyz, topology=topology, time=time,
                         unitcell_lengths=unitcell_lengths,
                         unitcell_angles=unitcell_angles)


@skipif(not HAVE_SCRIPTTEST, 'need "scripttest" module to test mdconvert (pip install scripttest)')
def test_mdconvert_0():
    """ensure that the xyz coordinates are preserved by a trip
       from python -> save in format X -> mdconvert to format Y -> python
    """
    env = TestFileEnvironment(test_dir)
    
    # save one copy of traj for use as a topology file
    topology_fn = os.path.join(staging_dir, 'topology.pdb')
    TRAJ.save(topology_fn)
    
    fns = ['traj.xtc', 'traj.dcd', 'traj.binpos', 'traj.trr', 'traj.nc', 'traj.pdb', 'traj.h5']
    
    for fn in fns:
        path = os.path.join(staging_dir, fn)
        TRAJ.save(path)
        
        for fn2 in filter(lambda e: e != fn, fns):
            if os.path.splitext(fn2)[1] in ['.pdb', '.h5']:
                generate = lambda : env.run('mdconvert', path, '-o', fn2, '-c 3', '-t', topology_fn, expect_stderr=True)
            else:
                generate = lambda : env.run('mdconvert', path, '-o', fn2, '-c 3', expect_stderr=True)
            
            generate.description = 'Converting %s -> %s' % (fn, fn2)
            yield generate
                
            # ensure that the xyz coordinates are preserved by a trip
            # from python -> save in format X -> mdconvert to format Y -> python
            if os.path.splitext(fn2)[1] in ['.pdb', '.h5']:
                check = lambda : eq(md.load(os.path.join(test_dir, fn2)).xyz, TRAJ.xyz)
            else:
                check = lambda : eq(md.load(os.path.join(test_dir, fn2), top=TRAJ.topology).xyz,
                    TRAJ.xyz)
            check.description = 'Checking %s -> %s' % (fn, fn2)
            yield check
             
            env.run('rm', '-f', fn2)        
        os.unlink(path)

@skipif(not HAVE_SCRIPTTEST, 'need "scripttest" module to test mdconvert (pip install scripttest)')
def test_mdconvert_1():
   """Test that atom_indices are converted correctly
   """
   pass