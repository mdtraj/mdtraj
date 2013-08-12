###############################################################################
# Imports
###############################################################################

import os
import shutil
import tempfile
import subprocess
from distutils.spawn import find_executable

import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif

import numpy as np
import scipy.sparse

###############################################################################
# Globals
###############################################################################

HBondDocStringTester = DocStringFormatTester(md.geometry.hbond)
HAVE_DSSP = find_executable('mkdssp')
tmpdir = None

def setup():
    global tmpdir
    tmpdir = tempfile.mkdtemp()

def teardown():
    shutil.rmtree(tmpdir)

@skipif(not HAVE_DSSP, "This tests required mkdssp to be installed, from http://swift.cmbi.ru.nl/gv/dssp/")
def test_hbonds():
    t = md.load(get_fn('2EQQ.pdb'))[0]
    pdb = os.path.join(tmpdir, 'f.pdb')
    dssp = os.path.join(tmpdir, 'f.pdb.dssp')
    t.save(pdb)

    cmd = ['mkdssp', '-i', pdb, '-o', dssp]
    subprocess.check_output(' '.join(cmd), shell=True)
    energy = scipy.sparse.lil_matrix((t.n_atoms, t.n_atoms))

    # read the dssp N-H-->O column from the output file
    with open(dssp) as f:
        # skip the lines until the description of each residue's hbonds
        while not f.readline().startswith('  #  RESIDUE AA STRUCTURE'):
            continue

        for i, line in enumerate(f):
            line = line.rstrip()
            offset0, e0 = map(float, line[39:50].split(','))
            offset1, e1 = map(float, line[61:72].split(','))            
            energy[i, int(i+offset0)] = e0
            energy[i, int(i+offset1)] = e1
            #print e0, e1

    energy = energy.todense()
    energy[energy > -0.49] = 0
    print "\nOur code"
    print md.geometry.hbond.kabsch_sander(t)[0]
    print "\nFrom DSSP"
    print scipy.sparse.csr_matrix(energy.T)
