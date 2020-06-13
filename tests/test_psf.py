##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Jason Swails
# Contributors: Robert McGibbon
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import os
import subprocess
from distutils.spawn import find_executable

import mdtraj as md
import pytest
from mdtraj.formats import psf
from mdtraj.testing import eq
from mdtraj.utils import enter_temp_directory

VMD = find_executable('vmd')
needs_vmd = pytest.mark.skipif(
    not find_executable('vmd'),
    reason='This test requires the VMD executable: http://www.ks.uiuc.edu/Research/vmd/'
)


def test_load_psf(get_fn):
    top = psf.load_psf(get_fn('ala_ala_ala.psf'))
    ref_top = md.load(get_fn('ala_ala_ala.pdb')).topology
    eq(top, ref_top)
    # Test the XPLOR psf format parsing
    top2 = psf.load_psf(get_fn('ala_ala_ala.xpsf'))
    eq(top2, ref_top)
    # Test segment_names are loaded properly
    assert next(top.residues).segment_id == "AAL", "Segment id is not being assigned correctly for ala_ala_ala.psf"
    assert next(top2.residues).segment_id == "AAL", "Segment id is not being assigned correctly for ala_ala_ala.xpsf"


def test_multichain_psf(get_fn):
    top = psf.load_psf(get_fn('3pqr_memb.psf'))
    # Check that each segment was turned into a chain
    eq(top.n_chains, 11)
    chain_lengths = [5162, 185, 97, 28, 24, 24, 45, 35742, 72822, 75, 73]
    for i, n in enumerate(chain_lengths):
        eq(top.chain(i).n_atoms, n)
    # There are some non-sequentially numbered residues across chain
    # boundaries... make sure resSeq is properly taken from the PSF
    eq(top.chain(0).residue(0).resSeq, 1)
    eq(top.chain(0).residue(-1).resSeq, 326)
    eq(top.chain(1).residue(0).resSeq, 340)
    eq(top.chain(1).residue(-1).resSeq, 350)
    eq(top.chain(2).residue(0).resSeq, 1)
    eq(top.chain(2).residue(-1).resSeq, 4)
    eq(top.chain(3).residue(0).resSeq, 1)
    eq(top.chain(3).residue(-1).resSeq, 1)
    eq(top.chain(4).residue(0).resSeq, 1)
    eq(top.chain(4).residue(-1).resSeq, 1)
    eq(top.chain(5).residue(0).resSeq, 1)
    eq(top.chain(5).residue(-1).resSeq, 1)
    eq(top.chain(6).residue(0).resSeq, 1)
    eq(top.chain(6).residue(-1).resSeq, 2)
    eq(top.chain(7).residue(0).resSeq, 1)
    eq(top.chain(7).residue(-1).resSeq, 276)
    eq(top.chain(8).residue(0).resSeq, 1)
    eq(top.chain(8).residue(-1).resSeq, 24274)
    eq(top.chain(9).residue(0).resSeq, 1)
    eq(top.chain(9).residue(-1).resSeq, 75)
    eq(top.chain(10).residue(0).resSeq, 1)
    eq(top.chain(10).residue(-1).resSeq, 73)


def test_load_mdcrd_with_psf(get_fn):
    traj = md.load(get_fn('ala_ala_ala.mdcrd'), top=get_fn('ala_ala_ala.psf'))
    ref_traj = md.load(get_fn('ala_ala_ala.pdb'))
    eq(traj.topology, ref_traj.topology)
    eq(traj.xyz, ref_traj.xyz)


@needs_vmd
@pytest.mark.parametrize('pdb', ['1vii.pdb', '2EQQ.pdb'])
def test_against_vmd(pdb, get_fn):
    pdb = get_fn(pdb)
    # this is probably not cross-platform compatible. I assume that the exact
    # path to this CHARMM topology that is included with VMD depends on
    # the install mechanism, especially for bundled mac or windows installers
    VMD_ROOT = os.path.join(os.path.dirname(os.path.realpath(VMD)), '..')
    top_paths = [os.path.join(r, f) for (r, _, fs) in os.walk(VMD_ROOT) for f in fs
                 if 'top_all27_prot_lipid_na.inp' in f]
    assert len(top_paths) >= 0
    top = os.path.abspath(top_paths[0]).replace(" ", "\\ ")

    TEMPLATE = '''
package require psfgen
topology %(top)s
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment U {pdb %(pdb)s}
coordpdb %(pdb)s U
guesscoord
writepdb out.pdb
writepsf out.psf
exit
    ''' % {'top': top, 'pdb': pdb}

    with enter_temp_directory():
        with open('script.tcl', 'w') as f:
            f.write(TEMPLATE)
        subprocess.check_call([VMD, '-startup', 'script.tcl', '-dispdev', 'none'])
        out_pdb = md.load('out.pdb')
        out_psf = md.load_psf('out.psf')

        # make sure the two topologies are equal
        eq(out_pdb.top, out_psf)
