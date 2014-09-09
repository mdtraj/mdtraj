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
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif
from mdtraj.formats import psf
from mdtraj.utils import enter_temp_directory
from subprocess import check_output
from distutils.spawn import find_executable

doc = DocStringFormatTester(psf)
VMD = find_executable('vmd')
VMD_MSG = 'This test requires the VMD executable: http://www.ks.uiuc.edu/Research/vmd/'

def test_load_psf():
    top = psf.load_psf(get_fn('ala_ala_ala.psf'))
    ref_top = md.load(get_fn('ala_ala_ala.pdb')).topology
    eq(top, ref_top)
    # Test the XPLOR psf format parsing
    top2 = psf.load_psf(get_fn('ala_ala_ala.xpsf'))
    eq(top2, ref_top)


def test_load_mdcrd_with_psf():
    traj = md.load(get_fn('ala_ala_ala.mdcrd'), top=get_fn('ala_ala_ala.psf'))
    ref_traj = md.load(get_fn('ala_ala_ala.pdb'))
    eq(traj.topology, ref_traj.topology)
    eq(traj.xyz, ref_traj.xyz)


@skipif(not VMD, VMD_MSG)
def test_vmd_psf():
    yield lambda: _test_against_vmd(get_fn('1vii.pdb'))
    yield lambda: _test_against_vmd(get_fn('2EQQ.pdb'))


def _test_against_vmd(pdb):
    top = os.path.join(os.path.dirname(VMD), '..', 'lib', 'plugins', 'noarch',
                       'tcl', 'readcharmmtop1.0', 'top_all27_prot_lipid_na.inp')
    assert os.path.exists(top)

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
    ''' % {'top': top, 'pdb' : pdb}

    with enter_temp_directory():
        with open('script.tcl', 'w') as f:
            f.write(TEMPLATE)
        check_output([VMD, '-e', 'script.tcl', '-dispdev', 'none'])
        out_pdb = md.load('out.pdb')
        out_psf = md.load_psf('out.psf')

        # make sure the two topologies are equal
        eq(out_pdb.top, out_psf)
