__author__ = 'CTK'
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
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


from __future__ import print_function
import numpy as np

import mdtraj as md
import mdtraj.geometry
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif, raises

"""
Reference values taken from gromacs 5.0.4 using the following commands on the
tip3_300K_1ATM.xtc trajectory.

echo -e "a O \n a H1 | a H2 \n q" | gmx make_ndx -f tip3p_300K_1ATM.pdb \
        -o tip3p_300K_1ATM.ndx

echo 3 3 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3_300K_1ATM_O-O_rdf.xvg
echo 3 4 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_O-H_rdf.xvg
echo 0 0 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_all-all_rdf.xvg

"""


def test_rdf_all_to_all():
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'),
                   top=get_fn('tip3p_300K_1ATM.pdb'))
    g_r, r = mdtraj.geometry.rdf.compute_rdf(traj, pair_names=None, r_range=[0, XXX], n_bins=XXX)
    g_r0, r0 = np.loadtxt(get_fn("tip3p_300K_1ATM_all-all_rdf.xvg"), )
    eq(r, r0)
    eq(g_r, g_r0)

def test_rdf_O_to_O():
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'),
                   top=get_fn('tip3p_300K_1ATM.pdb'))
    g_r, r = mdtraj.geometry.rdf.compute_rdf(traj, pair_names=('O', 'O'), r_range=[0, XXX], n_bins=XXX)
    g_r0, r0 = np.loadtxt(get_fn("tip3p_300K_1ATM_O-O_rdf.xvg"), )
    eq(r, r0)
    eq(g_r, g_r0)

def test_rdf_O_to_H():
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'),
                   top=get_fn('tip3p_300K_1ATM.pdb'))
    g_r, r = mdtraj.geometry.rdf.compute_rdf(traj, pair_names=('O', 'O'), r_range=[0, XXX], n_bins=XXX)
    g_r0, r0 = np.loadtxt(get_fn("tip3p_300K_1ATM_O-O_rdf.xvg"), )
    eq(r, r0)
    eq(g_r, g_r0)
