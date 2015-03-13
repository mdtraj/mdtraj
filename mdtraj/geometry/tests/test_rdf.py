##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christoph Klein
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
from mdtraj.testing import get_fn, eq, skipif

"""
Reference values taken from gromacs 5.0.4 using the following commands on the
`tip3_300K_1ATM.xtc` trajectory.

echo -e "a O \n a H1 | a H2 \n q" | gmx make_ndx -f tip3p_300K_1ATM.pdb \
        -o tip3p_300K_1ATM.ndx

echo 3 3 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3_300K_1ATM_O-O_rdf.xvg
echo 3 4 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_O-H_rdf.xvg
echo 0 0 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_all-all_rdf.xvg

"""
TRAJ = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
pairs = TRAJ.top.select_pairs('all', 'all')
R, RDF_ALL = mdtraj.geometry.rdf.compute_rdf(TRAJ, pairs)
pairs = TRAJ.top.select_pairs('name O', 'name O')
_, RDF_O_O = mdtraj.geometry.rdf.compute_rdf(TRAJ, pairs)
pairs = TRAJ.top.select_pairs('name O', "name =~ 'H.*'")
_, RDF_O_H = mdtraj.geometry.rdf.compute_rdf(TRAJ, pairs)


def test_rdf_norm():
    """Check if the RDF's tail is normalized to ~1.0 """
    assert eq(np.ones(20), RDF_ALL[-20:], decimal=1)
    assert eq(np.ones(20), RDF_O_O[-20:], decimal=1)
    assert eq(np.ones(20), RDF_O_H[-20:], decimal=1)


@skipif(True, 'Binning does not match up with gromacs currently.')
def test_compare_rdf_to_gromacs():
    compare_gromacs_xvg('tip3p_300K_1ATM_all-all_rdf.xvg', R, RDF_ALL)
    compare_gromacs_xvg('tip3p_300K_1ATM_O-O_rdf.xvg', R, RDF_O_O)
    compare_gromacs_xvg('tip3p_300K_1ATM_O-H_rdf.xvg', R, RDF_O_H)


def compare_gromacs_xvg(ref_filename, r, g_r):
    """Compare RDFs to a reference file produced by gromacs. """
    data = np.loadtxt(get_fn(ref_filename), skiprows=16)
    r0 = data[:, 0]
    g_r0 = data[:, 1]
    eq(r, r0, decimal=2)
    eq(g_r, g_r0, decimal=2)
