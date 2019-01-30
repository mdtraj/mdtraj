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
from mdtraj.testing import eq
import pytest

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


def test_rdf_norm(get_fn):
    # Check if the RDF's tail is normalized to ~1.0
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))

    pairs = traj.top.select_pairs('all', 'all')
    _, rdf_all = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_all[-20:], decimal=1)

    pairs = traj.top.select_pairs('name O', 'name O')
    _, rdf_O_O = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_O_O[-20:], decimal=1)

    pairs = traj.top.select_pairs('name O', "name =~ 'H.*'")
    _, rdf_O_H = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_O_H[-20:], decimal=1)


@pytest.mark.skip('Binning does not match up with gromacs currently.')
def test_compare_rdf_to_gromacs(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    pairs = traj.top.select_pairs('all', 'all')
    R, rdf_all = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_all-all_rdf.xvg'), R, rdf_all)

    pairs = traj.top.select_pairs('name O', 'name O')
    _, rdf_O_O = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_O-O_rdf.xvg'), R, rdf_O_O)

    pairs = traj.top.select_pairs('name O', "name =~ 'H.*'")
    _, rdf_O_H = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_O-H_rdf.xvg'), R, rdf_O_H)


def compare_gromacs_xvg(ref_filename, r, g_r):
    """Compare RDFs to a reference file produced by gromacs. """
    data = np.loadtxt(ref_filename, skiprows=16)
    r0 = data[:, 0]
    g_r0 = data[:, 1]
    eq(r, r0, decimal=2)
    eq(g_r, g_r0, decimal=2)

def test_rdf_t(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))[:4]
    ref_pairs = traj.top.select_pairs('all', 'all')
    pairs = ref_pairs
    times = np.vstack([np.zeros(np.shape(pairs)[0]), np.zeros(np.shape(pairs)[0])]).T

    for n in range(1, traj.n_frames):
        times = np.vstack([times, np.vstack([n*np.ones(np.shape(ref_pairs)[0]), n*np.zeros(np.shape(ref_pairs)[0])]).T])
        pairs = np.vstack([pairs, ref_pairs])
    r, g_r = md.compute_rdf_t(traj, pairs, times, n_bins=1000, periodic=True, opt=False)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogy(r, g_r)
    plt.ylim((0.01, 100))
    plt.savefig('tmp.pdf')
    assert False
