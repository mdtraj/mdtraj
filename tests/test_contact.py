##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Christian Schwantes
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

from __future__ import print_function

import itertools

import mdtraj as md
import numpy as np
from mdtraj.testing import eq


def test_contact_0(get_fn):
    pdb = md.load(get_fn('bpti.pdb'))
    contacts = np.loadtxt(get_fn('contacts.dat')).astype(int)

    ca, ca_pairs = md.compute_contacts(pdb, contacts, scheme='ca')
    closest, closest_pairs = md.compute_contacts(pdb, contacts, scheme='closest')
    closest_heavy, closest_heavy_pairs = md.compute_contacts(pdb, contacts, scheme='closest-heavy')
    sidechain, sidechain_pairs = md.compute_contacts(pdb, contacts, scheme='sidechain')
    sidechain_heavy, sidechain_heavy_pairs = md.compute_contacts(pdb, contacts, scheme='sidechain-heavy')

    ref_ca = np.loadtxt(get_fn('cc_ca.dat'))
    ref_closest = np.loadtxt(get_fn('cc_closest.dat'))
    ref_closest_heavy = np.loadtxt(get_fn('cc_closest-heavy.dat'))
    ref_sidechain = np.loadtxt(get_fn('cc_sidechain.dat'))
    ref_sidechain_heavy = np.loadtxt(get_fn('cc_sidechain-heavy.dat'))

    eq(ref_ca, ca.flatten())
    eq(ref_closest, closest.flatten())
    eq(ref_closest_heavy, closest_heavy.flatten())
    eq(ref_sidechain, sidechain.flatten())
    eq(ref_sidechain_heavy, sidechain_heavy.flatten())
    eq(contacts, ca_pairs)
    eq(contacts, closest_pairs)
    eq(contacts, closest_heavy_pairs)
    eq(contacts, sidechain_pairs)
    eq(contacts, sidechain_heavy_pairs)


def test_contact_1(get_fn):
    pdb = md.load(get_fn('bpti.pdb'))
    dists, pairs = md.compute_contacts(pdb)
    for r0, r1 in pairs:
        # are these valid residue indices?
        pdb.topology.residue(r0)
        pdb.topology.residue(r1)

        assert not (abs(r0 - r1) < 3)

    maps = md.geometry.squareform(dists, pairs)
    for i, (r0, r1) in enumerate(pairs):
        for t in range(pdb.n_frames):
            eq(maps[t, r0, r1], dists[t, i])


def test_contact_2(get_fn):
    pdb = md.load(get_fn('1vii_sustiva_water.pdb'))
    dists, pairs = md.compute_contacts(pdb, scheme='closest')
    for r0, r1 in pairs:
        assert pdb.topology.residue(r0).name != 'HOH'
        assert pdb.topology.residue(r1).name != 'HOH'

    # spot check one of the pairs
    r0, r1 = pairs[10]
    atoms_r0 = [a.index for a in pdb.topology.residue(r0).atoms]
    atoms_r1 = [a.index for a in pdb.topology.residue(r1).atoms]

    atomdist = md.compute_distances(pdb, list(itertools.product(atoms_r0, atoms_r1)))

    np.testing.assert_array_equal(dists[:, 10], np.min(atomdist, axis=1))

    maps = md.geometry.squareform(dists, pairs)
    for i, (r0, r1) in enumerate(pairs):
        for t in range(pdb.n_frames):
            eq(maps[t, r0, r1], dists[t, i])


def test_contact_3(get_fn):
    pdb = md.load(get_fn('bpti.pdb'))
    beta = 20
    dists, pairs = md.compute_contacts(pdb, soft_min=True, soft_min_beta=beta)

    maps = md.geometry.squareform(dists, pairs)
    for i, (r0, r1) in enumerate(pairs):
        for t in range(pdb.n_frames):
            assert np.allclose(beta / np.log(np.sum(np.exp(beta / maps[t, r0, r1]))), dists[t, i])
            
def test_contact_4(get_fn):
    pdb = md.load(get_fn('1am7_protein.pdb'))   # protonated and including at least one glycine residue (which has no heavy atoms in its sidechain)
    contacts = md.compute_contacts(pdb, contacts='all', scheme='sidechain-heavy')   # test passes if this doesn't raise an exception
