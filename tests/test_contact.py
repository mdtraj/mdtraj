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


import itertools

import numpy as np
import pytest

import mdtraj as md
from mdtraj.testing import eq


def test_contact_0(get_fn):
    pdb = md.load(get_fn("bpti.pdb"))
    contacts = np.loadtxt(get_fn("contacts.dat")).astype(int)

    ca, ca_pairs = md.compute_contacts(pdb, contacts, scheme="ca")
    closest, closest_pairs = md.compute_contacts(pdb, contacts, scheme="closest")
    closest_heavy, closest_heavy_pairs = md.compute_contacts(
        pdb,
        contacts,
        scheme="closest-heavy",
    )
    sidechain, sidechain_pairs = md.compute_contacts(pdb, contacts, scheme="sidechain")
    with pytest.warns(UserWarning, match="selected topology includes at least one glycine"):
        sidechain_heavy, sidechain_heavy_pairs = md.compute_contacts(
            pdb,
            contacts,
            scheme="sidechain-heavy",
        )

    ref_ca = np.loadtxt(get_fn("cc_ca.dat"))
    ref_closest = np.loadtxt(get_fn("cc_closest.dat"))
    ref_closest_heavy = np.loadtxt(get_fn("cc_closest-heavy.dat"))
    ref_sidechain = np.loadtxt(get_fn("cc_sidechain.dat"))
    ref_sidechain_heavy = np.loadtxt(get_fn("cc_sidechain-heavy.dat"))

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
    pdb = md.load(get_fn("bpti.pdb"))
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
    pdb = md.load(get_fn("1vii_sustiva_water.pdb"))
    dists, pairs = md.compute_contacts(pdb, scheme="closest")
    for r0, r1 in pairs:
        assert pdb.topology.residue(r0).name != "HOH"
        assert pdb.topology.residue(r1).name != "HOH"

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
    pdb = md.load(get_fn("bpti.pdb"))
    beta = 20
    dists, pairs = md.compute_contacts(pdb, soft_min=True, soft_min_beta=beta)
    for i, (r0, r1) in enumerate(pairs):
        atom_pairs = list(
            itertools.product(
                [a.index for a in pdb.top.residue(r0).atoms if a.element != md.element.hydrogen],
                [a.index for a in pdb.top.residue(r1).atoms if a.element != md.element.hydrogen],
            ),
        )
        # Compute the defining expression independently in double precision.
        atom_distances = md.compute_distances(pdb, atom_pairs).astype(np.float64)
        expected = beta / np.log(np.exp(beta / atom_distances).sum(axis=1))
        np.testing.assert_allclose(dists[:, i], expected, rtol=1e-6)


@pytest.mark.parametrize("scheme", ["closest", "closest-heavy", "sidechain", "sidechain-heavy"])
@pytest.mark.parametrize("beta", [0.1, 20.0, 1000.0])
def test_contact_soft_min_close_atoms(scheme, beta):
    top = md.Topology()
    chain = top.add_chain()
    left = top.add_residue("ALA", chain)
    right = top.add_residue("ALA", chain)
    top.add_atom("CB", md.element.carbon, left)
    top.add_atom("CG", md.element.carbon, left)
    top.add_atom("CB", md.element.carbon, right)
    traj = md.Trajectory([[[0, 0, 0], [0.25, 0, 0], [0.125, 0, 0]]], top)

    with np.errstate(over="raise", divide="raise", invalid="raise"):
        distances, pairs = md.compute_contacts(
            traj,
            [[0, 1]],
            scheme=scheme,
            soft_min=True,
            soft_min_beta=beta,
        )
    # Both cross-residue distances are exactly 1/8 nm, so the sum is analytic.
    expected = beta / (beta / 0.125 + np.log(2))
    np.testing.assert_allclose(distances, [[expected]], rtol=1e-6)
    np.testing.assert_array_equal(pairs, [[0, 1]])
    assert distances.dtype == np.float32


def test_contact_soft_min_single_pair():
    top = md.Topology()
    chain = top.add_chain()
    for _ in range(2):
        residue = top.add_residue("ALA", chain)
        top.add_atom("CB", md.element.carbon, residue)
    traj = md.Trajectory([[[0, 0, 0], [8, 0, 0]]], top)
    distances, _ = md.compute_contacts(
        traj,
        [[0, 1]],
        scheme="closest",
        soft_min=True,
        soft_min_beta=0.1,
    )
    # A single atom pair must return its input distance, also for small beta.
    np.testing.assert_array_equal(distances, [[8.0]])


def test_contact_4(get_fn):
    pdb = md.load(
        get_fn("1am7_protein.pdb"),
    )  # protonated and including at least one glycine residue (which has no heavy atoms in its sidechain)

    with pytest.warns(UserWarning, match="selected topology includes at least one glycine"):
        md.compute_contacts(
            pdb,
            contacts="all",
            scheme="sidechain-heavy",
        )  # test passes if this doesn't raise an exception
