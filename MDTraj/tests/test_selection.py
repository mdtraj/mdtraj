# #############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
# Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Matthew Harrigan
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
# #############################################################################

import logging

import mdtraj
import numpy as np
from mdtraj.core.selection import SelectionParser
from mdtraj.testing import eq, get_fn


# Conda v2.0.1 build py34_0 spews a bunch of DeprecationWarnings
# from pyparsing internal code
logging.captureWarnings(True)

ala = mdtraj.load(get_fn("alanine-dipeptide-explicit.pdb"))


def make_test_topology():
    t = mdtraj.Topology()
    c = t.add_chain()

    r1 = t.add_residue("ALA", c, resSeq=5)
    r2 = t.add_residue("HOH", c, resSeq=6)

    t.add_atom("CA", mdtraj.element.carbon, r1)
    t.add_atom("H", mdtraj.element.hydrogen, r1)

    t.add_atom("O", mdtraj.element.oxygen, r2)
    t.add_atom("H1", mdtraj.element.hydrogen, r2)
    t.add_atom("H2", mdtraj.element.hydrogen, r2)

    return t


tt = make_test_topology()


def test_unary_2():
    sp = SelectionParser('all')
    for a in tt.atoms:
        assert sp.filter(a)

    sp.parse('none')
    for a in tt.atoms:
        assert not sp.filter(a)


def test_unary_3():
    sp = SelectionParser('protein or water')
    for a in tt.atoms:
        assert sp.filter(a)

    sp.parse('protein and water')
    for a in tt.atoms:
        assert not sp.filter(a)

    sp.parse('not (protein and water)')
    for a in tt.atoms:
        assert sp.filter(a)

    sp.parse('not not (protein and water)')
    for a in tt.atoms:
        assert not sp.filter(a)


def test_binary_1():
    sp = SelectionParser('resname "ALA"')
    assert sp.filter(tt.atom(0))
    assert sp.filter(tt.atom(1))

    sp.parse('mass > 2')
    assert sp.filter(tt.atom(0))
    assert not sp.filter(tt.atom(1))
    assert sp.filter(tt.atom(2))

    sp.parse('name ne O')
    assert sp.filter(tt.atom(0))
    assert not sp.filter(tt.atom(2))

def test_binary_2():
    sp = SelectionParser('name O and mass > 2')
    assert sp.filter(tt.atom(2))
    assert not sp.filter(tt.atom(3))


def test_simple():
    sp = SelectionParser("protein")
    eq(sp.mdtraj_condition, "a.residue.is_protein")
    assert sp.filter(tt.atom(0))
    assert sp.filter(tt.atom(1))
    assert not sp.filter(tt.atom(2))


def test_alias():
    sp = SelectionParser("waters")
    eq(sp.mdtraj_condition, "a.residue.is_water")
    assert sp.filter(tt.atom(3))
    assert sp.filter(tt.atom(4))
    assert not sp.filter(tt.atom(0))


def test_bool():
    sp = SelectionParser("protein or water")
    eq(sp.mdtraj_condition, "(a.residue.is_protein or a.residue.is_water)")

    sp.parse("protein or water or nucleic")
    eq(sp.mdtraj_condition,
       "(a.residue.is_protein or a.residue.is_water or a.residue.is_nucleic)")

    sp.parse("protein and backbone")
    eq(sp.mdtraj_condition, "(a.residue.is_protein and a.residue.is_backbone)")

    sp.parse("protein && backbone")
    eq(sp.mdtraj_condition, "(a.residue.is_protein and a.residue.is_backbone)")


def test_nested_bool():
    sp = SelectionParser("protein and water or nucleic")
    eq(sp.mdtraj_condition,
       "((a.residue.is_protein and a.residue.is_water) or a.residue.is_nucleic)")

    sp.parse("protein and (water or nucleic)")
    eq(sp.mdtraj_condition,
       "(a.residue.is_protein and (a.residue.is_water or a.residue.is_nucleic))")


def test_values():
    sp = SelectionParser("resid 4")
    eq(sp.mdtraj_condition, "a.residue.index == 4")

    sp.parse("resid > 4")
    eq(sp.mdtraj_condition, "a.residue.index > 4")

    sp.parse("resid gt 4")
    eq(sp.mdtraj_condition, "a.residue.index > 4")

    sp.parse("resid 5 to 8")
    eq(sp.mdtraj_condition, "5 <= a.residue.index <= 8")


def test_element():
    sp = SelectionParser()

    sp.parse("element 'O'")
    eq(sp.mdtraj_condition, "a.element.symbol == 'O'")

    sp.parse("mass 5.5 to 12.3")
    eq(sp.mdtraj_condition, "5.5 <= a.element.mass <= 12.3")


def test_not():
    sp = SelectionParser("not protein")
    eq(sp.mdtraj_condition, "(not a.residue.is_protein)")

    sp.parse("not not protein")
    eq(sp.mdtraj_condition, "(not (not a.residue.is_protein))")

    sp.parse('!protein')
    eq(sp.mdtraj_condition, "(not a.residue.is_protein)")


def test_within():
    sp = SelectionParser("within 5 of (backbone or sidechain)")
    eq(sp.mdtraj_condition,
       "(a.within == 5 of (a.residue.is_backbone or a.residue.is_sidechain))")


def test_quotes():
    should_be = "(a.name == 'CA' and a.residue.name == 'ALA')"

    sp = SelectionParser("name CA and resname ALA")
    eq(sp.mdtraj_condition, should_be)
    assert sp.filter(tt.atom(0))

    sp.parse('name "CA" and resname ALA')
    eq(sp.mdtraj_condition, should_be)
    assert sp.filter(tt.atom(0))

    sp.parse("name 'CA' and resname ALA")
    eq(sp.mdtraj_condition, should_be)
    assert sp.filter(tt.atom(0))


def test_top():
    prot = ala.topology.select("protein")
    eq(np.asarray(prot), np.arange(22))

    wat = ala.topology.select("water")
    eq(np.asarray(wat), np.arange(22, 2269))


def test_top_2():
    expr = ala.topology.select_expression("name O and water")
    eq(expr,
       "[a.index for a in top.atoms if (a.name == 'O' and a.residue.is_water)]")
