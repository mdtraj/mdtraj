# #############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
# Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Matthew Harrigan
# Contributors: Carlos Xavier Hernandez
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

import ast
import sys

import numpy as np
import pytest
from pytest import MonkeyPatch

import mdtraj
from mdtraj.core.selection import parse_selection
from mdtraj.testing import eq


use_astunparse = False


def pnode(s):
    return ast.parse(s, mode="eval").body


@pytest.fixture()
def ala(get_fn):
    return mdtraj.load(get_fn("alanine-dipeptide-explicit.pdb"))


@pytest.fixture()
def gbp(get_fn):
    return mdtraj.load(get_fn("2EQQ.pdb"))


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


def test_range():
    eq(
        parse_selection("resSeq 1 to 10").astnode,
        pnode("1 <= atom.residue.resSeq <= 10"),
    )

    sp = parse_selection("resSeq 5 to 6")
    for a in tt.atoms:
        assert sp.expr(a)
    sp = parse_selection("resSeq 7 to 8")
    for a in tt.atoms:
        assert not sp.expr(a)


def test_unary_2():
    sp = parse_selection("all")
    for a in tt.atoms:
        assert sp.expr(a)

    sp = parse_selection("none")
    for a in tt.atoms:
        assert not sp.expr(a)


def test_unary_3():
    sp = parse_selection("protein or water")

    for a in tt.atoms:
        assert sp.expr(a)

    sp = parse_selection("protein and water")
    for a in tt.atoms:
        assert not sp.expr(a)

    sp = parse_selection("not (protein and water)")
    for a in tt.atoms:
        assert sp.expr(a)

    sp = parse_selection("not not (protein and water)")
    for a in tt.atoms:
        assert not sp.expr(a)

def test_binary_1():
    sp = parse_selection('resname "ALA"')
    assert sp.expr(tt.atom(0))
    assert sp.expr(tt.atom(1))

    sp = parse_selection("rescode A")
    assert sp.expr(tt.atom(0))
    assert sp.expr(tt.atom(1))

    sp = parse_selection("mass > 2")
    assert sp.expr(tt.atom(0))
    assert not sp.expr(tt.atom(1))
    assert sp.expr(tt.atom(2))

    sp = parse_selection("name ne O")
    assert sp.expr(tt.atom(0))
    assert not sp.expr(tt.atom(2))


def test_binary_2():
    sp = parse_selection("name O and mass > 2")
    assert sp.expr(tt.atom(2))
    assert not sp.expr(tt.atom(3))


def test_unary_1():
    eq(parse_selection("all").astnode, pnode("True"))
    eq(parse_selection("everything").astnode, pnode("True"))
    eq(parse_selection("none").astnode, pnode("False"))
    eq(parse_selection("nothing").astnode, pnode("False"))
    # eq(parse_selection('nucleic').astnode, pnode('atom.residue.is_nucleic'))
    # eq(parse_selection('is_nucleic').astnode, pnode('atom.residue.is_nucleic'))
    eq(parse_selection("protein").astnode, pnode("atom.residue.is_protein"))
    eq(parse_selection("is_protein").astnode, pnode("atom.residue.is_protein"))
    eq(parse_selection("water").astnode, pnode("atom.residue.is_water"))
    eq(parse_selection("is_water").astnode, pnode("atom.residue.is_water"))
    eq(parse_selection("waters").astnode, pnode("atom.residue.is_water"))


def test_binary_selection_operator():
    eq(parse_selection("name < 1").astnode, pnode("atom.name < 1"))
    eq(parse_selection("name lt 1").astnode, pnode("atom.name < 1"))
    eq(parse_selection("name > 1").astnode, pnode("atom.name > 1"))
    eq(parse_selection("name gt 1").astnode, pnode("atom.name > 1"))
    eq(parse_selection("name == 1").astnode, pnode("atom.name == 1"))
    eq(parse_selection("name eq 1").astnode, pnode("atom.name == 1"))
    eq(parse_selection("name != 1").astnode, pnode("atom.name != 1"))
    eq(parse_selection("name ne 1").astnode, pnode("atom.name != 1"))
    eq(parse_selection("name >= 1").astnode, pnode("atom.name >= 1"))
    eq(parse_selection("name ge 1").astnode, pnode("atom.name >= 1"))
    eq(parse_selection("name <= 1").astnode, pnode("atom.name <= 1"))
    eq(parse_selection("name le 1").astnode, pnode("atom.name <= 1"))

    eq(parse_selection("1 == name").astnode, pnode("1 == atom.name"))
    eq(parse_selection("1 eq name").astnode, pnode("1 == atom.name"))


def test_raises():
    pytest.raises(ValueError, lambda: parse_selection("or"))
    pytest.raises(ValueError, lambda: parse_selection("a <"))


def test_raises2():
    pytest.raises(ValueError, lambda: parse_selection("dog 5"))
    pytest.raises(ValueError, lambda: parse_selection("dog == 5"))
    pytest.raises(ValueError, lambda: parse_selection("dog frog"))
    pytest.raises(ValueError, lambda: parse_selection("not dog"))
    pytest.raises(ValueError, lambda: parse_selection("protein or dog"))
    pytest.raises(ValueError, lambda: parse_selection("dog 1 to 5"))
    pytest.raises(ValueError, lambda: parse_selection("dog"))


def test_top(ala):
    prot = ala.topology.select("protein")
    eq(np.asarray(prot), np.arange(22))

    wat = ala.topology.select("water")
    eq(np.asarray(wat), np.arange(22, 2269))


def test_backbone(gbp):
    ref_backbone = gbp.topology.select(
        "protein and (name C or name CA or " "name N or name O)",
    )
    backbone = gbp.topology.select("backbone")
    is_backbone = gbp.topology.select("is_backbone")

    eq(np.asarray(backbone), np.asarray(ref_backbone))
    eq(np.asarray(is_backbone), np.asarray(ref_backbone))


def test_sidechain(gbp):
    ref_sidechain = gbp.topology.select(
        "protein and not (name C or name CA " "or name N or name O or name H or name HA)",
    )
    sidechain = gbp.topology.select("sidechain")
    is_sidechain = gbp.topology.select("is_sidechain")

    eq(np.asarray(sidechain), np.asarray(ref_sidechain))
    eq(np.asarray(is_sidechain), np.asarray(ref_sidechain))


def test_literal(gbp):
    name_og1_0 = gbp.topology.select('name "OG1"')
    name_og1_1 = gbp.topology.select("name 'OG1'")
    name_og1_2 = gbp.topology.select("name OG1")

    ref_og1 = np.asarray([a.index for a in gbp.topology.atoms if a.name == "OG1"])
    eq(name_og1_0, ref_og1)
    eq(name_og1_1, ref_og1)
    eq(name_og1_2, ref_og1)


class TestSelectionUnparse():
    def test_alias(self):
        sp = parse_selection("waters")
        if use_astunparse:
            eq(sp.source, "atom.residue.is_water\n")
        else:
            eq(sp.source, "atom.residue.is_water")

        assert sp.expr(tt.atom(3))
        assert sp.expr(tt.atom(4))
        assert not sp.expr(tt.atom(0))

    def test_simple(self):
        sp = parse_selection("protein")
        if use_astunparse:
            eq(sp.source, "atom.residue.is_protein\n")
        else:
            eq(sp.source, "atom.residue.is_protein")

        assert sp.expr(tt.atom(0))
        assert sp.expr(tt.atom(1))
        assert not sp.expr(tt.atom(2))

    def test_bool(self):
        sp = parse_selection("protein or water")
        sp2 = parse_selection("protein or water or all\n")

        if use_astunparse:
            eq(sp.source, "(atom.residue.is_protein or atom.residue.is_water)\n")
            eq(
                sp2.source,
                "(atom.residue.is_protein or atom.residue.is_water or True)\n",
            )
        else:
            eq(sp.source, "atom.residue.is_protein or atom.residue.is_water")
            eq(
                sp2.source,
                "atom.residue.is_protein or atom.residue.is_water or True",
            )

    def test_nested_bool(self):
        sp = parse_selection("nothing and water or all")
        sp2 = parse_selection("nothing and (water or all)")

        if use_astunparse:
            eq(
                sp.source,
                "((False and atom.residue.is_water) or True)\n",
            )

            eq(
                sp2.source,
                "(False and (atom.residue.is_water or True))\n",
            )
        else:
            eq(
                sp.source,
                "False and atom.residue.is_water or True",
            )

            eq(
                sp2.source,
                "False and (atom.residue.is_water or True)",
            )

    def test_values(self):
        sp = parse_selection("resid 4")
        sp2 = parse_selection("chainid 4")
        sp3 = parse_selection("resid > 4")
        sp4 = parse_selection("resid gt 4")
        sp5 = parse_selection("resid 5 to 8")

        if use_astunparse:
            eq(sp.source, "(atom.residue.index == 4)\n")
            eq(sp2.source, "(atom.residue.chain.index == 4)\n")
            eq(sp3.source, "(atom.residue.index > 4)\n")
            eq(sp4.source, "(atom.residue.index > 4)\n")
            eq(sp5.source, "(5 <= atom.residue.index <= 8)\n")
        else:
            eq(sp.source, "atom.residue.index == 4")
            eq(sp2.source, "atom.residue.chain.index == 4")
            eq(sp3.source, "atom.residue.index > 4")
            eq(sp4.source, "atom.residue.index > 4")
            eq(sp5.source, "5 <= atom.residue.index <= 8")

    def test_element(self):
        sp = parse_selection("element 'O'")
        sp2 = parse_selection("mass 5.5 to 12.3")

        if use_astunparse:
            eq(sp.source, "(atom.element.symbol == 'O')\n")
            eq(sp2.astnode, pnode("(5.5 <= atom.element.mass <= 12.3)"))
        else:
            eq(sp.source, "atom.element.symbol == 'O'")
            eq(sp2.astnode, pnode("5.5 <= atom.element.mass <= 12.3"))

    def test_not(self):
        sp = parse_selection("not protein")
        sp2 = parse_selection("not not protein")
        sp3 = parse_selection("!protein")

        if use_astunparse:
            eq(sp.source, "(not atom.residue.is_protein)\n")
            eq(sp2.source, "(not (not atom.residue.is_protein))\n")
            eq(sp3.source, "(not atom.residue.is_protein)\n")
        else:
            eq(sp.source, "not atom.residue.is_protein")
            eq(sp2.source, "not not atom.residue.is_protein")
            eq(sp3.source, "not atom.residue.is_protein")

    def test_re(self):
        sp = parse_selection("name =~ 'C.*'")
        sp2 = parse_selection("(name =~ 'C.*') and all")

        if use_astunparse:
            eq(sp.source, "(re.match('C.*', atom.name) is not None)\n")
            eq(sp2.source, "((re.match('C.*', atom.name) is not None) and True)\n")
        else:
            eq(sp.source, "re.match('C.*', atom.name) is not None")
            eq(sp2.source, "re.match('C.*', atom.name) is not None and True")

    def test_quotes(self):
        if use_astunparse:
            should_be = "((atom.name == 'CA') and (atom.residue.name == 'ALA'))\n"
        else:
            should_be = "atom.name == 'CA' and atom.residue.name == 'ALA'"

        sp = parse_selection("name CA and resname ALA")
        eq(sp.source, should_be)
        assert sp.expr(tt.atom(0))

        sp = parse_selection('name "CA" and resname ALA')
        eq(sp.source, should_be)
        assert sp.expr(tt.atom(0))

        sp = parse_selection("name 'CA' and resname ALA")
        eq(sp.source, should_be)
        assert sp.expr(tt.atom(0))

    def test_top_2(self, ala):
        expr = ala.topology.select_expression("name O and water")

        if use_astunparse:
            eq(
                expr,
                "[atom.index for atom in topology.atoms if ((atom.name == 'O') and atom.residue.is_water)\n]",
            )
        else:
            eq(
                expr,
                "[atom.index for atom in topology.atoms if atom.name == 'O' and atom.residue.is_water]",
            )

    def test_in(self):
        sp = parse_selection("resname ALA ASP GLU")
        sp2 = parse_selection("resid 100 101 102")

        if use_astunparse:
            eq(sp.source, "(atom.residue.name in ['ALA', 'ASP', 'GLU'])\n")
            eq(sp2.source, "(atom.residue.index in [100, 101, 102])\n")
        else:
            eq(sp.source, "atom.residue.name in ['ALA', 'ASP', 'GLU']")
            eq(sp2.source, "atom.residue.index in [100, 101, 102]")

    @pytest.mark.skip('failed to parse')
    def test_within(self):
        sp = parse_selection("within 5 of (backbone or sidechain)")
        if use_astunparse:
            eq(sp.source,
                "(atom.within == 5 of (atom.residue.is_backbone or atom.residue.is_sidechain))")
        else:
            eq(sp.source,
                "atom.within == 5 of atom.residue.is_backbone or atom.residue.is_sidechain")

class TestSelectionASTunparse(TestSelectionUnparse):
    def setup_method(self, method):
        monkeypatch = MonkeyPatch()
        import astunparse
        monkeypatch.setitem(globals(), 'use_astunparse', True)
        monkeypatch.setattr('mdtraj.core.selection.unparse', astunparse.unparse)

    def teardown_method(self, method):
        monkeypatch = MonkeyPatch()
        monkeypatch.undo()
        use_astunparse = False
