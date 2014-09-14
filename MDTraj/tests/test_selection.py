__author__ = 'harrigan'


# TODO: Expose this better
from mdtraj.core.selection import SelectionParser
from mdtraj.testing import eq

import logging

# Conda v2.0.1 build py34_0 spews a bunch of DeprecationWarnings
# from pyparsing internal code
logging.captureWarnings(True)


def test_simple():
    sp = SelectionParser("protein")
    eq(sp.unambiguous, 'residue_protein')
    eq(sp.mdtraj_condition, "a.residue.is_protein")


def test_alias():
    sp = SelectionParser("waters")
    eq(sp.unambiguous, "residue_water")
    eq(sp.mdtraj_condition, "a.residue.is_water")


def test_bool():
    sp = SelectionParser("protein or water")
    eq(sp.unambiguous, "(residue_protein or residue_water)")
    eq(sp.mdtraj_condition, "(a.residue.is_protein or a.residue.is_water)")

    sp.parse("protein or water or nucleic")
    eq(sp.unambiguous, "(residue_protein or residue_water or residue_nucleic)")
    eq(sp.mdtraj_condition,
       "(a.residue.is_protein or a.residue.is_water or a.residue.is_nucleic)")

    sp.parse("protein and backbone")
    eq(sp.unambiguous, "(residue_protein and residue_backbone)")
    eq(sp.mdtraj_condition, "(a.residue.is_protein and a.residue.is_backbone)")

    sp.parse("protein && backbone")
    eq(sp.unambiguous, "(residue_protein and residue_backbone)")
    eq(sp.mdtraj_condition, "(a.residue.is_protein and a.residue.is_backbone)")


def test_nested_bool():
    sp = SelectionParser("protein and water or nucleic")
    eq(sp.unambiguous,
       "((residue_protein and residue_water) or residue_nucleic)")
    eq(sp.mdtraj_condition,
       "((a.residue.is_protein and a.residue.is_water) or a.residue.is_nucleic)")

    sp.parse("protein and (water or nucleic)")
    eq(sp.unambiguous,
       "(residue_protein and (residue_water or residue_nucleic))")
    eq(sp.mdtraj_condition,
       "(a.residue.is_protein and (a.residue.is_water or a.residue.is_nucleic))")


def test_values():
    sp = SelectionParser("resid 4")
    eq(sp.unambiguous, "residue_index == 4")
    eq(sp.mdtraj_condition, "a.residue.index == 4")

    sp.parse("resid > 4")
    eq(sp.unambiguous, "residue_index > 4")
    eq(sp.mdtraj_condition, "a.residue.index > 4")

    sp.parse("resid gt 4")
    eq(sp.unambiguous, "residue_index > 4")
    eq(sp.mdtraj_condition, "a.residue.index > 4")

    sp.parse("resid 5 to 8")
    eq(sp.unambiguous, "residue_index == range(5 to 8)")
    eq(sp.mdtraj_condition, "5 <= a.residue.index <= 8")


def test_not():
    sp = SelectionParser("not protein")
    eq(sp.unambiguous, "(not residue_protein)")
    eq(sp.mdtraj_condition, "(not a.residue.is_protein)")

    sp.parse("not not protein")
    eq(sp.unambiguous, "(not (not residue_protein))")
    eq(sp.mdtraj_condition, "(not (not a.residue.is_protein))")

    sp.parse('!protein')
    eq(sp.unambiguous, '(not residue_protein)')
    eq(sp.mdtraj_condition, "(not a.residue.is_protein)")


def test_within():
    sp = SelectionParser("within 5 of backbone or sidechain")
    eq(sp.unambiguous,
       "(atom_within == 5 of (residue_backbone or residue_sidechain))")
