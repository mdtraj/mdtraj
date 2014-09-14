__author__ = 'harrigan'


# TODO: Expose this better
from mdtraj.core.selection import SelectionParser
from mdtraj.testing import eq

import logging

# Conda v2.0.1 build py34_0 spews a bunch of DeprecationWarnings
# from pyparsing internal code
logging.captureWarnings(True)


def test_select_simple():
    sp = SelectionParser("protein")
    eq(sp.unambiguous, 'residue_protein')
    eq(sp.mdtraj_expression, "a.residue.is_protein")


def test_select_alias():
    sp = SelectionParser("waters")
    eq(sp.unambiguous, "residue_water")
    eq(sp.mdtraj_expression, "a.residue.is_water")


def test_select_bool():
    sp = SelectionParser("protein or water")
    eq(sp.unambiguous, "(residue_protein or residue_water)")
    # TODO: eq(sp.mdtraj_expression, "(a.residue.is_protein or a.residue.is_water")

    sp.parse("protein or water or nucleic")
    eq(sp.unambiguous, "(residue_protein or residue_water or residue_nucleic)")
    # TODO: eq(sp.mdtraj_expression, "(a.residue.is_protein or a.residue.is_water or a.residue.is_nucleic)")

    sp.parse("protein and backbone")
    eq(sp.unambiguous, "(residue_protein and residue_backbone)")
    # TODO: eq(sp.mdtraj_expression, "(a.residue.is_protein and a.residue.is_backbone)")


def test_select_nested_bool():
    sp = SelectionParser("protein and water or nucleic")
    eq(sp.unambiguous,
       "((residue_protein and residue_water) or residue_nucleic)")

    sp.parse("protein and (water or nucleic)")
    eq(sp.unambiguous,
       "(residue_protein and (residue_water or residue_nucleic))")


def test_select_values():
    sp = SelectionParser("resid 4")
    eq(sp.unambiguous, "residue_index == 4")

    sp.parse("resid > 4")
    eq(sp.unambiguous, "residue_index > 4")

    sp.parse("resid gt 4")
    eq(sp.unambiguous, "residue_index > 4")

    sp.parse("resid 5 to 8")
    eq(sp.unambiguous, "residue_index == range(5 to 8)")


def test_select_not():
    sp = SelectionParser("not protein")
    eq(sp.unambiguous, "(not residue_protein)")

    sp.parse("not not protein")
    eq(sp.unambiguous, "(not (not residue_protein))")

