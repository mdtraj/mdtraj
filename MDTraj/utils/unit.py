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

"""Automatic unit conversion using simtk.units behind the scenes.
"""
##############################################################################
# imports
#############################################################################

from __future__ import print_function, division
import ast
from mdtraj.utils import import_

# Debugging
# from meta.asttools import print_ast

__all__ = ['in_units_of', 'convert']

##############################################################################
# classes and functions
##############################################################################


def convert(quantity, in_unit, out_unit, inplace=False):
    """Convert between distance units. Note, inplace==True won't
    work unless the quantity is a numpy array. This is a simple
    version that doesn't require any external imports.
    """
    if quantity is None:
        return None

    factor = {('angstroms', 'angstroms'): 1,
              ('nanometers', 'nanometers'): 1,
              ('angstroms', 'nanometers'): 0.1,
              ('nanometers', 'angstroms'): 10}[(in_unit, out_unit)]
    if not inplace:
        return quantity * factor
    quantity *= factor



class _UnitContext(ast.NodeTransformer):
    """Node transformer for an AST hack that turns raw strings into
    complex simt.unit.Unit expressions. See _str_to_unit for how this
    is used -- it's not really meant to stand on its own
    """
    # we want to do some validation to ensure that the AST only
    # contains "safe" operations. These are operations that can reasonably
    # appear in unit expressions
    allowed_ops = [ast.Expression, ast.BinOp, ast.Name, ast.Attribute,
                   ast.Pow, ast.Div, ast.Mult, ast.Num]

    def visit(self, node):
        if not any(isinstance(node, a) for a in self.allowed_ops):
            raise ValueError('Invalid unit expression. Contains dissallowed '
                             'operation %s' % node.__class__.__name__)
        return super(_UnitContext, self).visit(node)

    def visit_Name(self, node):
        # we want to prefix all names to look like unit.nanometers instead
        # of just "nanometers", because I don't want to import * from
        # units into this module.
        units = import_('simtk.unit')
        if not (node.id == 'units' or hasattr(units, node.id)):
            # also, let's take this opporunity to check that the node.id
            # (which supposed to be the name of the unit, like "nanometers")
            # is actually an attribute in simtk.unit
            raise ValueError('%s is not a valid unit' % node.id)

        return ast.Attribute(value=ast.Name(id='units', ctx=ast.Load()),
                             attr=node.id, ctx=ast.Load())
_unit_context = _UnitContext()  # global instance of the visitor


def _str_to_unit(unit_string):
    """eval() based transformer that extracts a simtk.unit object
    from a string description.

    Parameters
    ----------
    unit_string : str
        string description of a unit. this may contain expressions with
        multiplication, division, powers, etc.

    Examples
    --------
    >>> type(_str_to_unit('nanometers**2/meters*gigajoules'))
    <class 'simtk.unit.unit.Unit'>
    >>> str(_str_to_unit('nanometers**2/meters*gigajoules'))
    'nanometer**2*gigajoule/meter'

    """
    units = import_('simtk.unit')
    # parse the string with the ast, and then run out unit context
    # visitor on it, which will basically change bare names like
    # "nanometers" into "unit.nanometers" and simulataniously check that
    # there's no nefarious stuff in the expression.


    node = _unit_context.visit(ast.parse(unit_string, mode='eval'))
    fixed_node = ast.fix_missing_locations(node)
    output = eval(compile(fixed_node, '<string>', mode='eval'))

    return output


def in_units_of(quantity, units_out, units_in=None):
    """Convert a quantity between unit systems

    Parameters
    ----------
    quantity : number, np.ndarray, or simtk.unit.Quantity
        quantity can either be a unitted quantity -- i.e. instance of
        simtk.unit.Quantity, or just a bare number or numpy array
    units_out : str
        A string description of the units you want out. This should look
        like "nanometers/picosecondsecond" or "nanometers**3" or whatever
    units_in : str
        If you supply a quantity that's not a simtk.unit.Quantity, you should
        tell me what units it is in. If you don't, i'm just going to echo you
        back your quantity without doing any unit checking.

    Examples
    --------
    >>> in_units_of(1*units.meter**2/units.second, 'nanometers**2/picosecond')
    1000000.0
    """
    units = import_('simtk.unit')

    if quantity is None:
        return quantity

    if isinstance(quantity, units.Quantity):
        return quantity.value_in_unit(_str_to_unit(units_out))
    else:
        if units_in is None:
            return quantity
        united_quantity = units.Quantity(quantity, _str_to_unit(units_in))
        return united_quantity.value_in_unit(_str_to_unit(units_out))
