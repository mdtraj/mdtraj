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
"""

Unit processing for MDTraj. This subpackage is a port of openmm.unit from
OpenMM. Unlike in openmm, the MDTraj library **does not pass around
"united quantities"**

The only publicly facing API from this package, for the purpose of MDTraj,
is "in_units_of", which does unit conversion of numbers or numpy arrays
where the input and output units are passed as strings.

"""

import ast
import sys

import numpy as np

from mdtraj.utils import import_
from mdtraj.utils.unit import unit_definitions
from mdtraj.utils.unit.quantity import Quantity

UNIT_DEFINITIONS = unit_definitions

try:
    import openmm.unit as openmm_unit
except ImportError:
    pass

__all__ = ["in_units_of"]


class _UnitContext(ast.NodeTransformer):
    """Node transformer for an AST hack that turns raw strings into
    complex simt.unit.Unit expressions. See _str_to_unit for how this
    is used -- it's not really meant to stand on its own
    """

    # we want to do some validation to ensure that the AST only
    # contains "safe" operations. These are operations that can reasonably
    # appear in unit expressions
    allowed_ops = [
        ast.Expression,
        ast.BinOp,
        ast.Name,
        ast.Attribute,
        ast.Pow,
        ast.Div,
        ast.Mult,
        ast.Constant,
    ]

    def visit(self, node):
        if not any(isinstance(node, a) for a in self.allowed_ops):
            raise ValueError(
                "Invalid unit expression. Contains dissallowed " "operation %s" % node.__class__.__name__,
            )
        return super().visit(node)

    def visit_Name(self, node):
        # we want to prefix all names to look like unit.nanometers instead
        # of just "nanometers", because I don't want to import * from
        # units into this module
        if not hasattr(unit_definitions, node.id):
            # also, let's take this opporunity to check that the node.id
            # (which supposed to be the name of the unit, like "nanometers")
            # is actually an attribute in openmm.unit
            raise ValueError("%s is not a valid unit" % node.id)

        return ast.Attribute(
            value=ast.Name(id="unit_definitions", ctx=ast.Load()),
            attr=node.id,
            ctx=ast.Load(),
        )


_unit_context = _UnitContext()  # global instance of the visitor


def _str_to_unit(unit_string, openmm=False):
    """eval() based transformer that extracts a openmm.unit object
    from a string description.

    Parameters
    ----------
    unit_string : str
        string description of a unit. this may contain expressions with
        multiplication, division, powers, etc.

    Examples
    --------
    >>> type(_str_to_unit('nanometers**2/meters*gigajoules'))
    <class 'openmm.unit.unit.Unit'>
    >>> str(_str_to_unit('nanometers**2/meters*gigajoules'))
    'nanometer**2*gigajoule/meter'

    """
    # parse the string with the ast, and then run out unit context
    # visitor on it, which will basically change bare names like
    # "nanometers" into "unit.nanometers" and simulataniously check that
    # there's no nefarious stuff in the expression.

    assert isinstance(unit_string, str)
    unit_definitions = UNIT_DEFINITIONS
    if openmm:
        unit_definitions = openmm_unit.unit_definitions
    parsed = ast.parse(unit_string, mode="eval")
    node = _unit_context.visit(parsed)
    fixed_node = ast.fix_missing_locations(node)
    output = eval(compile(fixed_node, "<string>", mode="eval"), {}, locals())
    return output


def in_units_of(quantity, units_in, units_out, inplace=False):
    """Convert a numerical quantity between unit systems.

    Parameters
    ----------
    quantity : {number, np.ndarray, openmm.unit.Quantity}
        quantity can either be a unitted quantity -- i.e. instance of
        openmm.unit.Quantity, or just a bare number or numpy array
    units_in : str
        If you supply a quantity that's not a openmm.unit.Quantity, you should
        tell me what units it is in. If you don't, i'm just going to echo you
        back your quantity without doing any unit checking.
    units_out : str
        A string description of the units you want out. This should look
        like "nanometers/picosecond" or "nanometers**3" or whatever
    inplace : bool
        Attempt to do the transformation inplace, by mutating the `quantity`
        argument and avoiding a copy. This is only possible if `quantity` is a
        writable numpy array.

    Returns
    -------
    rquantity : {number, np.ndarray}
        The resulting quantity, in the new unit system. If the function was
        called with `inplace=True` and `quantity` was a writable numpy array,
        `rquantity` will alias the same memory as the input `quantity`, which
        will have been changed inplace. Otherwise, if a copy was required,
        `rquantity` will point to new memory.

    Examples
    --------
    >>> in_units_of(1, 'meter**2/second', 'nanometers**2/picosecond')
    1000000.0
    """
    if quantity is None:
        return quantity

    if "openmm.unit" in sys.modules and isinstance(quantity, openmm_unit.Quantity):
        units_in = quantity.unit
        units_out = _str_to_unit(units_out, openmm=True)
        quantity = quantity._value
    elif isinstance(quantity, Quantity):
        units_in = quantity.unit
        units_out = _str_to_unit(units_out)
        quantity = quantity._value
    else:
        if units_in is None:
            return quantity
        units_in = _str_to_unit(units_in)
        units_out = _str_to_unit(units_out)

    if not units_in.is_compatible(units_out):
        raise TypeError(f'Unit "{units_in}" is not compatible with Unit "{units_out}".')

    factor = units_in.conversion_factor_to(units_out)
    if inplace and (isinstance(quantity, np.ndarray) and quantity.flags["WRITEABLE"]):
        quantity *= factor
        return quantity
    return quantity * factor
