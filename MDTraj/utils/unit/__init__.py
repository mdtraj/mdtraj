"""
Physical quantities with units for dimensional analysis and automatic unit conversion.
"""
__docformat__ = "epytext en"

__author__ = "Christopher M. Bruns"
__copyright__ = "Copyright 2010, Stanford University and Christopher M. Bruns"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Christopher M. Bruns"
__email__ = "cmbruns@stanford.edu"

import ast
from .unit import Unit, is_unit
from .quantity import Quantity, is_quantity
from .unit_math import *
import unit_definitions
from .unit_definitions import *
from .constants import *
from mdtraj.utils import import_, six
UNIT_DEFINITIONS = unit_definitions
try:
    import simtk.unit as simtk_unit
except ImportError:
    pass

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
        # units into this module
        if not hasattr(unit_definitions, node.id):
            # also, let's take this opporunity to check that the node.id
            # (which supposed to be the name of the unit, like "nanometers")
            # is actually an attribute in simtk.unit
            raise ValueError('%s is not a valid unit' % node.id)

        return ast.Attribute(value=ast.Name(id='unit_definitions', ctx=ast.Load()),
                             attr=node.id, ctx=ast.Load())
_unit_context = _UnitContext()  # global instance of the visitor


def _str_to_unit(unit_string, simtk=False):
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
    # parse the string with the ast, and then run out unit context
    # visitor on it, which will basically change bare names like
    # "nanometers" into "unit.nanometers" and simulataniously check that
    # there's no nefarious stuff in the expression.

    assert isinstance(unit_string, six.string_types)
    unit_definitions = UNIT_DEFINITIONS
    if simtk:
        unit_definitions = import_('simtk.unit').unit_definitions
    parsed = ast.parse(unit_string, mode='eval')
    node = _unit_context.visit(parsed)
    fixed_node = ast.fix_missing_locations(node)
    output = eval(compile(fixed_node, '<string>', mode='eval'), {}, locals())
    return output


def in_units_of(quantity, units_in, units_out, inplace=False):
    """Convert a quantity between unit systems

    Parameters
    ----------
    quantity : number, np.ndarray, or simtk.unit.Quantity
        quantity can either be a unitted quantity -- i.e. instance of
        simtk.unit.Quantity, or just a bare number or numpy array
    units_in : str
        If you supply a quantity that's not a simtk.unit.Quantity, you should
        tell me what units it is in. If you don't, i'm just going to echo you
        back your quantity without doing any unit checking.
    units_out : str
        A string description of the units you want out. This should look
        like "nanometers/picosecondsecond" or "nanometers**3" or whatever
    inplace : bool
        Do the transformation inplace. This will only work if the quantity
        is a mutable type, like a numpy array.

    Examples
    --------
    >>> in_units_of(1*units.meter**2/units.second, 'nanometers**2/picosecond')
    1000000.0
    """
    if quantity is None:
        return quantity

    if 'simtk.unit' in sys.modules and isinstance(quantity, simtk_unit.Quantity):
        units_in = quantity.unit
        units_out = _str_to_unit(units_out, simtk=True)
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

    print 'in, out', units_in, units_out

    factor = units_in.conversion_factor_to(units_out)
    print 'factor', factor, 'donefactor'

    if not inplace:
        return quantity * factor

    quantity *= factor

