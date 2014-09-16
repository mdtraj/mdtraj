# #############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
# Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Matthew Harrigan
# Contributors: Robert T. McGibbon
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
from mdtraj.utils import import_
from collections import namedtuple

NUMS = '.0123456789'
ATOM_NAME = '__ATOM__'

class _RewriteNames(ast.NodeTransformer):
    def visit_Name(self, node):
        if node.id == ATOM_NAME:
            # when we're building the AST, we refer to the current atom using
            # ATOM_NAME, but then before returning the ast to the user, we
            # rewrite the name as just 'atom'.
            return ast.Name(id='atom', ctx=ast.Load())
        # all other bare names are taken to be string literals. Thus something
        # like parse_selection('name CA') properly resolves CA as a string
        # literal, not a barename to be loaded from the global scope!
        return ast.Str(s=node.id)

_ParsedSelection = namedtuple('_ParsedSelection', ['expr', 'source', 'astnode'])

__all__ = ['parse_selection']

def _kw(*tuples):
    """Create a many-to-one dictionary.

    _kw((['one', '1'], 'one'))
    gives {'one': 'one', '1': 'one'}
    """
    dic = dict()
    for keys, val in tuples:
        for key in keys:
            dic[key] = val
    return dic


def _check_n_tokens(tokens, n_tokens, name):
    if not len(tokens) == n_tokens:
        err = "{} take 3 values. You gave {}"
        err = err.format(name, len(tokens))
        ParseException = import_('pyparsing').ParseException
        raise ParseException(err)


class UnarySelectionOperand(object):
    """Unary selections

    Examples
    --------
    'all' -> lambda atom: True
    'none' -> lambda atom: False
    'protein' -> lambda atom: atom.residue.is_protein()
    'water' -> lambda atom: atom.residue.is_water()
    """

    keyword_aliases = _kw(
        # Atom.<attribute>
        (['all', 'everything'], ast.Name(id='True', ctx=ast.Load())),
        (['none', 'nothing'], ast.Name(id='False', ctx=ast.Load())),
        # Atom.residue.<attribute>
        (['protein', 'is_protein'], ('residue', 'is_protein')),
        (['nucleic', 'is_nucleic'], ('residue', 'is_nucleic')),
        (['water', 'waters', 'is_water'], ('residue', 'is_water')),
    )

    def __init__(self, tokens):
        self._tokens = tokens
        _check_n_tokens(tokens, 1, 'Unary selectors')


    def ast(self):
        rhs = self.keyword_aliases[self._tokens[0]]
        if isinstance(rhs, ast.AST):
            return rhs

        # we structure the rhs to be an iterable, which we recursively
        # apply ast.Attribute() on. This transforms ('residue', 'is_protein'),
        # for example, into
        # Attribute(value=Attribute(value=Name(id=ATOM_NAME, ctx=Load()),
        #  attr='residue', ctx=Load()), attr='is_protein', ctx=Load())
        left = ast.Name(id=ATOM_NAME, ctx=ast.Load())
        for attr in rhs:
            left = ast.Attribute(value=left, attr=attr, ctx=ast.Load())

        return left


class BinarySelectionOperand(object):
    """Binary selections

    Examples
    --------
    'name CA' -> lambda atom: atom.name == 'CA'
    'index > 1' -> lambda atom: atom.index > 1
    'resname != TYR' -> lambda atom: atom.residue.name != 'TYR'
    """

    operator_aliases = _kw(
        (['<', 'lt'], ast.Lt),
        (['==', 'eq'], ast.Eq),
        (['<=', 'le'], ast.LtE),
        (['!=', 'ne'], ast.NotEq),
        (['>=', 'ge'], ast.GtE),
        (['>', 'gt'], ast.Gt),
    )

    keyword_aliases = _kw(
        (('name',),             ('name',)),
        (('index',),            ('index',)),
        (('numbonds',),         ('numbonds',)),
        (('type', 'element'),   ('element', 'symbol')),
        (('radius',),           ('element', 'radius')),
        (('mass',),             ('element', 'mass')),
        (('residue', 'resSeq'), ('residue', 'resSeq')),
        (('resname',),          ('residue', 'name')),
        (('resid',),            ('residue', 'index')),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        self._tokens = tokens
        _check_n_tokens(tokens, 3, 'Binary selectors')

    def ast(self):
        left = ast.Name(id=ATOM_NAME, ctx=ast.Load())
        for attr in self.keyword_aliases[self._tokens[0]]:
            left = ast.Attribute(value=left, attr=attr, ctx=ast.Load())

        ops = [self.operator_aliases[self._tokens[1]]()]
        comparators = [ast.parse(self._tokens[2], mode='eval').body]

        return ast.Compare(left=left, ops=ops, comparators=comparators)


class UnaryInfixOperand(object):
    n_terms = 1
    assoc = 'RIGHT'

    keyword_aliases = _kw(
        (['not', '!'], ast.Not()),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        _check_n_tokens(tokens, 2, 'Unary infix operators')
        self._op, self._value = tokens

    def ast(self):
        return ast.UnaryOp(op=self.keyword_aliases[self._op],
                           operand=self._value.ast())


class BinaryInfixOperand(object):
    n_terms = 2
    assoc = 'LEFT'

    keyword_aliases = _kw(
        (['and', '&&'], ast.And()),
        (['or', '||'], ast.Or()),
    )

    def __init__(self, tokens):
        tokes = tokens[0]
        if len(tokens) % 2 == 1:
            self._op = tokes[1]
            self._parts = tokes[::2]
        else:
            err = "Invalid number of infix expressions: {}"
            err = err.format(len(tokens))
            ParseException = import_('pyparsing').ParseException
            raise ParseException(err)

    def ast(self):
        return ast.BoolOp(op=self.keyword_aliases[self._op],
                          values=[e.ast() for e in self._parts])


class parse_selection(object):
    """Parse an atom selection expression

    Parameters
    ----------
    selection_string : str
        Selection string, a string in the MDTraj atom selection grammer.

    Returns
    -------
    expr : callable (atom -> bool)
        A callable object which accepts an MDTraj.core.topology.Atom object and
        returns a boolean value giving whether or not that particular atom
        satisfies the selection string.
    source : str
        Python source code corresponding to the expression ``expr``.
    astnode : ast.AST
        Python abstract syntax tree node containing the parsed expression

    Examples
    --------
    >>> expr, source, astnode = parse_selection('protein and type CA')
    >>> expr
    <function __main__.<lambda>>
    >>> source
    '(atom.residue.is_protein and (atom.element.symbol == CA))'
    >>> <_ast.BoolOp at 0x103969d50>
    """
    def __init__(self):
        self.is_initialized = False
        self.expression = None

    def _initialize(self):
        pp = import_('pyparsing')

        def keywords(klass):
            return pp.MatchFirst([pp.Keyword(kw) for kw in klass.keyword_aliases.keys()])

        def infix(klass):
            return [(kw, klass.n_terms, getattr(pp.opAssoc, klass.assoc), klass)
                for kw in klass.keyword_aliases.keys()]

        comparison_op = pp.oneOf(list(BinarySelectionOperand.operator_aliases.keys()))
        comparison_op = pp.Optional(comparison_op, '==')
        value = pp.Word(NUMS) | pp.quotedString | pp.Word(pp.alphas)

        unary = keywords(UnarySelectionOperand)
        unary.setParseAction(UnarySelectionOperand)

        binary = pp.Group(
            keywords(BinarySelectionOperand) + comparison_op + value)
        binary.setParseAction(BinarySelectionOperand)
        expression = pp.MatchFirst([unary, binary])

        logical_expr = pp.infixNotation(expression,
            infix(BinaryInfixOperand) + infix(UnaryInfixOperand))

        self.expression = logical_expr
        self.is_initialized = True

        self.transformer = _RewriteNames()

    def __call__(self, selection):
        if not self.is_initialized:
            self._initialize()

        parse_result = self.expression.parseString(selection)

        astnode = self.transformer.visit(parse_result[0].ast())


        func = ast.Expression(body=ast.Lambda(
            args=ast.arguments(args=[ast.Name(id=ATOM_NAME, ctx=ast.Param())],
                              vararg=None, kwarg=None, defaults=[]),
            body=astnode))
        expr = eval(
            compile(ast.fix_missing_locations(func), '<string>', mode='eval'))

        try:
            codegen = import_('codegen')
            source = codegen.to_source(astnode)
        except ImportError:
            source = None

        return _ParsedSelection(expr, source, astnode)

# Create the callable, and use it to overshadow the class. this way there's
# basically just one global instance of the "function", even thought its
# a callable class.
parse_selection = parse_selection()
