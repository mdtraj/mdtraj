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
from copy import deepcopy
from mdtraj.utils import import_
from collections import namedtuple

__all__ = ['parse_selection']

#############################################################################
# Globals
#############################################################################

NUMS = '.0123456789'
ATOM_NAME = '__ATOM__'
_ParsedSelection = namedtuple('_ParsedSelection', ['expr', 'source', 'astnode'])

#############################################################################
# Utils
#############################################################################

class _RewriteNames(ast.NodeTransformer):
    def visit_Name(self, node):
        _safe_names = {'None': None, 'True': True, 'False': False}
        if node.id in _safe_names:
            return node

        if node.id == ATOM_NAME:
            # when we're building the AST, we refer to the current atom using
            # ATOM_NAME, but then before returning the ast to the user, we
            # rewrite the name as just 'atom'.
            return ast.Name(id='atom', ctx=ast.Load())
        # all other bare names are taken to be string literals. Thus something
        # like parse_selection('name CA') properly resolves CA as a string
        # literal, not a barename to be loaded from the global scope!
        return ast.Str(s=node.id)


def _chained_atom_attr(*attrs):
    # This transforms, for example, ('residue', 'is_protein'),
    # into
    # Attribute(value=Attribute(value=Name(id=ATOM_NAME, ctx=Load()),
    #  attr='residue', ctx=Load()), attr='is_protein', ctx=Load())
    left = ast.Name(id=ATOM_NAME, ctx=ast.Load())
    for attr in attrs:
        left = ast.Attribute(value=left, attr=attr, ctx=ast.Load())
    return left


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

    This class implements one simplest component of the selection grammar,
    unary selections, enabling expressions like 'all', or 'protein'.

    Examples
    --------
    'all' -> lambda atom: True
    'none' -> lambda atom: False
    'protein' -> lambda atom: atom.residue.is_protein()
    'water' -> lambda atom: atom.residue.is_water()
    """

    keyword_aliases = _kw(
        # Atom.<attribute>
        (('all', 'everything'), ast.Name(id='True', ctx=ast.Load())),
        (('none', 'nothing'), ast.Name(id='False', ctx=ast.Load())),
        # Atom.residue.<attribute>
        (('protein', 'is_protein'), _chained_atom_attr('residue', 'is_protein')),
        (('nucleic', 'is_nucleic'), _chained_atom_attr('residue', 'is_nucleic')),
        (('water', 'waters', 'is_water'), _chained_atom_attr('residue', 'is_water')),
    )

    def __init__(self, tokens):
        # pyparsing constructs the instance while building the parse tree,
        # and gives us the set tokens. In this case, the tokens are
        self._tokens = tokens
        _check_n_tokens(tokens, 1, 'Unary selectors')
        assert tokens[0] in self.keyword_aliases

    def ast(self):
        return deepcopy(self.keyword_aliases[self._tokens[0]])


class BinarySelectionOperand(object):
    """Binary selections

    Examples
    --------
    'name CA' -> lambda atom: atom.name == 'CA'
    'index > 1' -> lambda atom: atom.index > 1
    'resname != TYR' -> lambda atom: atom.residue.name != 'TYR'
    """

    operator_aliases = _kw(
        (['<', 'lt'], ast.Lt()),
        (['==', 'eq'], ast.Eq()),
        (['<=', 'le'], ast.LtE()),
        (['!=', 'ne'], ast.NotEq()),
        (['>=', 'ge'], ast.GtE()),
        (['>', 'gt'], ast.Gt()),
    )

    keyword_aliases = _kw(
        (('name',),             _chained_atom_attr('name',)),
        (('index',),            _chained_atom_attr('index',)),
        (('numbonds',),         _chained_atom_attr('numbonds',)),
        (('type', 'element'),   _chained_atom_attr('element', 'symbol')),
        (('radius',),           _chained_atom_attr('element', 'radius')),
        (('mass',),             _chained_atom_attr('element', 'mass')),
        (('residue', 'resSeq'), _chained_atom_attr('residue', 'resSeq')),
        (('resname',),          _chained_atom_attr('residue', 'name')),
        (('resid',),            _chained_atom_attr('residue', 'index')),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        _check_n_tokens(tokens, 3, 'Binary selectors')
        self.keyword_token, self.op_token, self.comparator_token = tokens
        assert self.keyword_token in self.keyword_aliases
        assert self.op_token in self.operator_aliases

    def ast(self):
        left = self.keyword_aliases[self.keyword_token]
        ops = [self.operator_aliases[self.op_token]]
        comparators = [ast.parse(self.comparator_token, mode='eval').body]
        return deepcopy(ast.Compare(left=left, ops=ops, comparators=comparators))


class UnaryInfixOperand(object):
    n_terms = 1
    assoc = 'RIGHT'

    keyword_aliases = _kw(
        (['not', '!'], ast.Not()),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        _check_n_tokens(tokens, 2, 'Unary infix operators')
        self.op_token, self.value_token = tokens
        assert self.op_token in self.keyword_aliases

    def ast(self):
        return deepcopy(ast.UnaryOp(op=self.keyword_aliases[self.op_token],
                           operand=self.value_token.ast()))


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
            self.op_token = tokes[1]
            self.comparators = tokes[::2]
        else:
            err = "Invalid number of infix expressions: {}"
            err = err.format(len(tokens))
            ParseException = import_('pyparsing').ParseException
            raise ParseException(err)
        assert self.op_token in self.keyword_aliases

    def ast(self):
        return deepcopy(ast.BoolOp(op=self.keyword_aliases[self.op_token],
                                   values=[e.ast() for e in self.comparators]))


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