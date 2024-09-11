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
import re
import sys
from ast import unparse
from collections import namedtuple
from copy import deepcopy

from pyparsing import (
    Group,
    Keyword,
    MatchFirst,
    OneOrMore,
    ParseException,
    ParserElement,
    Word,
    alphanums,
    alphas,
    infixNotation,
    opAssoc,
    quotedString,
)


# this number arises from the current selection language, if the cache size is exceeded, it hurts performance a bit.
ParserElement.enablePackrat(cache_size_limit=304)

__all__ = ["parse_selection"]

# ############################################################################
# Globals
# ############################################################################

NUMS = ".0123456789"
THIS_ATOM = ast.Name(id="atom", ctx=ast.Load(), SINGLETON=True)
RE_MODULE = ast.Name(id="re", ctx=ast.Load(), SINGLETON=True)
SELECTION_GLOBALS = {"re": re}
_ParsedSelection = namedtuple("_ParsedSelection", ["expr", "source", "astnode"])

# ############################################################################
# Utils
# ############################################################################


class _RewriteNames(ast.NodeTransformer):
    def visit_Name(self, node):
        if hasattr(node, "SINGLETON"):
            return node

        _safe_names = {"None": None, "True": True, "False": False}
        if node.id in _safe_names:
            return ast.Constant(value=_safe_names[node.id], kind=None)

        # all other bare names are taken to be string literals. Thus something
        # like parse_selection('name CA') properly resolves CA as a string
        # literal, not a barename to be loaded from the global scope!
        return ast.Constant(value=node.id, kind=None)


def _chain(*attrs):
    """This transforms, for example, ('residue', 'is_protein'), into
    Attribute(value=Attribute(value=THIS_ATOM,
     attr='residue', ctx=Load()), attr='is_protein', ctx=Load())
    """
    left = THIS_ATOM
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
        err = "{} take {} values. You gave {}"
        err = err.format(name, n_tokens, len(tokens))
        raise ParseException(err)


class SelectionKeyword:
    keyword_aliases = _kw(
        # Atom.<attribute>
        (("all", "everything"), ast.Name(id="True", ctx=ast.Load())),
        (("none", "nothing"), ast.Name(id="False", ctx=ast.Load())),
        (("backbone", "is_backbone"), _chain("is_backbone")),
        (("sidechain", "is_sidechain"), _chain("is_sidechain")),
        # Atom.residue.<attribute>
        (("protein", "is_protein"), _chain("residue", "is_protein")),
        (("code", "rescode", "resc"), _chain("residue", "code")),
        # (('nucleic', 'is_nucleic'), _chain('residue', 'is_nucleic')),
        (("water", "waters", "is_water"), _chain("residue", "is_water")),
        (("name",), _chain("name")),
        (("index",), _chain("index")),
        (("n_bonds",), _chain("n_bonds")),
        (("residue", "resSeq"), _chain("residue", "resSeq")),
        (("resname", "resn"), _chain("residue", "name")),
        (("resid", "resi"), _chain("residue", "index")),
        (("segment_id", "segname"), _chain("segment_id")),
        # Atom.residue.chain.<attribute>
        (("chainid",), _chain("residue", "chain", "index")),
        # Atom.element.<attribute>
        (("type", "element", "symbol"), _chain("element", "symbol")),
        # (('radius',), _chain('element', 'radius')),
        (("mass",), _chain("element", "mass")),
    )

    def __init__(self, tokens):
        # pyparsing constructs the instance while building the parse tree,
        # and gives us the set tokens. In this case, the tokens are
        self._tokens = tokens
        _check_n_tokens(tokens, 1, "Unary selectors")
        assert tokens[0] in self.keyword_aliases

    def ast(self):
        return self.keyword_aliases[self._tokens[0]]


class Literal:
    def __init__(self, tokens):
        self.token = tokens[0]
        _check_n_tokens(tokens, 1, "literal")

    def ast(self):
        return ast.parse(self.token, mode="eval").body


class UnaryInfixOperand:
    n_terms = 1
    assoc = "RIGHT"

    keyword_aliases = _kw(
        (["not ", "!"], ast.Not()),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        _check_n_tokens(tokens, 2, "Unary infix operators")
        self.op_token, self.value_token = tokens
        assert self.op_token in self.keyword_aliases
        if isinstance(self.value_token, Literal):
            raise ValueError("Cannot use literals as booleans.")

    def ast(self):
        return ast.UnaryOp(
            op=self.keyword_aliases[self.op_token],
            operand=self.value_token.ast(),
        )


class RegexInfixOperand:
    n_terms = 2
    assoc = "LEFT"
    keyword_aliases = {"=~": "=~"}

    def __init__(self, tokens):
        self.tokens = tokens[0]
        _check_n_tokens(self.tokens, 3, "regex operator")
        self.string, op, self.pattern = self.tokens
        assert op == "=~"
        if isinstance(self.string, Literal):
            raise ValueError("Cannot do regex comparison on literal")

    def ast(self):
        pattern = self.tokens[2].ast()
        string = self.tokens[0].ast()
        return ast.Compare(
            left=ast.Call(
                func=ast.Attribute(
                    value=RE_MODULE,
                    attr="match",
                    ctx=ast.Load(),
                ),
                args=[pattern, string],
                keywords=[],
                starargs=None,
                kwargs=None,
            ),
            ops=[ast.IsNot()],
            comparators=[ast.Name(id="None", ctx=ast.Load())],
        )


class BinaryInfixOperand:
    n_terms = 2
    assoc = "LEFT"

    keyword_aliases = _kw(
        (["and", "&&"], ast.And()),
        (["or", "||"], ast.Or()),
        (["<", "lt"], ast.Lt()),
        (["==", "eq"], ast.Eq()),
        (["<=", "le"], ast.LtE()),
        (["!=", "ne"], ast.NotEq()),
        ([">=", "ge"], ast.GtE()),
        ([">", "gt"], ast.Gt()),
    )

    def __init__(self, tokens):
        tokens = tokens[0]
        if len(tokens) % 2 == 1:
            self.op_token = tokens[1]
            self.comparators = tokens[::2]
        else:
            err = "Invalid number of infix expressions: {}"
            err = err.format(len(tokens))
            raise ParseException(err)
        assert self.op_token in self.keyword_aliases

        # Check for too many literals and not enough keywords
        op = self.keyword_aliases[self.op_token]
        if isinstance(op, ast.boolop):
            if any(isinstance(c, Literal) for c in self.comparators):
                raise ValueError("Cannot use literals as truth")
        else:
            if all(isinstance(c, Literal) for c in self.comparators):
                raise ValueError("Cannot compare literals.")

    def ast(self):
        op = self.keyword_aliases[self.op_token]

        if isinstance(op, ast.boolop):
            # and and or use one type of AST node
            value = ast.BoolOp(op=op, values=[e.ast() for e in self.comparators])
        else:
            # remaining operators use another
            value = ast.Compare(
                left=self.comparators[0].ast(),
                ops=[op],
                comparators=[e.ast() for e in self.comparators[1:]],
            )
        return value


class RangeCondition:
    def __init__(self, tokens):
        tokens = tokens[0]
        _check_n_tokens(tokens, 4, "range condition")
        assert tokens[2] == "to"
        self._field, self._from, self._to = tokens[0], tokens[1], tokens[3]
        if isinstance(self._field, Literal):
            raise ValueError("Can't test literal in range.")

    def ast(self):
        return ast.Compare(
            left=self._from.ast(),
            ops=[ast.LtE(), ast.LtE()],
            comparators=[self._field.ast(), self._to.ast()],
        )


class InListCondition:
    def __init__(self, tokens):
        tokens = tokens[0]
        self._field = tokens[0]

        if len(tokens) == 2:
            # Implicit equality
            self.implicit_equality = True
        elif len(tokens) > 2:
            self.implicit_equality = False
        else:
            raise ValueError("Not enough tokens for `in` condition")

        self.compare_to = tokens[1:]

    def ast(self):
        if self.implicit_equality:
            return ast.Compare(
                left=self._field.ast(),
                ops=[ast.Eq()],
                comparators=[e.ast() for e in self.compare_to],
            )
        else:
            comparator = ast.List([e.ast() for e in self.compare_to], ast.Load())
            return ast.Compare(
                left=self._field.ast(),
                ops=[ast.In()],
                comparators=[comparator],
            )


class parse_selection:
    """Parse an atom selection expression

    Parameters
    ----------
    selection_string : str
        Selection string, a string in the MDTraj atom selection grammar.

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
        def keywords(klass):
            kws = sorted(klass.keyword_aliases.keys())
            return MatchFirst([Keyword(kw) for kw in kws])

        def infix(klass):
            kws = sorted(klass.keyword_aliases.keys())
            return [(kw, klass.n_terms, getattr(opAssoc, klass.assoc), klass) for kw in kws]

        # literals include words made of alphanumerics, numbers,
        # or quoted strings but we exclude any of the logical
        # operands (e.g. 'or') from being parsed literals
        literal = ~(keywords(BinaryInfixOperand) | keywords(UnaryInfixOperand)) + (
            Word(NUMS) | quotedString | Word(alphas, alphanums)
        )
        literal.setParseAction(Literal)

        # These are the other 'root' expressions,
        # the selection keywords (resname, resid, mass, etc)
        selection_keyword = keywords(SelectionKeyword)
        selection_keyword.setParseAction(SelectionKeyword)
        base_expression = MatchFirst([selection_keyword, literal])

        # range condition matches expressions such as 'mass 1 to 20'
        range_condition = Group(
            selection_keyword + literal + Keyword("to") + literal,
        )
        range_condition.setParseAction(RangeCondition)

        # matches expression such as `resname GLU ASP ARG` and also
        # handles implicit equality `resname ALA`
        in_list_condition = Group(selection_keyword + OneOrMore(literal))
        in_list_condition.setParseAction(InListCondition)

        expression = range_condition | in_list_condition | base_expression
        logical_expr = infixNotation(
            expression,
            infix(UnaryInfixOperand) + infix(BinaryInfixOperand) + infix(RegexInfixOperand),
        )

        self.expression = logical_expr
        self.is_initialized = True

        self.transformer = _RewriteNames()

    def __call__(self, selection):
        if not self.is_initialized:
            self._initialize()

        try:
            parse_result = self.expression.parseString(selection, parseAll=True)
        except ParseException as e:
            msg = str(e)
            lines = [
                f"{msg}: {selection}",
                " " * (12 + len("%s: " % msg) + e.loc) + "^^^",
            ]
            raise ValueError("\n".join(lines))

        # Change __ATOM__ in function bodies. It must bind to the arg
        # name specified below (i.e. 'atom')
        astnode = self.transformer.visit(deepcopy(parse_result[0].ast()))

        # Special check for a single literal

        if isinstance(astnode, ast.Constant) and astnode.value not in {True, False, None}:
            raise ValueError(
                "Cannot use a single literal as a boolean. " f"Choked on node with value {astnode.value}",
            )

        args = [ast.arg(arg="atom", annotation=None)]
        signature = ast.arguments(
            args=args,
            vararg=None,
            kwarg=None,
            posonlyargs=[],
            kwonlyargs=[],
            defaults=[],
            kw_defaults=[],
        )

        func = ast.Expression(body=ast.Lambda(signature, astnode))
        source = unparse(astnode)

        expr = eval(
            compile(ast.fix_missing_locations(func), "<string>", mode="eval"),
            SELECTION_GLOBALS,
        )
        return _ParsedSelection(expr, source, astnode)


# Create the callable, and use it to overshadow the class. this way there's
# basically just one global instance of the "function", even thought its
# a callable class.
parse_selection = parse_selection()

if __name__ == "__main__":
    exp = parse_selection(sys.argv[1])
    print(exp.source)
    print(ast.dump(exp.astnode))
