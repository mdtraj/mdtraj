from pyparsing import (Word, alphas, nums, oneOf, Group, infixNotation, opAssoc,
                       Keyword, MatchFirst, ParseException,
                       Optional, quotedString)


class Operand(object):
    keywords = []
    keyword_aliases = {}

    def mdtraj_condition(self):
        """Build an mdtraj-compatible piece of python code."""
        raise NotImplementedError

    @classmethod
    def get_keywords(cls):
        raise NotImplementedError


class SelectionOperand(Operand):
    @classmethod
    def get_keywords(cls):
        return MatchFirst([Keyword(kw) for kw in cls.keywords])


    def get_top_item(self):
        """Get the name of the appropriate topology field for mdtraj.

        E.g.: for atoms it is a.whatever and for residues it
        is a.residue.whatever
        """
        raise NotImplementedError

    def get_top_name(self):
        raise NotImplementedError


class UnaryOperand(SelectionOperand):
    """Single keyword selection identifier."""

    def __init__(self, tokes):
        tokes = tokes[0]
        self.value = tokes

    def __str__(self):
        return "{top_type}_{value}".format(top_type=self.get_top_name(),
                                           value=self.keyword_aliases[
                                               self.value])

    __repr__ = __str__


class AtomUnaryOperand(UnaryOperand):
    keywords = [
        'all', 'none', 'everything', 'nothing'
    ]

    keyword_aliases = {'all': True, 'everything': True,
                       'none': False, 'nothing': False}

    def mdtraj_condition(self):
        return "{}".format(self.keyword_aliases[self.value])

    def get_top_item(self):
        return 'a'

    def get_top_name(self):
        return 'atom'


class ResidueUnaryOperand(UnaryOperand):
    """Single keyword selections that specify residues"""

    keywords = [
        'protein', 'nucleic', 'backbone', 'sidechain',
        'water', 'waters'
    ]

    keyword_aliases = dict([(v, v) for v in keywords])
    keyword_aliases.update(waters='water')

    def mdtraj_condition(self):
        return "a.residue.is_{}".format(self.keyword_aliases[self.value])

    def get_top_item(self):
        return 'a.residue'

    def get_top_name(self):
        return 'residue'


def _quote_value(value):
    """Put quotes around something if it's not a number"""
    try:
        # Is it an int?
        val = int(value)
    except ValueError:
        # Not an int. What about a float?
        try:
            val = float(value)
        except ValueError:
            # Not a float. What about a complex?
            try:
                val = complex(value)
            except ValueError:
                # Not any sort of number. Let's put quotes around it.
                val = value.strip("\"'")
                val = "'{}'".format(val)
    return val


class BinaryOperand(SelectionOperand):
    """Selection of the form: field operator value"""

    operators = ['<', '==', '<=', '!=', '>', '>=', 'eq', 'gt', 'lt', 'ne']
    operator_aliases = dict([(v, v) for v in operators])
    operator_aliases.update(eq='==', gt='>', lt='<', ne='!=')

    # Override in subclasses:
    keyword_aliases = []

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 3:
            self.key, self.in_op, self.value = tokes
        else:
            raise ParseException(
                "Binary selectors take 3 values. You gave {}".format(
                    len(tokes)))

    @property
    def operator(self):
        """Return an unambiguous operator."""
        return self.operator_aliases[self.in_op]

    def __str__(self):
        fmt_dict = dict(top_type=self.get_top_name(),
                        key=self.keyword_aliases[self.key],
                        operator=self.operator, value=self.value)
        return "{top_type}_{key} {operator} {value}".format(**fmt_dict)

    def mdtraj_condition(self):
        field = self.keyword_aliases[self.key]
        if isinstance(self.value, RangeOperand):
            # Special case for dealing with a range
            full_field = "{top}.{field}".format(top=self.get_top_item(),
                                                field=field)
            return self.value.range_condition(full_field, self.operator)

        else:
            fmt_string = "{top}.{field} {op} {value}"
            fmt_dict = dict(top=self.get_top_item(), field=field,
                            op=self.operator,
                            value=_quote_value(str(self.value)))
            return fmt_string.format(**fmt_dict)

    __repr__ = __str__


class AtomBinaryOperand(BinaryOperand):
    keywords = [
        'name', 'type', 'index', 'id', 'numbonds', 'radius', 'mass', 'within',
        'element.symbol', 'element.mass', 'element.name', 'element.number'
    ]
    keyword_aliases = dict([(v, v) for v in keywords])
    keyword_aliases.update(id='index', type='element.symbol',
                           radius='element.radius', mass='element.mass')

    def get_top_name(self):
        return 'atom'

    def get_top_item(self):
        return 'a'


class ResidueBinaryOperand(BinaryOperand):
    """Selections that specify residues."""

    keywords = [
        'residue', 'resname', 'resid'
    ]
    keyword_aliases = dict([(v, v) for v in keywords])
    keyword_aliases.update(residue='resSeq', resid='index',
                           resname='name')  # TODO: Are these correct?

    def get_top_name(self):
        return 'residue'

    def get_top_item(self):
        return 'a.residue'


class RangeOperand:
    """Values of the form: first to last"""

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 3:
            self.first, self.to, self.last = tokes

    def __str__(self):
        return "range({first} {to} {last})".format(**self.__dict__)

    __repr__ = __str__

    def range_condition(self, field, op):
        if self.to == 'to' and op == '==':
            pass
            return "{first} <= {field} <= {last}".format(
                first=self.first, field=field, last=self.last
            )
        else:
            # We may want to be able to do more fancy things later on
            # For example: "mass > 5 to 10" could be parsed (even though
            # it's kinda stupid)
            raise ParseException("Incorect use of ranged value")


class InfixOperand(Operand):
    # Overwrite the following in base classes.
    num_terms = -1
    assoc = None

    @classmethod
    def get_keywords(cls):
        # Prepare tuples for pyparsing.infixNotation
        return [(kw, cls.num_terms, cls.assoc, cls) for kw in cls.keywords]


class BinaryInfix(InfixOperand):
    """Deal with binary infix operators: and, or, etc"""

    num_terms = 2
    assoc = opAssoc.LEFT

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) % 2 == 1:
            self.in_op = tokes[1]
            self.parts = tokes[::2]
        else:
            raise ParseException(
                "Invalid number of infix expressions: {}".format(len(tokes)))


    def __str__(self):
        # Join all the parts
        middle = " {} ".format(self.operator).join(
            ["{}".format(p) for p in self.parts])

        # And put parenthesis around it
        return "({})".format(middle)

    __repr__ = __str__

    def mdtraj_condition(self):

        # Join all the parts
        middle = " {} ".format(self.operator).join(
            [p.mdtraj_condition() for p in self.parts])

        # Put parenthesis around it
        return "({})".format(middle)

    @property
    def operator(self):
        return self.keyword_aliases[self.in_op]


class AndInfix(BinaryInfix):
    keywords = ['and', '&&']
    keyword_aliases = dict([(kw, 'and') for kw in keywords])


class OrInfix(BinaryInfix):
    keywords = ['or', '||']
    keyword_aliases = dict([(kw, 'or') for kw in keywords])


class OfInfix(BinaryInfix):
    keywords = ['of']
    keyword_aliases = dict([(kw, 'of') for kw in keywords])


class NotInfix(InfixOperand):
    """Deal with unary infix operators: not"""

    num_terms = 1
    assoc = opAssoc.RIGHT

    keywords = ['not', '!']
    keyword_aliases = dict([(kw, 'not') for kw in keywords])

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 2:
            self.in_op, self.value = tokes
        else:
            raise ParseException(
                "Invalid number of expressions for unary operation: {}".format(
                    len(tokes)))

    def __str__(self):
        return "({} {})".format(self.operator, self.value)

    __repr__ = __str__

    def mdtraj_condition(self):
        return "({} {})".format(self.operator, self.value.mdtraj_condition())

    @property
    def operator(self):
        return self.keyword_aliases[self.in_op]


class SelectionParser(object):
    def __init__(self, select_string=None):
        self.parser = _make_parser()

        if select_string is not None:
            self.last_parse = self.parser.parseString(select_string)[0]
        else:
            # Default to all atoms
            self.last_parse = self.parser.parseString("all")[0]


    def parse(self, select_string):
        """Parse a selection string."""

        # We need to select element zero of the result of parseString
        # I don't know what would ever go in the rest of the list.
        self.last_parse = self.parser.parseString(select_string)[0]
        return self

    @property
    def mdtraj_condition(self):
        return self.last_parse.mdtraj_condition()

    @property
    def unambiguous(self):
        return str(self.last_parse)


def _make_parser():
    """Build an expression that can parse atom selection strings."""

    # Single Keyword
    # - Atom
    single_atom = AtomUnaryOperand.get_keywords()
    single_atom.setParseAction(AtomUnaryOperand)
    # - Residue
    single_resi = ResidueUnaryOperand.get_keywords()
    single_resi.setParseAction(ResidueUnaryOperand)

    # Key Value Keywords
    # - A value is a number, word, or range
    numrange = Group(Word(nums) + "to" + Word(nums))
    numrange.setParseAction(RangeOperand)
    value = numrange | Word(nums) | quotedString | Word(alphas)

    # - Operators
    comparison_op = oneOf(BinaryOperand.operators)
    comparison_op = Optional(comparison_op, '==')

    # - Atom
    binary_atom = Group(
        AtomBinaryOperand.get_keywords() + comparison_op + value)
    binary_atom.setParseAction(AtomBinaryOperand)
    # - Residue
    binary_resi = Group(
        ResidueBinaryOperand.get_keywords() + comparison_op + value)
    binary_resi.setParseAction(ResidueBinaryOperand)

    # Put it together
    expression = MatchFirst([
        single_atom, single_resi, binary_atom, binary_resi
    ])

    # And deal with logical expressions
    logical_expr = infixNotation(expression, (
        AndInfix.get_keywords() + OrInfix.get_keywords() +
        NotInfix.get_keywords() + OfInfix.get_keywords()
    ))
    return logical_expr