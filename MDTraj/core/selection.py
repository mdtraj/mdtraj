from pyparsing import (Word, alphas, nums, oneOf, Group, infixNotation, opAssoc,
                       Keyword, MatchFirst, ParseException,
                       Optional)


class Operand:
    _keywords = []

    @classmethod
    def get_keywords(cls):
        return MatchFirst([Keyword(kw) for kw in cls._keywords])


class UnaryOperand(Operand):
    """Single keyword selection identifier."""

    def __init__(self, tokes):
        tokes = tokes[0]
        self.value = tokes
        self.top_type = ''

    def __str__(self):
        return "{top_type}_{value}".format(**self.__dict__)

    __repr__ = __str__


class AtomUnaryOperand(UnaryOperand):
    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'atom'

    _keywords = [
        'all', 'none'
    ]


class ResidueUnaryOperand(UnaryOperand):
    """Single keyword selections that specify residues"""

    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'residue'

    _keywords = [
        'protein', 'nucleic', 'backbone', 'sidechain',
        'water'
    ]


class BinaryOperand(Operand):
    """Selection of the form: field operator value"""

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 3:
            self.key, self.operator, self.value = tokes
        else:
            raise ParseException(
                "Binary selectors take 3 values. You gave {}".format(
                    len(tokes)))

    def __str__(self):
        return "{top_type}_{key} {operator} {value}".format(**self.__dict__)

    __repr__ = __str__


class AtomBinaryOperand(BinaryOperand):
    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'atom'

    _keywords = [
        'name', 'type', 'index', 'numbonds', 'radius', 'mass', 'within'
    ]


class ResidueBinaryOperand(BinaryOperand):
    """Selections that specify residues."""

    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'residue'


    _keywords = [
        'residue', 'resname', 'resid'
    ]


class RangeOperand:
    """Values of the form: first to last"""

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 3:
            self.first, self.operator, self.last = tokes

    def __str__(self):
        return "range({first} {operator} {last})".format(**self.__dict__)

    __repr__ = __str__


class BinaryInfix:
    """Deal with binary infix operators: and, or, etc"""

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) % 2 == 1:
            self.operator = tokes[1]
            self.parts = tokes[::2]
        else:
            raise ParseException(
                "Invalid number of infix expressions: {}".format(len(tokes)))


    def __str__(self):
        middle = " {} ".format(self.operator).join(
            ["{}".format(p) for p in self.parts])
        return "({})".format(middle)

    __repr__ = __str__


class AndInfix(BinaryInfix):
    pass


class OrInfix(BinaryInfix):
    pass


class OfInfix(BinaryInfix):
    pass


class UnaryInfix:
    """Deal with unary infix operators: not"""

    def __init__(self, tokes):
        tokes = tokes[0]
        if len(tokes) == 2:
            self.operator, self.value = tokes
        else:
            raise ParseException(
                "Invalid number of expressions for unary operation: {}".format(
                    len(tokes)))

    def __str__(self):
        return "({} {})".format(self.operator, self.value)

    __repr__ = __str__


class NotInfix(UnaryInfix):
    pass


def _make_parser():
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
    value = Word(alphas) | Word(nums) | numrange

    # - Operators
    comparison_op = oneOf("< == > <= >= != eq gt lt ne")
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
    logical_expr = infixNotation(expression, [
        ('not', 1, opAssoc.RIGHT, NotInfix),
        ('and', 2, opAssoc.LEFT, AndInfix),
        ('or', 2, opAssoc.LEFT, OrInfix),
        ('of', 2, opAssoc.LEFT, OfInfix)
    ])
    return logical_expr