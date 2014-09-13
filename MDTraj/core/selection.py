from pyparsing import (Word, alphas, nums, oneOf, Group, infixNotation, opAssoc,
                       Keyword, MatchFirst, ParseException,
                       Optional)


class UnaryOperand:
    """Single keyword selection identifier."""

    def __init__(self, tokes):
        tokes = tokes[0]
        self.value = tokes

    def __str__(self):
        return "{top_type}_{value}".format(**self.__dict__)

    __repr__ = __str__


class ResidueUnaryOperand(UnaryOperand):
    """Single keyword selections that specify residues"""

    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'residue'


class BinaryOperand:
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


class ResidueBinaryOperand(BinaryOperand):
    """Selections that specify residues."""

    def __init__(self, tokes):
        super().__init__(tokes)
        self.top_type = 'residue'


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


def _get_parser():
    # Support for numeric ranges
    numrange = Group(Word(nums) + "to" + Word(nums))
    numrange.setParseAction(RangeOperand)
    numeric = numrange | Word(nums)

    # Single keyword selections
    single_keyword = MatchFirst([
        Keyword('protein'), Keyword('solvent'), Keyword('water')
    ])
    single_keyword.setParseAction(ResidueUnaryOperand)

    # Key value keywords
    residue_keyword = MatchFirst([
        Keyword('resname'), Keyword('resid'), Keyword('within')
    ])

    # Comparison operators
    comparison_op = oneOf("< == > <= >= != eq gt lt ne")
    comparison_expr = Group(residue_keyword + Optional(comparison_op, '==') + (
        Word(alphas) | numeric))
    comparison_expr.setParseAction(ResidueBinaryOperand)

    tot_expr = single_keyword | comparison_expr

    logicalExpr = infixNotation(tot_expr,
                                [
                                    ('not', 1, opAssoc.RIGHT, NotInfix),
                                    ('and', 2, opAssoc.LEFT, AndInfix),
                                    ('or', 2, opAssoc.LEFT, OrInfix),
                                    ('of', 2, opAssoc.LEFT, OfInfix)
                                ])
    return logicalExpr