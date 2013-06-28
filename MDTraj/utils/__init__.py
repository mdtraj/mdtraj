from delay_import import import_
from arrays import ensure_type
from unit import in_units_of


def ilen(iterable):
    """Length of an iterator. Note, this consumes the iterator

    Parameters
    ----------
    iterable : iterable
        An iterable, such as a generator, list, etc.

    Returns
    -------
    length  : int
        The number of elements in the iterable
    """
    return sum(1 for _ in iterable)
