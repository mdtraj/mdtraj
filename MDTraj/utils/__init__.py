from __future__ import print_function, division
from .delay_import import import_
from .arrays import ensure_type
from .unit import in_units_of
from .unitcell import (lengths_and_angles_to_box_vectors,
                       box_vectors_to_lengths_and_angles)

__all__ = ["ensure_type", "import_", "in_units_of",
    "lengths_and_angles_to_box_vectors", "box_vectors_to_lengths_and_angles",
    "ilen"]


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
