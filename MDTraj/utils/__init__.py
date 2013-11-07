from __future__ import print_function, division
import time
from .delay_import import import_
from .arrays import ensure_type
from .unit import in_units_of
from .unitcell import (lengths_and_angles_to_box_vectors,
                       box_vectors_to_lengths_and_angles)

__all__ = ["ensure_type", "import_", "in_units_of",
    "lengths_and_angles_to_box_vectors", "box_vectors_to_lengths_and_angles",
    "ilen", "timing"]


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


class timing(object):
    """A timing context manager

    Example
    -------
    >>> with timing('long function'):
    >>>    long_function()
    long function: 0.500 seconds
    """
    def __init__(self, name='block'):
        self.name = name
    
    def __enter__(self):
        self.start = time.time()
    
    def __exit__(self, ty, val, tb):
        end = time.time()
        print("%s : %0.3f seconds" % (self.name, end - self.start))
        return False
