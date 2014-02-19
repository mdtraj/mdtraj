from __future__ import print_function, division
import time
from mdtraj.utils.delay_import import import_
from mdtraj.utils.arrays import ensure_type, cast_indices
from mdtraj.utils.unit import in_units_of, convert
from mdtraj.utils.unitcell import (lengths_and_angles_to_box_vectors,
                       box_vectors_to_lengths_and_angles)

__all__ = ["ensure_type", "import_", "in_units_of",
    "lengths_and_angles_to_box_vectors", "box_vectors_to_lengths_and_angles",
    "ilen", "timing", "cast_indices", "convert"]


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
    >>> long_function = lambda : None
    >>> with timing('long_function'):
    ...     long_function()
    long_function: 0.000 seconds
    """
    def __init__(self, name='block'):
        self.name = name
        self.time = 0
        self.start = None
        self.end = None
    
    def __enter__(self):
        self.start = time.time()
        return self
    
    def __exit__(self, ty, val, tb):
        self.end = time.time()
        self.time = self.end - self.start
        print("%s: %0.3f seconds" % (self.name, self.time))
        return False
