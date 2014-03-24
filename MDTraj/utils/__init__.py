from __future__ import print_function, division
import time
from mdtraj.utils.delay_import import import_
from mdtraj.utils.validation import ensure_type, cast_indices, check_random_state
from mdtraj.utils.unit import in_units_of
from mdtraj.utils.rotation import rotation_matrix_from_quaternion, uniform_quaternion
from mdtraj.utils.unitcell import (lengths_and_angles_to_box_vectors,
                       box_vectors_to_lengths_and_angles)
from mdtraj.utils.contextmanagers import timing, enter_temp_directory

__all__ = ["ensure_type", "import_", "in_units_of",
           "lengths_and_angles_to_box_vectors",
           "box_vectors_to_lengths_and_angles",
           "ilen", "timing", "cast_indices", "check_random_state",
           "rotation_matrix_from_quaternion", "uniform_quaternion",
           "enter_temp_directory", "timing"]


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

