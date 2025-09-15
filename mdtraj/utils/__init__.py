import time
import warnings

from mdtraj.utils.contextmanagers import enter_temp_directory, timing
from mdtraj.utils.delay_import import import_
from mdtraj.utils.rotation import rotation_matrix_from_quaternion, uniform_quaternion
from mdtraj.utils.unit import in_units_of
from mdtraj.utils.unitcell import (
    box_vectors_to_lengths_and_angles,
    lengths_and_angles_to_box_vectors,
    lengths_and_angles_to_tilt_factors,
)
from mdtraj.utils.validation import cast_indices, check_random_state, ensure_type
from mdtraj.utils.zipped import open_maybe_zipped

__all__ = [
    "ensure_type",
    "import_",
    "in_units_of",
    "lengths_and_angles_to_box_vectors",
    "box_vectors_to_lengths_and_angles",
    "tilt_factors_to_angles",
    "ilen",
    "timing",
    "cast_indices",
    "check_random_state",
    "rotation_matrix_from_quaternion",
    "uniform_quaternion",
    "enter_temp_directory",
    "timing",
    "deprecated",
]


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


class deprecated:
    """Decorator to mark a function or class as deprecated.

    Issue a warning when the function is called/the class is instantiated and
    adds a warning to the docstring.

    The optional extra argument will be appended to the deprecation message
    and the docstring. Note: to use this with the default value for extra, put
    in an empty of parentheses:

    >>> from sklearn.utils import deprecated
    >>> deprecated() # doctest: +ELLIPSIS
    <sklearn.utils.deprecated object at ...>

    >>> @deprecated()
    ... def some_function(): pass
    """

    # Copied from scikit-learn: sklearn/utils/__init__.py
    # Adapted from http://wiki.python.org/moin/PythonDecoratorLibrary,
    # but with many changes.

    def __init__(self, extra=""):
        """
        Parameters
        ----------
        extra: string
          to be added to the deprecation messages

        """
        self.extra = extra

    def __call__(self, obj):
        if isinstance(obj, type):
            return self._decorate_class(obj)
        else:
            return self._decorate_fun(obj)

    def _decorate_class(self, cls):
        msg = "Class %s is deprecated" % cls.__name__
        if self.extra:
            msg += "; %s" % self.extra

        # FIXME: we should probably reset __new__ for full generality
        init = cls.__init__

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return init(*args, **kwargs)

        cls.__init__ = wrapped

        wrapped.__name__ = "__init__"
        wrapped.__doc__ = self._update_doc(init.__doc__)
        wrapped.deprecated_original = init

        return cls

    def _decorate_fun(self, fun):
        """Decorate function fun"""

        msg = "Function %s is deprecated" % fun.__name__
        if self.extra:
            msg += "; %s" % self.extra

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return fun(*args, **kwargs)

        wrapped.__name__ = fun.__name__
        wrapped.__dict__ = fun.__dict__
        wrapped.__doc__ = self._update_doc(fun.__doc__)

        return wrapped

    def _update_doc(self, olddoc):
        newdoc = "DEPRECATED"
        if self.extra:
            newdoc = f"{newdoc}: {self.extra}"
        if olddoc:
            newdoc = f"{newdoc}\n\n{olddoc}"
        return newdoc
