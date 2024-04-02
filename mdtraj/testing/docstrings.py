##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
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
##############################################################################


import importlib

##############################################################################
# Imports
##############################################################################
import pkgutil
import warnings
from inspect import (
    getargs,
    getdoc,
    getmembers,
    getmodule,
    isbuiltin,
    isclass,
    isfunction,
    ismethod,
    ismodule,
)

from mdtraj.testing.docscrape import NumpyDocString

__all__ = ["docstring_verifiers", "import_all_modules"]

##############################################################################
# functions
##############################################################################


def docstring_verifiers(module, error_on_none=False):
    """Yield an iterable of tests that verify the docstrings (for
    adherance to the NumpyDoc format) of all functions/classes defined
    in a module.

    Parameters
    ----------
    module : module
        The module to test
    error_on_none : bool, default=False
        Throw an error if no docstring is defined
    """

    # These are the types that we want to check
    # currently, the docstring on classes are not being checked, since
    # isclass() is not in the list
    acceptors = [isfunction, ismethod, isbuiltin]

    def accept(f):
        return any([acc(f) for acc in acceptors])

    functions = [f for f in walk(module) if accept(f)]

    def format(f):
        """
        Format a method/function/class as a string

        Parameters
        ----------
        f : function, method, class

        Returns
        -------
        repr : string
            A string represntation
        """
        if ismethod(f):
            return ".".join([getmodule(f).__name__, f.im_class.__name__, f.__name__])
        if isfunction(f) or isbuiltin(f):
            return ".".join([getmodule(f).__name__, f.__name__])
        if isclass(f):
            return f.__name__
        return "Error"

    def check_docstring(f):
        """
        Ensure the docstring of `f` is in accordance with the numpy standard

        Currently, only the Parameters section of the docstring is checked.

        Parameters
        ----------
        f : function, method, class

        Returns
        -------
        repr : string
            A string represntation
        """

        doc = getdoc(f)
        if doc is None:
            if error_on_none:
                raise ValueError("no docstring for %s" % format(f))
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                parsed = NumpyDocString(doc)

            param_names = {e[0] for e in parsed["Parameters"]}

            if isbuiltin(f):
                # You can't get the arglist from a builtin, which
                # is how cython functions turn up

                # but you can, hackily, get the number of arguments it wants
                # by parseing the error hen you supply too many
                import re

                try:
                    f(*list(range(100)))
                except TypeError as e:
                    m = re.search(r"takes at most (\d+) positional arguments", str(e))
                    if not m:
                        return
                    n_args = int(m.group(1))

                if len(param_names) != n_args:
                    raise ValueError(
                        "In %s, number of arguments, %d, doesn't "
                        " match the length of the Parameters in the "
                        "docstring, %d" % (format(f), n_args, len(param_names)),
                    )
                return

            args = set(getargs(f.__code__).args)

            if "self" in args:
                args.remove("self")
            if "cls" in args:
                args.remove("cls")

            if args != param_names:
                raise ValueError(
                    f"In {format(f)}, arguments {list(args)} don't " f"match Parameters list {list(param_names)}",
                )

    for f in functions:

        def qq():
            return check_docstring(f)

        qq.description = f"NumpyDoc: {module.__name__}.{f.__name__}"
        qq.fname = f.__name__
        yield qq


def ispackage(obj):
    """
    Check if obj is a package. Simply look for whether its a module whose
    filename is __init__.py(c)

    Parameters
    ----------
    obj : module
    """
    if ismodule(obj):
        return obj.__file__.endswith("__init__.pyc") or obj.__file__.endswith(
            "__init__.py",
        )
    return False


def walk(module):
    """
    Get all of the functions, classes and their methods defined within
    a python module

    Parameters
    ----------
    module : module

    Returns
    -------
    funcs : list
        List of functions, classes and methods
    """
    assert ismodule(module)
    if ispackage(module):
        raise ValueError("Sorry, you need to supply me a module, not a package")

    def is_valid(obj):
        if getmodule(obj) == module:
            # cython specific stuff
            if module.__file__.endswith(".so"):
                if isbuiltin(obj):
                    return not obj.__name__.startswith("_")

            if ismethod(obj) or isfunction(obj):
                return not obj.__name__.startswith("_")
            if isclass(obj):
                return True
        return False

    instack = [v for k, v in getmembers(module) if is_valid(v)]
    outstack = []

    while True:
        try:
            item = instack.pop()
        except IndexError:
            break
        outstack.append(item)

        if isclass(item):
            instack.extend([v for k, v in getmembers(item) if is_valid(v)])

    return outstack


def import_all_modules(pkg):
    result = []
    for _, modname, ispkg in pkgutil.iter_modules(pkg.__path__):
        c = f"{pkg.__name__}.{modname}"
        if modname.startswith("test_"):
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=DeprecationWarning)
                mod = importlib.import_module(c)
            if ispkg:
                result.extend(import_all_modules(mod))
            else:
                result.append(mod)
        except ImportError as e:
            print("e", e)
            continue

    return result
