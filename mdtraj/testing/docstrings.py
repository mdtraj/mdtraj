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


##############################################################################
# Imports
##############################################################################
from __future__ import print_function, division
import sys
import types
import warnings
from inspect import (isclass, ismodule, isfunction, ismethod,
                     getmembers, getdoc, getmodule, getargs, isbuiltin)
from mdtraj.testing.docscrape import NumpyDocString
from mdtraj.utils.six import get_function_code, get_function_closure, PY2

__all__ = ['DocStringFormatTester']

##############################################################################
# functions
##############################################################################


def DocStringFormatTester(module, error_on_none=False):
    """
    Create a class that tests the docstrings of a python module for adherance
    to the numpy docstring format, for use with Nose.

    Parameters
    ----------
    module : module
        The module to test
    error_on_none : bool
        Throw an error if no docstring is defined

    Notes
    -----
    For example, test_trajectory.py in mdtraj contains the line:

    TestDocstrings = DocStringFormatTester(mdtraj.trajectory)

    TestDocstrings is now a class with one method defined per class/function in
    mdtraj.trajectory. Each method, when called, validates a single docstring.
    When nosetests runs that file (in verbose mode), you see:

    NumpyDoc: mdtraj.trajectory.load_xtc ... ok
    NumpyDoc: mdtraj.trajectory.load_pdb ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save_binpos ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save ... ok
    NumpyDoc: mdtraj.trajectory.load ... ok
    NumpyDoc: mdtraj.trajectory.load_hdf ... ok
    NumpyDoc: mdtraj.trajectory.load_dcd ... ok
    NumpyDoc: mdtraj.trajectory.load_binpos ... ok
    NumpyDoc: mdtraj.trajectory.load ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save_xtc ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save_pdb ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save_hdf ... ok
    NumpyDoc: mdtraj.trajectory.Trajectory.save_dcd ... ok
    """

    # These are the types that we want to check
    # currently, the docstring on classes are not being checked, since
    # isclass() is not in the list
    acceptors = [isfunction, ismethod, isbuiltin]
    accept = lambda f: any([acc(f) for acc in acceptors])
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
            if PY2:
                return '.'.join([getmodule(f).__name__, f.im_class.__name__, f.__name__])
            else:
                return '.'.join([getmodule(f).__name__, f.__self__.__class__.__name__, f.__name__])
        if isfunction(f) or isbuiltin(f):
            return '.'.join([getmodule(f).__name__, f.__name__])
        if isclass(f):
            return f.__name__
        return 'Error'

    def check_docstring(self, f):
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
                raise ValueError('no docstring for %s' % format(f))
        else:
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                parsed = NumpyDocString(doc)

            param_names = set([e[0] for e in parsed['Parameters']])

            if isbuiltin(f):
                # You can't get the arglist from a builtin, which
                # is how cython functions turn up

                # but you can, hackily, get the number of arguments it wants
                # by parseing the error hen you supply too many
                import re
                try:
                    f(*list(range(100)))
                except TypeError as e:
                    m = re.search('takes at most (\d+) positional arguments', str(e))
                    if not m:
                        return
                    n_args = int(m.group(1))

                if len(param_names) != n_args:
                    raise ValueError("In %s, number of arguments, %d, doesn't "
                        " match the length of the Parameters in the "
                        "docstring, %d" % (format(f), n_args, len(param_names)))
                return

            args = set(getargs(get_function_code(f)).args)
            if 'self' in args:
                args.remove('self')

            if args != param_names:
                raise ValueError("In %s, arguments %s don't "
                    "match Parameters list %s" % (format(f),
                        list(args), list(param_names)))

    funcdict = {}
    # populate the func dict before calling type()
    for i, f in enumerate(functions):
        name = 'test_%s' % i  # this is the name we give the method

        # create the method. this is a little complicated, but the basic
        # idea is that NoseTests checks the func_name, so we need to create
        # the method correctly. Just giving it a pointer to a closure we define
        # won't set func_name correctly. Instead, we make a function where
        # func_code is set to the func_code of check_docstring, but we add
        # set second argument of that function (after self) to be a default
        # arg -- the f that we're iterating over.
        method = types.FunctionType(get_function_code(check_docstring), globals(), name,
            (f,), get_function_closure(check_docstring))

        # give the method a short docstring
        if PY2:
            method.func_doc = 'NumpyDoc: ' + format(f)
        else:
            method.__doc__ = 'NumpyDoc: ' + format(f)
        funcdict[name] = method

    return type('TestDoc', (), funcdict)


def ispackage(obj):
    """
    Check if obj is a package. Simply look for whether its a module whose
    filename is __init__.py(c)

    Parameters
    ----------
    obj : module
    """
    if ismodule(obj):
        return obj.__file__.endswith("__init__.pyc") or \
            obj.__file__.endswith("__init__.py")
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
        raise ValueError('Sorry, you need to supply me a module, not a package')

    def is_valid(obj):
        if getmodule(obj) == module:
            # cython specific stuff
            if module.__file__.endswith('.so'):
                if isbuiltin(obj):
                    return not obj.__name__.startswith('_')

            if ismethod(obj) or isfunction(obj):
                return not obj.__name__.startswith('_')
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


TestModule = DocStringFormatTester(sys.modules[__name__])
