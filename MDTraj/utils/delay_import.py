"""
Code to delay the import of a moldule, and give a nice error message if
the module is not installed. for dealing with dependencies.
"""
##############################################################################
# imports
##############################################################################

import os
import sys
import inspect
import importlib
import textwrap

__all__ = ['import_']

##############################################################################
# globals
##############################################################################

MESSAGES = {
    'networkx': '''
    The code at {filename}:{line_number} requires the python module
    NetworkX, which is a software package for the creation, manipulation, and study of
    the structure, dynamics, and functions of complex networks.

    NetworkX can be downloaded from https://pypi.python.org/pypi/networkx/, or
    installed with the python "easy_install" or "pip" package managers using:

    easy_install networkx
    or
    pip install networkx
    ''',

    'tables': '''
    The code at {filename}:{line_number} requires the python module PyTables,
    which is a package for managing hierarchical datasets and designed to
    efficiently and easily cope with extremely large amounts of data.

    PyTables can be downloaded from http://www.pytables.org, or installed with
    the python "easy_install" or "pip" package managers using:

    easy_install tables
    or
    pip install tables

    PyTables also depends on the numexpr package, as well as the C-language
    HDF5 library. For detailed installation instructions, visit
    http://pytables.github.io/usersguide/installation.html
    ''',

    'netCDF4': '''
    The code at {filename}:{line_number} requires the netcdf4-python module,
    which is a python interface to the NetCDF software libraries and self-describing,
    machine-independent data formats that support the creation, access, and
    sharing of array-oriented scientific data.

    netcdf4-python can be downloaded from https://pypi.python.org/pypi/netCDF,
    or installed with the python "easy_install" or "pip" package managers using:

    easy_install netCDF4
    or
    pip install netCDF4

    netcdf4-python also depends on the C-language HDF5 and NetCDF libraries.
    For detailed installation instructions, visit
    http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    ''',

    'simtk.unit': '''
    The code at {filename}:{line_number} requires the simtk.unit module,
    which is a python package for unit conversion.

    simtk.unit is installed with OpenMM, which is available at http://openmm.org
    It's also installable as a separate standalone package from
    https://github.com/rmcgibbo/simtk.unit, and can be installed with the python
    "pip" package mangers using:

    pip install git+git://github.com/rmcgibbo/simtk.unit
    '''
}

##############################################################################
# functions
##############################################################################

def import_(module):
    """Import a module, and issue a nice message to stderr if the module isn't installed.

    Currently, this function will print nice error messages for networkx,
    tables, netCDF4, and simtk.unit, which are optional MDTraj dependencies.

    Parameters
    ----------
    module : str
        The module you'd like to import, as a string

    Returns
    -------
    module : {module, object}
        The module object

    Examples
    --------
    >>> # the following two lines are equivalent. the difference is that the
    >>> # second will check for an ImportError and print you a very nice
    >>> # user-facing message about what's wrong (where you can install the
    >>> # module from, etc) if the import fails
    >>> import tables
    >>> tables = import_('tables')
    """
    try:
        return importlib.import_module(module)
    except ImportError:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = 'The code at {filename}:{line_number} requires the ' + module + ' package'

        frame,filename,line_number,function_name,lines,index = \
            inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(filename=filename, line_number=line_number)
        m = textwrap.dedent(m)

        bar = '\033[91m' + '#' * max(len(line) for line in m.split(os.linesep)) + '\033[0m'

        print >> sys.stderr, bar
        print >> sys.stderr, m
        print >> sys.stderr, bar
        sys.exit(1)
