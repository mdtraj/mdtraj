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


# Portions of this code were copied from the Numpy source
# those portions are owned by the numpy developers, and released
# under the following license:

# Copyright 2005-2012, NumPy Developers.
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# Neither the name of the NumPy Developers nor the names of any contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

"""
mdtraj.io provides IO for generic arrays via the functions ``mdtraj.io.loadh``
and ``mdtraj.io.saveh``, for loading and saving generic arrays to disk.
These functions act like ``numpy.savez`` and ``numpy.loadz``, but use a
PyTables HDF5 for superior performance and compression.

It also provides open_maybe_zipped for selecting the correct reader based on the
filename.

Nothing in this module is specific to molecular dynamics trajectories.

Examples
--------
>>> import numpy as np
>>> from mdtraj import io
>>> x = np.random.randn(100, 3)
>>> io.saveh('file.hdf5', x=x)
>>> np.all(x == io.loadh('file.hdf5')['x'])
True

Functions
---------
"""

import os
import warnings

import numpy as np

from mdtraj.utils import import_

tables = import_("tables")

tables = import_("tables")
TABLES2 = tables.__version__ < "3.0.0"

__all__ = ["saveh", "loadh"]

try:
    COMPRESSION = tables.Filters(complevel=9, complib="zlib", shuffle=True)
except Exception:  # type?
    warnings.warn("Missing Zlib; no compression will used.")
    COMPRESSION = tables.Filters()


# Note to developers: This module is pseudo-deprecated. It provides (loadh, saveh)
# which are useful functions (and we want to maintain them), but aren't really
# within the scope of MDTraj as we now understand it.
#
# With that said, many people use these functions and no good would come from getting
# rid of them. But we shouldn't add any new functions or new features to this file.
#
# One potential landmine is that this file _requires_ the `tables` package, which is
# only an _optional_ dependency in MDTraj. So if you import this file (or anything in
# it) from another file that is imported at startup (on a user running `import mdtraj`)
# then tables ceases to be an optional dependency and becomes a strict requirement.
#
# So add new features to a different file, and there shouldn't be any reason for
# any files inside MDTraj to `import mdtraj.io`.
#
# See github issue #852.


def saveh(file, *args, **kwargs):
    """Save several numpy arrays into a single file in compressed ``.hdf`` format.

    If arguments are passed in with no keywords, the corresponding variable
    names, in the ``.hdf`` file, are 'arr_0', 'arr_1', etc. If keyword arguments
    are given, the corresponding variable names, in the ``.hdf`` file will
    match the keyword names.

    Parameters
    ----------
    file : str or tables.File
        Either the file name (string) or an open pytables file
        (file-like object opened with tables.openFile(...))
        where the data will be saved.
    args : Arguments, optional
        Arrays to save to the file. Since it is not possible for Python to
        know the names of the arrays outside `savez`, the arrays will be saved
        with names "arr_0", "arr_1", and so on. These arguments can be any
        expression.
    kwds : Keyword arguments, optional
        Arrays to save to the file. Arrays will be saved in the file with the
        keyword names.

    Notes
    -----
    `saveh` will overwrite files by default. If you have an hdf5 that contains the
    arrays `arr_0` and `arr_1` and you attempt to save a new array `x`, it will
    go in side by side. But if you save a new `arr_0`, it will overwrite your
    previous array.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        When arrays are of an unsupported type

    See Also
    --------
    numpy.savez : Save several arrays into a single file in uncompressed ``.npz`` format.
    """

    if isinstance(file, str):
        handle = tables.open_file(file, "a")
        own_fid = True
    else:
        if not isinstance(file, tables.File):
            raise TypeError(
                "file must be either a string " "or an open tables.File: %s" % file,
            )
        handle = file
        own_fid = False

    # name all the arrays
    namedict = kwargs
    for i, val in enumerate(args):
        key = "arr_%d" % i
        if key in namedict.keys():
            if own_fid:
                handle.close()
            raise ValueError(
                "Cannot use un-named variables " " and keyword %s" % key,
            )
        namedict[key] = val

    # ensure that they don't already exist
    current_nodes = [e.name for e in handle.list_nodes(where="/")]

    for key in namedict.keys():
        if key in current_nodes:
            handle.remove_node("/", name=key)
            # per discussion on github, https://github.com/rmcgibbo/mdtraj/issues/5
            # silent overwriting appears to be the desired functionality
            # raise IOError('Array already exists in file: %s' % key)

    # save all the arrays
    try:
        for key, val in namedict.items():
            if not isinstance(val, np.ndarray):
                raise TypeError(
                    "Only numpy arrays can " f"be saved: type({key}) is {type(val)}",
                )
            try:
                atom = tables.Atom.from_dtype(val.dtype)
            except ValueError:
                raise TypeError(
                    "Arrays of this dtype " "cannot be saved: %s" % val.dtype,
                )

            node = handle.create_carray(
                where="/",
                name=key,
                atom=atom,
                shape=val.shape,
                filters=COMPRESSION,
            )

            node[:] = val

    except Exception:
        handle.close()
        if own_fid:
            os.unlink(file)
        raise

    handle.flush()
    if own_fid:
        handle.close()


def loadh(file, name=Ellipsis, deferred=True):
    """Load one or more array(s) from HDF5 format files

    Parameters
    ----------
    file : string or tables.File
        The file to read. It must be either a string, or a an open PyTables
        file handle
    name : string, optional
        The name of a single to read from the file. If not supplied, all arrays
        will be read
    deferred : bool, optional
        If true, and you did not request just a single name, the result will
        be lazyily loaded.

    Returns
    -------
    result : array or dict-like
        If name is a single string, a single array will be returned. Otherwise,
        the return value is a dict-like mapping the name(s) to the array(s) of
        data.

    Raises
    ------
    IOError
        If file does not exist
    KeyError
        If the request name does not exist

    See Also
    --------
    numpy.load : Load an array(s) or pickled objects from .npy, .npz, or pickled files.
    """

    if isinstance(file, str):
        handle = tables.open_file(file, mode="r")
        own_fid = True
    else:
        if not isinstance(file, tables.File):
            raise TypeError(
                "file must be either a string " "or an open tables.File: %s" % file,
            )
        handle = file
        own_fid = False

    # if name is a single string, deferred loading is not used
    if isinstance(name, str):
        try:
            node = handle.get_node(where="/", name=name)
        except tables.NoSuchNodeError:
            raise KeyError(
                f'Node "{name}" does not exist ' f"in file {file}",
            )

        return_value = np.array(node[:])
        if own_fid:
            handle.close()
        return return_value

    if not deferred:
        result = {}
        iterator = handle.walk_nodes(where="/")
        for node in iterator:
            if isinstance(node, tables.Array):
                # note that we want to strip off the leading "/"
                # also, we're skipping Tables and other hdf5 structures
                result[node._v_pathname[1:]] = node[:]
        if own_fid:
            handle.close()
        return result

    return DeferredTable(handle, own_fid)


class DeferredTable:
    def __init__(self, handle, own_fid):
        self._handle = handle

        # get the paths of all of the nodes that are arrays (note that)
        # we're skipping Tables
        self._node_names = [
            node._v_pathname[1:] for node in handle.walk_nodes(where="/") if isinstance(node, tables.Array)
        ]

        self._loaded = {}
        self._own_fid = own_fid

        repr_strings = []
        for name in self._node_names:
            repr_strings.append(
                "  {}: [shape={}, dtype={}]".format(
                    name,
                    handle.get_node(where="/", name=name).shape,
                    handle.get_node(where="/", name=name).dtype,
                ),
            )
        self._repr_string = "{\n%s\n}" % ",\n".join(repr_strings)

    def __repr__(self):
        return self._repr_string

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, "_own_fid") and self._own_fid:
            self._handle.close()

    def __getitem__(self, key):
        if key not in self._node_names:
            raise KeyError(f"{key} not in {self._node_names}")
        if key not in self._loaded:
            self._loaded[key] = self._handle.get_node(where="/", name=key)[:]
        return self._loaded[key]

    def iteritems(self):
        for name in self._node_names:
            yield (name, getattr(self, name))

    def keys(self):
        return self._node_names

    def iterkeys(self):
        return iter(self._node_names)

    def __contains__(self, key):
        return self._node_names.__contains__(key)
