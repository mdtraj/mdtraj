##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Christoph Klein
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

from datetime import date
import itertools
import os

import numpy as np

from mdtraj.formats.registry import _FormatRegistry
from mdtraj.io import open_maybe_zipped
from mdtraj.utils import cast_indices, in_units_of, ensure_type
from mdtraj.utils.six import string_types
from mdtraj.utils.six.moves import xrange
from mdtraj.version import version

__all__ = ['XYZTrajectoryFile', 'load_xyz']


class _EOF(IOError):
    pass


@_FormatRegistry.register_loader('.xyz')
@_FormatRegistry.register_loader('.xyz.gz')
def load_xyz(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load a xyz trajectory file.

    While there is no universal standard for this format, this plugin adheres
    to the same format as the VMD plugin:

    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

    Most notably, units are in angstroms and anything past the 'z' field is
    ignored.

    Parameters
    ----------
    filename : str
        String filename of xyz trajectory file.
    top : {str, Trajectory, Topology}
        The xyz format does not contain topology information. Pass in
        either the path to a pdb file, a trajectory, or a topology to supply
        this information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file.
    frame : int, optional
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.XYZTrajectoryFile :  Low level interface to xyz files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    # We make `top` required. Although this is a little weird, its good because
    # this function is usually called by a dispatch from load(), where top comes
    # from **kwargs. So if its not supplied, we want to give the user an
    # informative error message.
    if top is None:
        raise ValueError('"top" argument is required for load_xyz')

    if not isinstance(filename, string_types):
        raise TypeError('filename must be of type string for load_xyz. '
                        'you supplied %s'.format(type(filename)))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)
    if atom_indices is not None:
        topology = topology.subset(atom_indices)

    with XYZTrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            xyz = f.read(n_frames=1, atom_indices=atom_indices)
        else:
            xyz = f.read(stride=stride, atom_indices=atom_indices)
        in_units_of(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)

    time = np.arange(len(xyz))
    if frame is not None:
        time += frame
    elif stride is not None:
        time *= stride

    t = Trajectory(xyz=xyz, topology=topology, time=time)
    return t


@_FormatRegistry.register_fileobject('.xyz')
@_FormatRegistry.register_fileobject('.xyz.gz')
class XYZTrajectoryFile(object):
    """Interface for reading and writing to xyz files.

    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for
        write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?
    """

    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=True):
        """Open a xyz file for reading/writing. """
        self._is_open = False
        self._filename = filename
        self._mode = mode
        self._frame_index = 0
        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0

        if mode == 'r':
            self._fh = open_maybe_zipped(filename, 'r')
            self._is_open = True
        elif mode == 'w':
            self._fh = open_maybe_zipped(filename, 'w', force_overwrite)
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "{0}"'.format(mode))

    def close(self):
        """Close the xyz file. """
        if self._is_open:
            self._fh.close()
            self._is_open = False

    def __del__(self):
        self.close()

    def __enter__(self):
        """Support the context manager protocol. """
        return self

    def __exit__(self, *exc_info):
        """Support the context manager protocol. """
        self.close()

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a xyz file.

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates
            from the file.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
        """
        if not self._mode == 'r':
            raise ValueError('read() is only available when file is opened '
                             'in mode="r"')

        if n_frames is None:
            frame_counter = itertools.count()
        else:
            frame_counter = xrange(n_frames)

        if stride is None:
            stride = 1

        all_coords = []
        for i in frame_counter:
            try:
                frame_coords = self._read()
                if atom_indices is not None:
                    frame_coords = frame_coords[atom_indices, :]
            except _EOF:
                break

            all_coords.append(frame_coords)

            for j in range(stride - 1):
                # throw away these frames
                try:
                    self._read()
                except _EOF:
                    break

        all_coords = np.array(all_coords)
        return all_coords

    def _read(self):
        """Read a single frame. """

        first = self._fh.readline()  # Number of atoms.
        if first == '':
            raise _EOF()
        else:
            self._n_atoms = int(first)
        self._fh.readline()  # Comment line.
        self._line_counter += 2

        xyz = np.empty(shape=(self._n_atoms, 3))
        types = np.empty(shape=self._n_atoms, dtype=str)

        for i in xrange(self._n_atoms):
            line = self._fh.readline()
            if line == '':
                raise _EOF()
            split_line = line.split()
            try:
                types[i] = split_line[0]
                xyz[i] = [float(x) for x in split_line[1:4]]
            except Exception:
                raise IOError('xyz parse error on line {0:d} of "{1:s}". '
                              'This file does not appear to be a valid '
                              'xyz file.'.format(
                        self._line_counter,  self._filename))
            self._line_counter += 1
        # --- end body ---

        self._frame_index += 1
        return xyz

    def write(self, xyz, types=None):
        """Write one or more frames of data to a xyz file.

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write.
        types : np.ndarray, shape(3, )
            The type of each particle.
        """

        if not self._mode == 'w':
            raise ValueError('write() is only available when file is opened '
                             'in mode="w"')

        if not types:
            # Make all particles the same type.
            types = ['X' for _ in xrange(xyz.shape[1])]
        xyz = ensure_type(xyz, np.float32, 3, 'xyz', can_be_none=False,
                        shape=(None, None, 3), warn_on_cast=False,
                        add_newaxis_on_deficient_ndim=True)
        in_units_of(xyz, 'nanometers', self.distance_unit, inplace=True)

        for i in range(xyz.shape[0]):
            self._fh.write('{0}\n'.format(xyz.shape[1]))
            self._fh.write("Created with MDTraj {0}, {1}\n".format(version, str(date.today())))

            for j, coord in enumerate(xyz[i]):
                self._fh.write('{0} {1:8.3f} {2:8.3f} {3:8.3f}\n'.format(
                    types[j], coord[0], coord[1], coord[2]))

    def seek(self, offset, whence=0):
        """Move to a new file position.

        Parameters
        ----------
        offset : int
            A number of frames.
        whence : {0, 1, 2}
            0: offset from start of file, offset should be >=0.
            1: move relative to the current position, positive or negative
            2: move relative to the end of file, offset should be <= 0.
            Seeking beyond the end of a file is not supported
        """
        if self._mode == 'r':
            advance, absolute = None, None
            if whence == 0 and offset >= 0:
                if offset >= self._frame_index:
                    advance = offset - self._frame_index
                else:
                    absolute = offset
            elif whence == 1 and offset >= 0:
                advance = offset
            elif whence == 1 and offset < 0:
                absolute = offset + self._frame_index
            elif whence == 2 and offset <= 0:
                raise NotImplementedError('offsets from the end are not supported yet')
            else:
                raise IOError('Invalid argument')

            if advance is not None:
                for i in range(advance):
                    self._read()  # advance and throw away these frames
            elif absolute is not None:
                self._fh.close()
                self._fh = open(self._filename, 'r')
                self._frame_index = 0
                self._line_counter = 0
                for i in range(absolute):
                    self._read()
            else:
                raise RuntimeError()

        else:
            raise NotImplementedError('offsets in write mode are not supported yet')

    def tell(self):
        """Current file position.

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return int(self._frame_index)

    def __len__(self):
        """Number of frames in the file. """
        raise NotImplementedError()
