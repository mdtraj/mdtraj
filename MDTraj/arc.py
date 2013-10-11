# Copyright 2013-present mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import os
import itertools
import numpy as np
from mdtraj.utils import ensure_type

##############################################################################
# Classes
##############################################################################

class _EOF(IOError):
    pass


class ArcTrajectoryFile(object):
    """Interface for reading and writing to an TINKER archive files.
    (Note that the TINKER .xyz format is identical to this.)  This is
    a file-like object, that both reading or writing depending on the
    `mode` flag. It implements the context manager protocol, so you
    can also use it with the python 'with' statement.

    The conventional units in the arc file is angstrom. The format only
    supports storing the cartesian coordinates and box lengths.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r'}
        The mode in which to open the file, only 'r' for read is supported.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?
    """

    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=True):
        """Open an TINKER.arc file for reading/writing.
        """
        self._is_open = False
        self._filename = filename
        self._mode = mode

        if mode == 'w':
            raise ValueError('Writing TINKER .arc files is not supported at this time')

        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0

        if mode == 'r':
            # if n_atoms is None:
            #     raise ValueError('To open a mdcrd file in mode="r", you must '
            #                      'supply the number of atoms, "n_atoms"')
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self._fh = open(filename, 'r')
            self._is_open = True
        else:
            raise ValueError('mode must be "r". '
                             'you supplied "%s"' % mode)

    def close(self):
        """Close the .arc file"""
        if self._is_open:
            self._fh.close()
            self._is_open = False

    def __del__(self):
        self.close()

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a TINKER .arc file.  Note that only the
        Cartesian coordinates are read in.  The .arc file also
        contains TINKER-specific numeric atom types and some bonding
        information, which we do not read in.

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
            The cartesian coordinates, in angstroms
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

        coords = []
        for i in frame_counter:
            try:
                coord = self._read()
                if atom_indices is not None:
                    coord = coord[atom_indices, :]
            except _EOF:
                break

            coords.append(coord)

            for j in range(stride - 1):
                # throw away these frames
                self._read()

        coords = np.array(coords)

        return coords

    def _read(self):
        "Read a single frame"
        i = 0

        # Read in the number of atoms.
        line = self._fh.readline()
        if line == '':
            raise _EOF()
        
        self._n_atoms = int(line.split()[0])
        self._line_counter += 1

        coords = np.empty((self._n_atoms, 3), dtype=np.float32)
        while i < self._n_atoms:
            line = self._fh.readline()
            s = line.split()
            coords[i,:] = [float(s[pos]) for pos in [2, 3, 4]]
            i += 1
            self._line_counter += 1
        return coords

    def write(self, xyz):
        """ The ArcTrajectoryFile does not have a write method,
        because TINKER .arc files have special numerical atom types
        which are not shared by any other trajectory file format.

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. 

        """
        
        raise RuntimeError('write() is not available for .arc files')

