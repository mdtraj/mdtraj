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

import os
import itertools
import numpy as np


##############################################################################
# Classes
##############################################################################

class _EOF(IOError):
    pass


class MDCRDTrajectoryFile(object):
    """Interface for reading and writing to an AMBER mdcrd files.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The conventional units in the mdcrd file are angstroms. The format only
    supports storing the cartesian coordinates and box lengths.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    n_atoms : int
        The number of atoms in the system. This is _required_ when mode == 'r'
        and irrelivant when mode == 'w'.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for
        write.
    has_box = 'detect'
        Does the mdcrd file contain box length information? This is optional
        when mode == 'r' (and irrelvant when mode == 'w'). The presence or
        absense of box information can generally be infered from the file,
        but there might be corner cases in which this is not possible,
        because of limitations in the mdcrd format.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?
    """

    distance_unit = 'angstroms'

    def __init__(self, filename,  n_atoms=None, mode='r', has_box='detect',
                 force_overwrite=True):
        """Open an AMBER mdcrd file for reading/writing.
        """
        self._is_open = False
        self._filename = filename
        self._n_atoms = n_atoms
        self._mode = mode
        self._has_box = has_box

        if has_box not in [True, False, "detect"]:
            raise ValueError('has_box must be one of [True, False, "detect"]')

        if mode == 'r':
            if n_atoms is None:
                raise ValueError('To open a mdcrd file in mode="r", you must '
                                 'supply the number of atoms, "n_atoms"')
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self._fh = open(filename, 'r')
            self._is_open = True
            self._fh.readline()  # read comment
        elif mode == 'w':
            if os.path.exists(filename) and not force_overwrite:
                raise IOError("The file '%s' already exists" % filename)
            self._fh = open(filename, 'w')
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "%s"' % mode)

    def close(self):
        """Close the mdcrd file"""
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
        """Read data from a mdcrd file

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
        unitcell_lengths : {np.ndarray, None}
            If the file contains unitcell lengths, they will be returned as an
            array of shape=(n_frames, 3). Otherwise, unitcell_angles will be
            None.
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

        coords, boxes = [], []
        for i in frame_counter:
            try:
                coord, box = self._read()
                if atom_indices is not None:
                    coord = coord[atom_indices, :]
            except _EOF:
                break

            coords.append(coord)
            boxes.append(box)

            for j in range(stride - 1):
                # throw away these frames
                self._read()

        coords = np.array(coords)
        if all(b is None for b in boxes):
            # if there was no box information in any frame, that's cool
            return coords, None

        if not all(b is not None for b in boxes):
            # but if some of them had box information and others didn't
            # that probably means there was a bug in the parsing.
            raise IOError('Inconsistent box information. Try manually '
                          'setting has_box? Your mdcrd file might be '
                          'corrupt.')

        return coords, np.array(boxes, dtype=np.float32)

    def _read(self):
        "Read a single frame"
        i = 0
        coords = np.empty(self._n_atoms*3, dtype=np.float32)
        box = None

        while i < self._n_atoms * 3:
            line = self._fh.readline()
            if line == '':
                raise _EOF()

            items = [float(elem) for elem in line.split()]
            length = len(items)

            if i + length > len(coords):
                raise IOError('Incorrct buffer shape. The number of atoms '
                              'you specified is likely incorrect.')

            coords[i:i+length] = items
            i += length

            if i == self._n_atoms * 3:
                if self._has_box is False:
                    break

                # peek ahead for box
                here = self._fh.tell()
                line = self._fh.readline()
                peek = [float(elem) for elem in line.split()]
                if len(peek) == 3:
                    box = peek
                else:
                    if self._has_box is True:
                        raise IOError('Box information not found in file.')
                    self._fh.seek(here)
                break

        return coords.reshape(self._n_atoms, 3), box


f = MDCRDTrajectoryFile('/Users/rmcgibbo/local/mdtraj/frame0.mdcrd', n_atoms=22)
print f.read()[0][-1]

f = MDCRDTrajectoryFile('/Users/rmcgibbo/local/mdtraj/frame0.mdcrdbox', n_atoms=22)
print f.read()[0][-1]
