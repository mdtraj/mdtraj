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
import os
import itertools
import numpy as np
from mdtraj.utils import ensure_type, cast_indices, in_units_of
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils.six import string_types, PY3
from mdtraj.utils.six.moves import xrange

__all__ = ['MDCRDTrajectoryFile', 'load_mdcrd']


##############################################################################
# Classes
##############################################################################

class _EOF(IOError):
    pass

@FormatRegistry.register_loader('.mdcrd')
@FormatRegistry.register_loader('.crd')
def load_mdcrd(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load an AMBER mdcrd file.

    Parameters
    ----------
    filename : path-like
        Path of AMBER mdcrd file.
    top : {str, Trajectory, Topology}
        The BINPOS format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
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
    mdtraj.MDCRDTrajectoryFile :  Low level interface to MDCRD files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_mdcrd')

    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_mdcrd. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with MDCRDTrajectoryFile(filename, topology.n_atoms) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None
        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)


@FormatRegistry.register_fileobject('.mdcrd')
@FormatRegistry.register_fileobject('.crd')
class MDCRDTrajectoryFile(object):
    """Interface for reading and writing to an AMBER mdcrd files.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The conventional units in the mdcrd file are angstroms. The format only
    supports storing the cartesian coordinates and box lengths.

    Parameters
    ----------
    filename : path-like
        The filename to open. A path to a file on disk.
    n_atoms : int
        The number of atoms in the system. This is _required_ when mode == 'r'
        and irrelevant when mode == 'w'.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for
        write.
    has_box = 'detect'
        Does the mdcrd file contain box length information? This is optional
        when mode == 'r' (and irrelevant when mode == 'w'). The presence or
        absence of box information can generally be inferred from the file,
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
        self._w_has_box = None
        self._frame_index = 0
        self._has_box = has_box
        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0

        if has_box not in [True, False, "detect"]:
            raise ValueError('has_box must be one of [True, False, "detect"]')

        if mode == 'r':
            if n_atoms is None:
                raise ValueError('To open a mdcrd file in mode="r", you must '
                                 'supply the number of atoms, "n_atoms"')
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self._fh = open(filename, 'rb')
            self._is_open = True
            self._fh.readline()  # read comment
            self._line_counter += 1
        elif mode == 'w':
            if os.path.exists(filename) and not force_overwrite:
                raise IOError('"%s" already exists' % filename)
            self._fh = open(filename, 'wb')
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

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from a mdcrd file

        Parameters
        ----------
        topology : Topology
            The system topology
        n_frames : int, optional
            If positive, then read only the next `n_frames` frames. Otherwise read all
            of the frames in the file.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        trajectory : Trajectory
            A trajectory object containing the loaded portion of the file.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        initial = int(self._frame_index)
        xyz, cell_lengths = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(cell_lengths, self.distance_unit, Trajectory._distance_unit, inplace=True)

        if cell_lengths is None:
            cell_angles = None
        else:
            # Assume that its a rectilinear box
            cell_angles = 90.0 * np.ones_like(cell_lengths)

        if stride is None:
            stride = 1
        time = (stride*np.arange(len(xyz))) + initial

        t = Trajectory(xyz=xyz, topology=topology, time=time)
        t.unitcell_lengths = cell_lengths
        t.unitcell_angles = cell_angles
        return t

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
        cell_lengths : {np.ndarray, None}
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
                try:
                    self._read()
                except _EOF:
                    break

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
            self._line_counter += 1

            if line == b'':
                raise _EOF()
            try:
                items = [float(line[j:j+8])
                         for j in range(0, len(line.rstrip()), 8)]
                assert 0 < len(items) <= 10
            except Exception:
                raise IOError('mdcrd parse error on line %d of "%s". This file '
                              'does not appear to be a valid mdcrd file.' % \
                              (self._line_counter,  self._filename))

            length = len(items)

            if i + length > len(coords):
                raise IOError(
                    'mdcrd parse error: specified n_atoms (%d) is likely incorrect. '
                    'Incorrect buffer size encountered on line=%d' % (
                    self._n_atoms, self._line_counter))

            coords[i:i+length] = items
            i += length

            if i == self._n_atoms * 3:
                if self._has_box is False:
                    break

                # peek ahead for box
                here = self._fh.tell()
                line = self._fh.readline()
                peek = [float(elem) for elem in line.strip().split()]
                if len(peek) == 3:
                    box = peek
                else:
                    if self._has_box is True:
                        raise IOError('Box information not found in file.')
                    self._fh.seek(-len(line), 1)
                    self._fh.seek(here)
                break

        self._frame_index += 1
        return coords.reshape(self._n_atoms, 3), box

    def write(self, xyz, cell_lengths=None):
        """Write one or more frames of data to a mdcrd file

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention, the
            lengths should be in units of angstroms.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
            The length of the periodic box in each frame, in each direction,
            `a`, `b`, `c`. By convention the lengths should be in units
            of angstroms.
        """
        if not self._mode == 'w':
            raise ValueError('write() is only available when file is opened '
                             'in mode="w"')

        xyz = ensure_type(xyz, np.float32, 3, 'xyz', can_be_none=False,
                shape=(None, None, 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, np.float32, 2, 'cell_lengths',
                can_be_none=True, shape=(len(xyz), 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)

        if self._w_has_box is None:
            # this is the first write()
            self._n_atoms = xyz.shape[1]
            comment = 'TITLE : Created by MDTraj with %d atoms\n' % self._n_atoms
            if PY3:
                comment = comment.encode('ascii')
            self._fh.write(comment)

            if cell_lengths is None:
                self._w_has_box = False
            else:
                self._w_has_box = True
        elif self._w_has_box is True:
            if cell_lengths is None:
                raise ValueError('This mdcrd file must contain unitcell '
                                 'information')
        elif self._w_has_box is False:
            if cell_lengths is not None:
                raise ValueError('This mdcrd file must not contain unitcell '
                                 'information')
        else:
            raise RuntimeError()

        for i in range(xyz.shape[0]):
            for j, coord in enumerate(xyz[i].reshape(-1)):
                lfdone = False
                out = "%8.3f" % coord
                if len(out) > 8:
                    raise ValueError('Overflow error')
                if PY3:
                    out = out.encode('ascii')
                self._fh.write(out)
                if (j+1) % 10 == 0:
                    self._fh.write(b"\n")
                    lfdone = True

            if not lfdone:
                self._fh.write(b"\n")

            if cell_lengths is not None:
                line = "%8.3f %8.3f %8.3f\n" % tuple(cell_lengths[i])
                if PY3:
                    line = line.encode('ascii')
                self._fh.write(line)

    def seek(self, offset, whence=0):
        """Move to a new file position

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
                self._fh = open(self._filename, 'rb')
                self._fh.readline()  # read comment
                self._frame_index = 0
                self._line_counter = 1
                for i in range(absolute):
                    self._read()
            else:
                raise RuntimeError()

        else:
            raise NotImplementedError('offsets in write mode are not supported yet')

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return int(self._frame_index)

    def __len__(self):
        "Number of frames in the file"
        raise NotImplementedError()
