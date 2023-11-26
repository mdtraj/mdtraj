##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors: Robert T. McGibbon
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

from mdtraj.utils import (ensure_type, cast_indices, in_units_of,
                          lengths_and_angles_to_box_vectors)
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils.six import string_types
from mdtraj.utils.six.moves import xrange

__all__ = ['LAMMPSTrajectoryFile', 'load_lammpstrj']


class _EOF(IOError):
    pass


@FormatRegistry.register_loader('.lammpstrj')
def load_lammpstrj(filename, top=None, stride=None, atom_indices=None,
                   frame=None, unit_set='real'):
    """Load a LAMMPS trajectory file.

    Parameters
    ----------
    filename : path-like
        Path of LAMMPS trajectory file.
    top : {str, Trajectory, Topology}
        The lammpstrj format does not contain topology information. Pass in
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
    unit_set : str, optional
        The LAMMPS unit set that the simulation was performed in. See
        http://lammps.sandia.gov/doc/units.html for options. Currently supported
        unit sets: 'real'.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.LAMMPSTrajectoryFile :  Low level interface to lammpstrj files
    """
    from mdtraj.core.trajectory import _parse_topology

    # We make `top` required. Although this is a little weird, its good because
    # this function is usually called by a dispatch from load(), where top comes
    # from **kwargs. So if its not supplied, we want to give the user an
    # informative error message.
    if top is None:
        raise ValueError('"top" argument is required for load_lammpstrj')
    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_lammpstrj. '
                        'you supplied %s'.format(type(filename)))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with LAMMPSTrajectoryFile(filename) as f:
        # TODO: Support other unit sets.
        if unit_set == 'real':
            f.distance_unit == 'angstroms'
        else:
            raise ValueError('Unsupported unit set specified: {0}.'.format(unit_set))
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, stride=stride, atom_indices=atom_indices)


@FormatRegistry.register_fileobject('.lammpstrj')
class LAMMPSTrajectoryFile(object):
    """Interface for reading and writing to a LAMMPS lammpstrj files.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    Parameters
    ----------
    filename : path-like
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
        """Open a LAMMPS lammpstrj file for reading/writing. """
        self._is_open = False
        self._filename = filename
        self._mode = mode
        self._frame_index = 0
        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0

        if mode == 'r':
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self._fh = open(filename, 'r')
            self._is_open = True
        elif mode == 'w':
            if not force_overwrite and os.path.exists(filename):
                raise IOError('"%s" already exists' % filename)
            self._fh = open(filename, 'w')
            self._is_open = True
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied "{0}"'.format(mode))

    def close(self):
        """Close the lammpstrj file. """
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

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from a lammpstrj file

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

        See Also
        --------
        read : Returns the raw data from the file

        Notes
        -----
        If coordinates are specified in more than one style, the first complete
        trio of x/y/z coordinates will be read in according to the following
        order:
            1) x,y,z (unscaled coordinates)
            2) xs,ys,zs (scaled atom coordinates)
            3) xu,yu,zu (unwrapped atom coordinates)
            4) xsu,ysu,zsu (scaled unwrapped atom coordinates)

        E.g., if the file contains x, y, z, xs, ys, zs then x, y, z will be used.
              if the file contains x, y, xs, ys, zs then xs, ys, zs will be used.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        initial = int(self._frame_index)
        xyz, cell_lengths, cell_angles = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(cell_lengths, self.distance_unit, Trajectory._distance_unit, inplace=True)

        if stride is None:
            stride = 1
        time = (stride*np.arange(len(xyz))) + initial

        t = Trajectory(xyz=xyz, topology=topology, time=time)
        t.unitcell_lengths = cell_lengths
        t.unitcell_angles = cell_angles
        return t

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a lammpstrj file.

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
        cell_lengths : np.ndarray, None
            The lengths (a,b,c) of the unit cell for each frame, or None if
            the information is not present in the file.
        cell_angles : np.ndarray, None
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame, or None if  the information is not present in the file.
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

        all_coords, all_lengths, all_angles = [], [], []
        for _ in frame_counter:
            try:
                frame_coords, frame_lengths, frame_angles = self._read()
                if atom_indices is not None:
                    frame_coords = frame_coords[atom_indices, :]
            except _EOF:
                break

            all_coords.append(frame_coords)
            all_lengths.append(frame_lengths)
            all_angles.append(frame_angles)

            for j in range(stride - 1):
                # throw away these frames
                try:
                    self._read()
                except _EOF:
                    break

        all_coords = np.array(all_coords)
        all_lengths = np.array(all_lengths, dtype=np.float32)
        all_angles = np.array(all_angles, dtype=np.float32)
        return all_coords, all_lengths, all_angles

    def parse_box(self, style):
        """Extract lengths and angles from a frame.

        Parameters
        ----------
        style : str
            Type of box, 'triclinic' or 'orthogonal'.

        Returns
        -------
            lengths : ndarray
            angles : ndarray

        Notes
        -----
        For more info on how LAMMPS defines boxes:
        http://lammps.sandia.gov/doc/Section_howto.html#howto_12
        """
        box = np.empty(shape=(3, 2))
        if style == 'triclinic':
            factors = np.empty(3)
            for i in range(3):
                line = self._fh.readline().split()
                box[i] = line[:2]
                factors[i] = line[2]
            xy, xz, yz = factors

            xlo = box[0, 0] - np.min([0.0, xy, xz, xy+xz])
            xhi = box[0, 1] - np.max([0.0, xy, xz, xy+xz])
            ylo = box[1, 0] - np.min([0.0, yz])
            yhi = box[1, 1] - np.max([0.0, yz])
            zlo = box[2, 0]
            zhi = box[2, 1]

            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo

            a = lx
            b = np.sqrt(ly**2 + xy**2)
            c = np.sqrt(lz**2 + xz**2 + yz**2)
            alpha = np.arccos((xy*xz + ly*yz) / (b*c))
            beta = np.arccos(xz / c)
            gamma = np.arccos(xy / b)

            lengths = np.array([a, b, c])
            angles = np.degrees(np.array([alpha, beta, gamma]))
        elif style == 'orthogonal':
            box[0] = self._fh.readline().split()  # x-dim of box
            box[1] = self._fh.readline().split()  # y-dim of box
            box[2] = self._fh.readline().split()  # z-dim of box
            lengths = np.diff(box, axis=1).reshape(1, 3)[0]  # box lengths
            angles = np.empty(3)
            angles.fill(90.0)
        return lengths, angles

    def _read(self):
        """Read a single frame. """

        # --- begin header ---
        first = self._fh.readline()  # ITEM: TIMESTEP
        if first == '':
            raise _EOF()
        self._fh.readline()  # timestep
        self._fh.readline()  # ITEM: NUMBER OF ATOMS
        self._n_atoms = int(self._fh.readline())  # num atoms

        box_header = self._fh.readline().split()  # ITEM: BOX BOUNDS
        self._line_counter += 5
        if len(box_header) == 9:
            lengths, angles = self.parse_box('triclinic')
        elif len(box_header) == 6:
            lengths, angles = self.parse_box('orthogonal')
        else:
            raise IOError('lammpstrj parse error on line {0:d} of "{1:s}". '
                          'This file does not appear to be a valid '
                          'lammpstrj file.'.format(
                    self._line_counter,  self._filename))

        column_headers = self._fh.readline().split()[2:]  # ITEM: ATOMS ...
        if self._frame_index == 0:
            # Detect which columns the atom index, type and coordinates are.
            columns = {header: idx for idx, header in enumerate(column_headers)}

            # Make sure the file contains an x, y, and z-coordinate of the same
            # style.
            coord_keywords = [('x', 'y', 'z'),  # unscaled
                              ('xs', 'ys', 'zs'),  # scaled
                              ('xu', 'yu', 'zu'),  # unwrapped
                              ('xsu', 'ysu', 'zsu')]  # scaled and unwrapped
            for keywords in coord_keywords:
                if set(keywords).issubset(column_headers):
                    break
            else:
                raise IOError('Invalid .lammpstrj file. Must contain x, y, and '
                              'z coordinates that all adhere to the same style.')

            try:
                self._atom_index_column = columns['id']
                self._atom_type_column = columns['type']
                self._xyz_columns = [columns[keywords[0]], columns[keywords[1]], columns[keywords[2]]]
            except KeyError:
                raise IOError("Invalid .lammpstrj file. Must contain 'id', "
                              "'type', 'x*', 'y*' and 'z*' entries.")
        self._line_counter += 4
        # --- end header ---

        xyz = np.empty(shape=(self._n_atoms, 3))
        types = np.empty(shape=self._n_atoms, dtype='int')

        # --- begin body ---
        for _ in xrange(self._n_atoms):
            line = self._fh.readline()
            if line == '':
                raise _EOF()
            split_line = line.split()
            try:
                atom_index = int(split_line[self._atom_index_column])
                types[atom_index - 1] = int(split_line[self._atom_type_column])
                xyz[atom_index - 1] = [float(split_line[column]) for column in self._xyz_columns]
            except Exception:
                raise IOError('lammpstrj parse error on line {0:d} of "{1:s}". '
                              'This file does not appear to be a valid '
                              'lammpstrj file.'.format(
                        self._line_counter,  self._filename))
            self._line_counter += 1
        # --- end body ---

        self._frame_index += 1
        return xyz, lengths, angles

    def write_box(self, lengths, angles, mins):
        """Write the box lines in the header of a frame.

        Parameters
        ----------
        lengths : np.ndarray, dtype=np.double, shape=(3, )
            The lengths (a,b,c) of the unit cell for each frame.
        angles : np.ndarray, dtype=np.double, shape=(3, )
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame.
        mins : np.ndarray, dtype=np.double, shape=(3, )
            The minimum coordinates in the x-, y- and z-directions.
        """
        if np.allclose(angles, np.array([90, 90, 90])):
            self._fh.write('ITEM: BOX BOUNDS pp pp pp\n')
            self._fh.write('{0} {1}\n'.format(mins[0], mins[0] + lengths[0]))
            self._fh.write('{0} {1}\n'.format(mins[1], mins[1] + lengths[1]))
            self._fh.write('{0} {1}\n'.format(mins[2], mins[2] + lengths[2]))
        else:
            a, b, c = lengths
            alpha, beta, gamma = np.radians(angles)

            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b*c*np.cos(alpha) - xy*xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

            xlo = mins[0]
            xhi = xlo + lx
            ylo = mins[1]
            yhi = ylo + ly
            zlo = mins[2]
            zhi = zlo + lz

            xlo_bound = xlo + np.min([0.0, xy, xz, xy+xz])
            xhi_bound = xhi + np.max([0.0, xy, xz, xy+xz])
            ylo_bound = ylo + np.min([0.0, yz])
            yhi_bound = yhi + np.max([0.0, yz])
            zlo_bound = zlo
            zhi_bound = zhi
            self._fh.write('ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
            self._fh.write('{0} {1} {2}\n'.format(xlo_bound, xhi_bound, xy))
            self._fh.write('{0} {1} {2}\n'.format(ylo_bound, yhi_bound, xz))
            self._fh.write('{0} {1} {2}\n'.format(zlo_bound, zhi_bound, yz))

    def write(self, xyz, cell_lengths, cell_angles=None, types=None, unit_set='real'):
        """Write one or more frames of data to a lammpstrj file.

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention,
            the lengths should be in units of angstroms.
        cell_lengths : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The lengths (a,b,c) of the unit cell for each frame. By convention,
            the lengths should be in units of angstroms.
        cell_angles : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame. (Units of degrees).
        types : np.ndarray, shape(3, ), dtype=int
            The numeric type of each particle.
        unit_set : str, optional
            The LAMMPS unit set that the simulation was performed in. See
            http://lammps.sandia.gov/doc/units.html for options. Currently supported
            unit sets: 'real'.
        """
        if not self._mode == 'w':
            raise ValueError('write() is only available when file is opened '
                             'in mode="w"')

        xyz = ensure_type(xyz, np.float32, 3, 'xyz', can_be_none=False,
                shape=(None, None, 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, np.float32, 2, 'cell_lengths',
                can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        if cell_angles is None:
            cell_angles = np.empty_like(cell_lengths)
            cell_angles.fill(90)
        cell_angles = ensure_type(cell_angles, np.float32, 2, 'cell_angles',
                can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        if not types:
            # Make all particles the same type.
            types = np.ones(shape=(xyz.shape[1]))
        types = ensure_type(types, int, 1, 'types', can_be_none=True,
                shape=(xyz.shape[1], ), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=False)

        # TODO: Support other unit sets.
        if unit_set == 'real':
            self.distance_unit == 'angstroms'
        else:
            raise ValueError('Unsupported unit set specified: {0}.'.format(unit_set))

        for i in range(xyz.shape[0]):
            # --- begin header ---
            self._fh.write('ITEM: TIMESTEP\n')
            self._fh.write('{0}\n'.format(i))  # TODO: Write actual time if known.
            self._fh.write('ITEM: NUMBER OF ATOMS\n')
            self._fh.write('{0}\n'.format(xyz.shape[1]))
            self.write_box(cell_lengths[i], cell_angles[i], xyz[i].min(axis=0))
            # --- end header ---

            # --- begin body ---
            self._fh.write('ITEM: ATOMS id type xu yu zu\n')
            for j, coord in enumerate(xyz[i]):
                self._fh.write('{0:d} {1:d} {2:8.3f} {3:8.3f} {4:8.3f}\n'.format(
                    j+1, types[j], coord[0], coord[1], coord[2]))
            # --- end body ---

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
        "Number of frames in the file"
        raise NotImplementedError()


