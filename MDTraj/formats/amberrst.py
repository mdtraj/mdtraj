##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Jason Swails
# Contributors:
#
# This code for reading Amber restart and inpcrd files was taken from ParmEd,
# which is released under the GNU Lesser General Public License
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

"""
This module provides the ability to read Amber inpcrd/restart files as well as
Amber NetCDF restart files. This code was taken from ParmEd and simplified by
removing the functionality that is not needed.
"""

from __future__ import print_function, division
from math import ceil
import os

import numpy as np
from mdtraj import version
from mdtraj.formats.registry import _FormatRegistry
from mdtraj.utils import ensure_type, import_, in_units_of, cast_indices

__all__ = ['AmberRestartFile', 'load_restrt', 'AmberNetCDFRestartFile',
           'load_ncrestrt']

try:
    xrange
except NameError:
    xrange = range

@_FormatRegistry.register_loader('.rst7')
@_FormatRegistry.register_loader('.restrt')
@_FormatRegistry.register_loader('.inpcrd')
def load_restrt(filename, top=None, atom_indices=None):
    """Load an AMBER ASCII restart/inpcrd file. Since this file doesn't contain
    information to specify the topology, you need to supply a topology

    Parameters
    ----------
    filename : str
        name of the AMBER restart file
    top : {str, Trajectory, Topology}
        Pass in either the path to a file containing topology information (e.g.,
        a PDB, an AMBER prmtop, or certain types of Trajectory objects) to
        supply the necessary topology information that is not present in these
        files
    atom_indices : array_like, optional
        If not None, then read only a subset of the atoms coordinates from the
        file.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object

    See Also
    --------
    mdtraj.AmberRestartFile : Low level interface to AMBER restart files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)
    if atom_indices is not None:
        topology = topology.subset(atom_indices)

    with AmberRestartFile(filename) as f:
        xyz, time, cell_lengths, cell_angles = f.read(atom_indices=atom_indices)
        xyz = in_units_of(xyz, f.distance_unit, Trajectory._distance_unit,
                          inplace=True)
        cell_lengths = in_units_of(cell_lengths, f.distance_unit,
                                   Trajectory._distance_unit, inplace=True)

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time,
                            unitcell_lengths=cell_lengths,
                            unitcell_angles=cell_angles)
    return trajectory

@_FormatRegistry.register_fileobject('.rst7')
@_FormatRegistry.register_fileobject('.restrt')
@_FormatRegistry.register_fileobject('.inpcrd')
class AmberRestartFile(object):
    """Interface for reading and writing AMBER NetCDF files. This is a file-like
    object, that supports both reading and writing depending on the `mode` flag.
    It implements the context manager protocol, so you can also use it with the
    python 'with' statement.

    Parameters
    ----------
    filename : str
        The name of the file to open
    mode : {'r', 'w'}, default='r'
        The mode in which to open the file. Valid options are 'r' or 'w' for
        'read' or 'write'
    force_overwrite : bool, default=False
        In write mode, if a file named `filename` already exists, clobber it and
        overwrite it
    """
    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=True):
        self._closed = True
        self._mode = mode
        self._filename = filename

        if mode not in ('r', 'w'):
            raise ValueError("mode must be one of ['r', 'w']")

        if mode == 'w' and not force_overwrite and os.path.exists(filename):
            raise IOError('"%s" already exists' % filename)

        if mode == 'w':
            self._needs_initialization = True
            self._handle = open(filename, mode)
            self._closed = False
        elif mode == 'r':
            with open(filename, mode) as f:
                f.readline()
                words = f.readline().split()
                try:
                    self._n_atoms = int(words[0])
                except (IndexError, ValueError):
                    raise TypeError('"%s" is not a recognized Amber restart' %
                                    filename)
            self._needs_initialization = False
        else:
            raise RuntimeError()

    @property
    def n_atoms(self):
        self._validate_open()
        if self._needs_initialization:
            raise IOError('The file is uninitialized')
        return self._n_atoms

    def _validate_open(self):
        if self._closed:
            raise IOError('The file is closed.')

    def _parse(self, lines):
        """ Parses the file """
        self._time = None
        try:
            words = lines[1].split()
            self._n_atoms = natom = int(words[0])
        except (IndexError, ValueError):
            raise TypeError('not a recognized Amber restart')

        time = None
        if len(words) > 2:
            time = float(words[1])

        lines_per_frame = int(ceil(natom / 2))
        if len(lines) == lines_per_frame + 2:
            hasbox = hasvels = False
        elif natom in (1, 2) and len(lines) == 4:
            # This is the _only_ case where line counting does not work -- there
            # is either 1 or 2 atoms and there are 4 lines. The 1st 3 lines are
            # the title, natom/time, and coordinates. The 4th are almost always
            # velocities since it's hard to have a periodic system this small.
            # However, velocities (which are scaled down by 20.445) have a ~0%
            # chance of being 60+, so we can pretty easily tell if the last line
            # has box dimensions and angles or velocities. I cannot envision a
            # plausible scenario where the detection here will ever fail
            line = lines[3]
            if natom == 1:
                tmp = [line[i:i+12] for i in range(0, 72, 12) if
                        line[i:i+12].strip()]
                if len(tmp) == 3:
                    hasvels = True
                    hasbox = False
                elif len(tmp) == 6:
                    hasbox = True
                    hasvels = False
                else:
                    raise TypeError('not a recognized Amber restart')
            else:
                # Ambiguous case
                tmp = [float(line[i:i+12]) >= 60.0 for i in range(0, 72, 12)]
                if any(tmp):
                    hasbox = True
                    hasvels = False
                else:
                    hasvels = True
                    hasbox = False
        elif len(lines) == lines_per_frame + 3:
            hasbox = True
            hasvels = False
        elif len(lines) == 2*lines_per_frame + 2:
            hasbox = False
            hasvels = True
        elif len(lines) == 2*lines_per_frame + 3:
            hasbox = hasvels = True
        else:
            raise TypeError('Badly formatted restart file. Has %d lines for '
                            '%d atoms' % (len(lines), natom))

        coordinates = np.zeros((1, 3, natom))
        if time is None:
            time = np.zeros(1)
        else:
            time = np.asarray((time,))

        # Fill the coordinates
        for i in xrange(lines_per_frame):
            line = lines[i]
            i2 = i * 2
            coordinates[0,i2,:] = [float(line[j:j+12]) for j in xrange(0,36,12)]
            i2 += 1
            if i2 < natom:
                coordinates[0,i2,:] = [float(line[j:j+12]) for j in
                                       xrange(36,72,12)]
        if hasbox:
            cell_lengths = np.zeros((1,3))
            cell_angles = np.zeros((1,3))
            line = lines[-1]
            cell_lengths[0,:] = [float(line[i:i+12]) for i in xrange(0,36,12)]
            cell_angles[0,:] = [float(line[i:i+12]) for i in xrange(36,72,12)]
        else:
            cell_lengths = cell_angles = None

        return coordinates, time, cell_lengths, cell_angles

    def read(self, atom_indices=None):
        """Read data from an AMBER ASCII restart file

        Parameters
        ----------
        atom_indices : np.ndarray, dtype=int, optional
            The specific indices of the atoms you'd like to retrieve. If not
            supplied, all of the atoms will be retrieved.

        Returns
        -------
        coordinates : np.ndarray, shape=(1, n_atoms, 3)
            The cartesian coordinates of the atoms, in units of angstroms. These
            files only ever contain 1 frame
        time : np.ndarray, None
            The time corresponding to the frame, in units of picoseconds, or
            None if no time information is present
        cell_lengths : np.ndarray, None
            The lengths (a, b, c) of the unit cell for the frame in angstroms,
            or None if the information is not present in the file
        cell_angles : np.ndarray, None
            The angles (\alpha, \beta, \gamma) defining the unit cell for each
            frame, or None if the information is not present in the file.
        """
        if self._mode != 'r':
            raise IOError('The file was opened in mode=%s. Reading is not '
                          'allowed.' % self._mode)

        with open(self._filename, 'r') as f:
            lines = f.readlines()
            
        coordinates, time, cell_lengths, cell_angles = self._parse(lines)

        if atom_indices is not None:
            atom_slice = ensure_type(atom_indices, dtype=np.int, ndim=1,
                                     name='atom_indices', warn_on_cast=False)
            if not np.all(atom_slice) >= 0:
                raise ValueError('Entries in atom_slice must be >= 0')
            coordinates = coordinates[:, atom_slice, :]

        return coordinates, time, cell_lengths, cell_angles

    def write(self, coordinates, time=None, cell_lengths=None,
              cell_angles=None):
        """Write one frame of a MD trajectory to disk in the AMBER ASCII restart
        file format.

        Parameters
        ----------
        coordinates : np.ndarray, dtype=np.float32, shape=([1,] n_atoms, 3)
            The cartesian coordinates of each atom, in units of angstroms. Must
            be only a single frame (shape can be (1,N,3) or (N,3) where N is
            the number of atoms)
        time : array-like with 1 element or float, optional
            The time corresponding to this frame. If not specified, a place
            holder of 0 will be written
        cell_lengths : np.ndarray, dtype=np.double, shape=([1,] 3)
            The lengths (a,b,c) of the unit cell for the frame in Angstroms
        cell_angles : np.ndarray, dtype=np.double, shape=([1,] 3)
            The angles between the unit cell vectors for the frame in Degrees
        """
        if self._mode != 'w':
            raise IOError('The file was opened in mode=%s. Writing not allowed.'
                          % self._mode)
        if not self._needs_initialization:
            # Must have already been written -- can only write once
            raise RuntimeError('restart file has already been written -- can '
                               'only write one frame to restart files.')
        coordinates = in_units_of(coordinates, None, 'angstroms')
        time = in_units_of(time, None, 'picoseconds')
        cell_lengths = in_units_of(cell_lengths, None, 'angstroms')
        cell_angles = in_units_of(cell_angles, None, 'degrees')

        # typecheck all of the input arguments rigorously
        coordinates = ensure_type(coordinates, np.float32, 3, 'coordinates',
                                  length=None, can_be_none=False,
                                  shape=(1,None,3), warn_on_cast=False,
                                  add_newaxis_on_deficient_ndim=True)
        n_frames, self._n_atoms = coordinates.shape[0], coordinates.shape[1]
        if n_frames != 1:
            raise ValueError('Can only write 1 frame to a restart file!')
        if time is not None:
            try:
                time = float(time)
            except TypeError:
                raise TypeError('Can only provide a single time')
        else:
            time = 0.0
        cell_lengths = ensure_type(cell_lengths, np.float64, 2, 'cell_lengths',
                                   length=1, can_be_none=True,
                                   warn_on_cast=False,
                                   add_newaxis_on_deficient_ndim=True)
        cell_angles = ensure_type(cell_lengths, np.float64, 2, 'cell_angles',
                                  length=1, can_be_none=True,
                                  warn_on_cast=False,
                                  add_newaxis_on_deficient_ndim=True)
        if ((cell_lengths is None and cell_angles is not None) or
            (cell_lengths is not None and cell_angles is None)):
            prov, negl = 'cell_lengths', 'cell_angles'
            if cell_lengths is None:
                prov, negl = negl, prov
            raise ValueError('You provided the variable "%s" but did not '
                             'provide "%s". Either provide both or neither -- '
                             'one without the other is meaningless.' %
                             (prov, negl))

        self._handle.write('Amber restart file (without velocities) written by '
                           'MDTraj\n')
        self._handle.write('%5d%15.7e\n' % (self._n_atoms, time))
        fmt = '%12.7f%12.7f%12.7f'
        for i in xrange(self._n_atoms):
            acor = coordinates[0, i, :]
            self._handle.write(fmt % (acor[0], acor[1], acor[2]))
            if i % 2 == 1: self._handle.write('\n')
        if self._n_atoms % 2 == 1: self._handle.write('\n')
        if cell_lengths is not None:
            self._handle.write(fmt % (cell_lengths[0,0], cell_lengths[0,1],
                                      cell_lengths[0,2]))
            self._handle.write(fmt % (cell_angles[0,0], cell_angles[0,1],
                                      cell_angles[0,2]))
        self._handle.flush()
