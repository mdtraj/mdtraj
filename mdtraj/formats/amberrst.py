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
from distutils.version import StrictVersion
from math import ceil
import os
import warnings

import numpy as np
from mdtraj import version
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import ensure_type, import_, in_units_of, cast_indices, six

__all__ = ['AmberRestartFile', 'load_restrt', 'AmberNetCDFRestartFile',
           'load_ncrestrt']

range = six.moves.range

@FormatRegistry.register_loader('.rst7')
@FormatRegistry.register_loader('.restrt')
@FormatRegistry.register_loader('.inpcrd')
def load_restrt(filename, top=None, atom_indices=None):
    """Load an AMBER ASCII restart/inpcrd file. Since this file doesn't contain
    information to specify the topology, you need to supply a topology

    Parameters
    ----------
    filename : path-like
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
    from mdtraj.core.trajectory import _parse_topology

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with AmberRestartFile(filename) as f:
        return f.read_as_traj(topology, atom_indices=atom_indices)


@FormatRegistry.register_fileobject('.rst7')
@FormatRegistry.register_fileobject('.restrt')
@FormatRegistry.register_fileobject('.inpcrd')
class AmberRestartFile(object):
    """Interface for reading and writing AMBER ASCII restart files. This is a
    file-like object, that supports both reading and writing depending on the
    `mode` flag.  It implements the context manager protocol, so you can also
    use it with the python 'with' statement.

    Parameters
    ----------
    filename : path-like
        The name of the file to open
    mode : {'r', 'w'}, default='r'
        The mode in which to open the file. Valid options are 'r' or 'w' for
        'read' or 'write'
    force_overwrite : bool, default=False
        In write mode, if a file named `filename` already exists, clobber it and
        overwrite it

    See Also
    --------
    md.AmberNetCDFRestartFile : Low level interface to AMBER NetCDF-format restart files
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

    @property
    def n_frames(self):
        return 1 # always 1 frame

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
        if len(words) >= 2:
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

        coordinates = np.zeros((1, natom, 3))
        if time is None:
            time = np.zeros(1)
        else:
            time = np.asarray((time,))

        # Fill the coordinates
        for i in range(lines_per_frame):
            line = lines[i+2]  # Skip first two lines
            i2 = i * 2
            coordinates[0,i2,:] = [float(line[j:j+12]) for j in range(0,36,12)]
            i2 += 1
            if i2 < natom:
                coordinates[0,i2,:] = [float(line[j:j+12]) for j in
                                       range(36,72,12)]
        if hasbox:
            cell_lengths = np.zeros((1,3))
            cell_angles = np.zeros((1,3))
            line = lines[-1]
            cell_lengths[0,:] = [float(line[i:i+12]) for i in range(0,36,12)]
            cell_angles[0,:] = [float(line[i:i+12]) for i in range(36,72,12)]
        else:
            cell_lengths = cell_angles = None

        return coordinates, time, cell_lengths, cell_angles

    def read_as_traj(self, topology, atom_indices=None):
        """Read an AMBER ASCII restart file as a trajectory.

        Parameters
        ----------
        topology : Topology
            The system topology
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        trajectory : Trajectory
            A trajectory object with 1 frame created from the file.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        xyz, time, cell_lengths, cell_angles = self.read(atom_indices=atom_indices)
        xyz = in_units_of(xyz, self.distance_unit, Trajectory._distance_unit,
                          inplace=True)
        cell_lengths = in_units_of(cell_lengths, self.distance_unit,
                                   Trajectory._distance_unit, inplace=True)

        return Trajectory(xyz=xyz, topology=topology, time=time,
                          unitcell_lengths=cell_lengths,
                          unitcell_angles=cell_angles)

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
            atom_slice = ensure_type(atom_indices, dtype=int, ndim=1,
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
        # These are no-ops.
        # coordinates = in_units_of(coordinates, None, 'angstroms')
        # time = in_units_of(time, None, 'picoseconds')
        # cell_lengths = in_units_of(cell_lengths, None, 'angstroms')
        # cell_angles = in_units_of(cell_angles, None, 'degrees')

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
        cell_angles = ensure_type(cell_angles, np.float64, 2, 'cell_angles',
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
        for i in range(self._n_atoms):
            acor = coordinates[0, i, :]
            self._handle.write(fmt % (acor[0], acor[1], acor[2]))
            if i % 2 == 1: self._handle.write('\n')
        if self._n_atoms % 2 == 1: self._handle.write('\n')
        if cell_lengths is not None:
            self._handle.write(fmt % (cell_lengths[0,0], cell_lengths[0,1],
                                      cell_lengths[0,2]))
            self._handle.write(fmt % (cell_angles[0,0], cell_angles[0,1],
                                      cell_angles[0,2]) + '\n')
        self._handle.flush()

    def __enter__(self):
        return self
    def __exit__(self, *exc_info):
        self.close()
    def close(self):
        if not self._closed and hasattr(self, '_handle'):
            self._handle.close()
        self._closed = True
    def __del__(self):
        self.close()
    def __len__(self):
        return 1 # All restarts have only 1 frame

@FormatRegistry.register_loader('.ncrst')
def load_ncrestrt(filename, top=None, atom_indices=None):
    """Load an AMBER NetCDF restart/inpcrd file. Since this file doesn't
    contain information to specify the topology, you need to supply a topology

    Parameters
    ----------
    filename : path-like
        name of the AMBER restart file
    top : {path-like, Trajectory, Topology}
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
    from mdtraj.core.trajectory import _parse_topology

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with AmberNetCDFRestartFile(filename) as f:
        return f.read_as_traj(topology, atom_indices=atom_indices)


@FormatRegistry.register_fileobject('.ncrst')
class AmberNetCDFRestartFile(object):
    """Interface for reading and writing AMBER NetCDF files. This is a file-like
    object, that supports both reading and writing depending on the `mode` flag.
    It implements the context manager protocol, so you can also use it with the
    python 'with' statement.

    Parameters
    ----------
    filename : path-like
        The name of the file to open
    mode : {'r', 'w'}, default='r'
        The mode in which to open the file. Valid options are 'r' or 'w' for
        'read' or 'write'
    force_overwrite : bool, default=False
        In write mode, if a file named `filename` already exists, clobber it and
        overwrite it
    """
    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=False):
        self._closed = True
        self._mode = mode
        if StrictVersion(import_('scipy.version').short_version) < StrictVersion('0.12.0'):
            raise ImportError('MDTraj NetCDF support requires scipy>=0.12.0. '
                              'You have %s' % import_('scipy.version').short_version)
        netcdf = import_('scipy.io').netcdf_file

        if mode not in ('r', 'w'):
            raise ValueError("mode must be one of ['r', 'w']")

        if mode == 'w' and not force_overwrite and os.path.exists(filename):
            raise IOError('"%s" already exists' % filename)

        # AMBER uses the NetCDF3 format, with 64 bit encodings, which for
        # scipy.io.netcdf_file is "version=2"
        self._handle = netcdf(filename, mode=mode, version=2)
        self._closed = False
        if mode == 'w':
            self._needs_initialization = True
        elif mode == 'r':
            self._needs_initialization = False
        else:
            raise RuntimeError()

    @property
    def n_atoms(self):
        self._validate_open()
        if self._needs_initialization:
            raise IOError('The file is uninitialized')
        return self._handle.dimensions['atom']

    @property
    def n_frames(self):
        return 1 # always 1 frame

    def _validate_open(self):
        if self._closed:
            raise IOError('The file is closed.')

    def read_as_traj(self, topology, atom_indices=None):
        """Read an AMBER ASCII restart file as a trajectory.

        Parameters
        ----------
        topology : Topology
            The system topology
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        trajectory : Trajectory
            A trajectory object with 1 frame created from the file.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        xyz, time, cell_lengths, cell_angles = self.read(atom_indices=atom_indices)
        xyz = in_units_of(xyz, self.distance_unit, Trajectory._distance_unit,
                          inplace=True)
        cell_lengths = in_units_of(cell_lengths, self.distance_unit,
                                   Trajectory._distance_unit, inplace=True)

        return Trajectory(xyz=xyz, topology=topology, time=time,
                          unitcell_lengths=cell_lengths,
                          unitcell_angles=cell_angles)

    def read(self, atom_indices=None):
        """Read data from an AMBER NetCDF restart file

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

        Notes
        -----
        If the file is not a NetCDF file with the appropriate convention, a
        TypeError is raised. If variables that are needed do not exist or if
        illegal values are passed in for parameters, ValueError is raised. If
        I/O errors occur, IOError is raised.
        """
        if self._mode != 'r':
            raise IOError('The file was opened in mode=%s. Reading is not '
                          'allowed.' % self._mode)

        if self._closed:
            raise IOError("The file has been closed!")

        if 'coordinates' not in self._handle.variables:
            raise ValueError('No coordinates found in the NetCDF file.')

        # Check that conventions are correct
        try:
            conventions = self._handle.Conventions.decode('ascii')
        except UnicodeDecodeError:
            raise TypeError('NetCDF file does not have correct Conventions')
        try:
            convention_version = self._handle.ConventionVersion.decode('ascii')
        except UnicodeDecodeError:
            raise ValueError('NetCDF file does not have correct ConventionVersion')
        except AttributeError:
            raise TypeError('NetCDF file does not have ConventionVersion')
        if (not hasattr(self._handle, 'Conventions') or
                conventions != 'AMBERRESTART'):
            raise TypeError('NetCDF file does not have correct Conventions')
        if convention_version != '1.0':
            raise ValueError('NetCDF restart has ConventionVersion %s. Only '
                             'Version 1.0 is supported.' % convention_version)
        if atom_indices is not None:
            atom_slice = ensure_type(atom_indices, dtype=int, ndim=1,
                                     name='atom_indices', warn_on_cast=False)
            if not np.all(atom_slice) >= 0:
                raise ValueError('Entries in atom_slice must be >= 0')
            coordinates = self._handle.variables['coordinates'][atom_slice, :]
        else:
            coordinates = self._handle.variables['coordinates'][:, :]

        # Get unit cell parameters
        if 'cell_lengths' in self._handle.variables:
            cell_lengths = self._handle.variables['cell_lengths'][:]
        else:
            cell_lengths = None
        if 'cell_angles' in self._handle.variables:
            cell_angles = self._handle.variables['cell_angles'][:]
        else:
            cell_angles = None

        if cell_lengths is None and cell_angles is not None:
            warnings.warn('cell_lengths were found, but no cell_angles')
        if cell_lengths is not None and cell_angles is None:
            warnings.warn('cell_angles were found, but no cell_lengths')

        if 'time' in self._handle.variables:
            time = self._handle.variables['time'].getValue()
        else:
            time = None

        # scipy.io.netcdf variables are mem-mapped, and are only backed by valid
        # memory while the file handle is open. This is _bad_ because we need to
        # support the user opening the file, reading the coordinates, and then
        # closing it, and still having the coordinates be a valid memory
        # segment.
        # https://github.com/mdtraj/mdtraj/issues/440
        if coordinates is not None and not coordinates.flags['WRITEABLE']:
            coordinates = np.array(coordinates, copy=True)
        if cell_lengths is not None and not cell_lengths.flags['WRITEABLE']:
            cell_lengths = np.array(cell_lengths, copy=True)
        if cell_angles is not None and not cell_angles.flags['WRITEABLE']:
            cell_angles = np.array(cell_angles, copy=True)

        # The leading frame dimension is missing on all of these arrays since
        # restart files have only one frame. Reshape them to add this extra
        # dimension
        coordinates = coordinates[np.newaxis,:]
        if cell_lengths is not None:
            cell_lengths = cell_lengths[np.newaxis,:]
        if cell_angles is not None:
            cell_angles = cell_angles[np.newaxis,:]
        if time is not None:
            time = np.asarray([time,])

        return coordinates, time, cell_lengths, cell_angles

    def write(self, coordinates, time=None, cell_lengths=None,
              cell_angles=None):
        """Write one frame of a MD trajectory to disk in the AMBER NetCDF
        restart file format.

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

        Notes
        -----
        You must only have one frame to write to this file.
        """
        if self._mode != 'w':
            raise IOError('The file was opened in mode=%s. Writing not allowed.'
                          % self._mode)
        if not self._needs_initialization:
            # Must have already been written -- can only write once
            raise RuntimeError('NetCDF restart file has already been written '
                               '-- can only write one frame to restart files.')
        # these are no-ops
        # coordinates = in_units_of(coordinates, None, 'angstroms')
        # time = in_units_of(time, None, 'picoseconds')
        # cell_lengths = in_units_of(cell_lengths, None, 'angstroms')
        # cell_angles = in_units_of(cell_angles, None, 'degrees')

        # typecheck all of the input arguments rigorously
        coordinates = ensure_type(coordinates, np.float32, 3, 'coordinates',
                                  length=None, can_be_none=False,
                                  shape=(1,None,3), warn_on_cast=False,
                                  add_newaxis_on_deficient_ndim=True)
        n_frames, n_atoms = coordinates.shape[0], coordinates.shape[1]
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
        cell_angles = ensure_type(cell_angles, np.float64, 2, 'cell_angles',
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

        self._initialize_headers(n_atoms=n_atoms,
                                 set_coordinates=True,
                                 set_time=(time is not None),
                                 set_cell=(cell_lengths is not None))
        self._needs_initialization = False

        # Write the time, coordinates, and box info
        if time is not None:
            self._handle.variables['time'][0] = float(time)
        self._handle.variables['coordinates'][:,:] = coordinates[0,:,:]
        if cell_lengths is not None:
            self._handle.variables['cell_angles'][:] = cell_angles[0,:]
            self._handle.variables['cell_lengths'][:] = cell_lengths[0,:]
        self.flush()

    def _initialize_headers(self, n_atoms, set_coordinates, set_time, set_cell):
        """Initialize the headers and convention properties of the NetCDF
        restart file
        """
        ncfile = self._handle
        ncfile.Conventions = 'AMBERRESTART'
        ncfile.ConventionVersion = "1.0"
        ncfile.title = 'NetCDF Restart file written by MDTraj w/out velocities'
        ncfile.application = 'Omnia'
        ncfile.program = 'MDTraj'
        ncfile.programVersion = version.short_version
        # Dimensions
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', n_atoms)
        if set_cell:
            ncfile.createDimension('cell_spatial', 3)
            ncfile.createDimension('label', 5)
            ncfile.createDimension('cell_angular', 3)
        if set_time:
            ncfile.createDimension('time', 1)
        # Variables
        v = ncfile.createVariable('spatial', 'c', ('spatial',))
        v[:] = np.asarray(list('xyz'))
        v = ncfile.createVariable('coordinates', 'd', ('atom', 'spatial'))
        v.units = 'angstrom'
        if set_cell:
            v = ncfile.createVariable('cell_angular', 'c',
                                      ('cell_angular', 'label'))
            v[0] = np.asarray(list('alpha'))
            v[1] = np.asarray(list('beta '))
            v[2] = np.asarray(list('gamma'))
            v = ncfile.createVariable('cell_spatial', 'c', ('cell_spatial',))
            v[:] = np.asarray(list('abc'))
            v = ncfile.createVariable('cell_lengths', 'd', ('cell_spatial',))
            v.units = 'angstrom'
            v = ncfile.createVariable('cell_angles', 'd', ('cell_angular',))
            v.units = 'degree'
        if set_time:
            v = ncfile.createVariable('time', 'd', ('time',))
            v.units = 'picosecond'
        self.flush()

    def __enter__(self):
        return self
    def __exit__(self, *exc_info):
        self.close()
    def close(self):
        if not self._closed and hasattr(self, '_handle'):
            self._handle.close()
        self._closed = True
    def __del__(self):
        self.close()
    def __len__(self):
        return 1 # All restarts have only 1 frame
    def flush(self):
        self._validate_open()
        if self._mode != 'w':
            raise IOError('Cannot flush a file opened for reading')
        self._handle.flush()
