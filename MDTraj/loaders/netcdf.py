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


"""
This module provides the ability to read and write AMBER NetCDF trajectories.

The code is heavily based on amber_netcdf_trajectory_tools.py by John Chodera.
"""

##############################################################################
# imports
##############################################################################

from __future__ import print_function, division
# stdlib
import os
import warnings

import numpy as np
from mdtraj import version
from mdtraj.registry import _FormatRegistry
from mdtraj.utils import ensure_type, import_, in_units_of, convert, cast_indices

__all__ = ['NetCDFTrajectoryFile', 'load_netcdf']

##############################################################################
# classes
##############################################################################

@_FormatRegistry.register_loader('.nc')
@_FormatRegistry.register_loader('.netcdf')
def load_netcdf(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load an AMBER NetCDF file. Since the NetCDF format doesn't contain
    information to specify the topology, you need to supply a topology

    Parameters
    ----------
    filename : str
        filename of AMBER NetCDF file.
    top : {str, Trajectory, Topology}
        The NetCDF format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it
        requires an extra copy, but will save memory.
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
    mdtraj.NetCDFTrajectoryFile :  Low level interface to NetCDF files
    """
    from mdtraj.trajectory import _parse_topology, Trajectory

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)
    if atom_indices is not None:
        topology = topology.subset(atom_indices)

    with NetCDFTrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            xyz, time, cell_lengths, cell_angles = f.read(n_frames=1, atom_indices=atom_indices)
        else:
            xyz, time, cell_lengths, cell_angles = f.read(stride=stride, atom_indices=atom_indices)

        convert(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)
        convert(cell_lengths, f.distance_unit, Trajectory._distance_unit, inplace=True)

    if isinstance(time, np.ma.masked_array) and np.all(time.mask):
        # if time is a masked array and all the entries are masked
        # then we just tread it as if we never found it
        time = None
    if isinstance(cell_lengths, np.ma.masked_array) and np.all(cell_lengths.mask):
        cell_lengths = None
    if isinstance(cell_angles, np.ma.masked_array) and np.all(cell_angles.mask):
        cell_angles = None

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time,
                            unitcell_lengths=cell_lengths,
                            unitcell_angles=cell_angles)
    return trajectory


@_FormatRegistry.register_fileobject('.nc')
@_FormatRegistry.register_fileobject('.netcdf')
class NetCDFTrajectoryFile(object):
    """Interface for reading and writing to AMBER NetCDF files. This is a
    file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    Parameters
    ----------
    filename : str
        The name of the file to open
    mode : {'r', 'w', 'a', 'ws', 'as'}, default='r'
        The mode in which to open the file. Valid options are 'r', 'w',
        and 'a' for 'read', 'write', and 'append' respectively. The modes
        'w' and 'a' may also be prefixed with 's' which turns off buffering.
    force_overwrite : bool, default=False
        In write mode, if a file named `filename` already exists, clobber
        it and overwrite it.
    """
    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=True):
        self._closed = True   # is the file currently closed?
        self._mode = mode      # what mode were we opened in
        netcdf = import_('netCDF4')

        if mode not in ['r', 'w', 'a', 'ws', 'as']:
            raise ValueError(("mode must be one of ['r', 'w', 'a', 'ws', 'as']"
                " 'r' indicates read, 'w' indicates write, and 'a' indicates"
                " append. 'a' and 'w' can be appended with 's', which turns "
                " off buffering"))

        if mode in ['w', 'ws'] and not force_overwrite and os.path.exists(filename):
            raise IOError('"%s" already exists')

        # AMBER uses the NetCDF3 format, with 64 bit encodings
        self._handle = netcdf.Dataset(filename, mode=mode, format='NETCDF3_64BIT',
                                      clobber=force_overwrite)
        self._closed = False

        # self._frame_index is the current frame that we're at in the
        #     file
        # self._needs_initialization indicates whether we need to set the
        #     global properties of the file. This is required before the first
        #     write operation on a new file
        # self._n_atoms is the number of atoms in the file

        if mode in ['a', 'as']:
            self._frame_index = len(self._handle.dimensions['frame'])
            self._n_atoms = len(self._handle.dimensions['atom'])
            self._needs_initialization = False
        elif mode in ['w', 'ws']:
            self._frame_index = 0
            self._n_atoms = None
            # self._n_atoms will be set during _initialize_headers call
            self._needs_initialization = True
        elif mode == 'r':
            self._frame_index = 0
            self._n_atoms = len(self._handle.dimensions['atom'])
            self._needs_initialization = False
        else:
            raise RuntimeError()

    @property
    def n_atoms(self):
        self._validate_open()
        return len(self._handle.dimensions['atom'])

    @property
    def n_frames(self):
        self._validate_open()
        return len(self._handle.dimensions['frame'])

    def _validate_open(self):
        if self._closed:
            raise IOError('The file is closed.')

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a molecular dynamics trajectory in the AMBER NetCDF
        format.

        Parameters
        ----------
        n_frames : int, optional
            If n_frames is not None, the next n_frames of data from the file
            will be read. Otherwise, all of the frames in the file will be read.
        stride : int, optional
            If stride is not None, read only every stride-th frame from disk.
        atom_indices : np.ndarray, dtype=int, optional
            The specific indices of the atoms you'd like to retreive. If not
            supplied, all of the atoms will be retreived.

        Returns
        -------
        coordinates : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms, in units of angstroms.
        time : np.ndarray, None
            The time corresponding to each frame, in units of picoseconds, or
            None if no time information is present in the trajectory.
        cell_lengths : np.ndarray, None
            The lengths (a,b,c) of the unit cell for each frame, or None if
            the information is not present in the file.
        cell_angles : np.ndarray, None
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame, or None if  the information is not present in the file.
        """
        self._validate_open()
        if self._mode != 'r':
            raise IOError('The file was opened in mode=%s. Reading is not allowed.' % self._mode)

        if n_frames is None:
            n_frames = np.inf

        total_n_frames = len(self._handle.dimensions['frame'])
        frame_slice = slice(self._frame_index, self._frame_index + min(n_frames, total_n_frames), stride)
        if self._frame_index >= total_n_frames:
            # just return something that'll look like len(xyz) == 0
            # this is basically just an alternative to throwing an indexerror
            return np.array([]), None, None, None

        if atom_indices is None:
            # get all of the atoms
            atom_slice = slice(None)
        else:
            atom_slice = ensure_type(atom_indices, dtype=np.int, ndim=1,
                                     name='atom_indices', warn_on_cast=False)
            if not np.all(atom_slice < self._n_atoms):
                raise ValueError('As a zero-based index, the entries in '
                    'atom_indices must all be less than the number of atoms '
                    'in the trajectory, %d' % self._n_atoms)
            if not np.all(atom_slice >= 0):
                raise ValueError('The entries in atom_indices must be greater '
                    'than or equal to zero')

        if 'coordinates' in self._handle.variables:
            coordinates = self._handle.variables['coordinates'][frame_slice, atom_slice, :]
        else:
            raise ValueError('No coordinates found in the NetCDF file. The only'
                            'variables in the file were %s' % self._handle.variables.keys())

        if 'time' in self._handle.variables:
            time = self._handle.variables['time'][frame_slice]
        else:
            warnings.warn('No time information found in the NetCDF file')
            time = None

        if 'cell_lengths' in self._handle.variables:
            cell_lengths = self._handle.variables['cell_lengths'][frame_slice]
        else:
            cell_lengths = None

        if 'cell_angles' in self._handle.variables:
            cell_angles = self._handle.variables['cell_angles'][frame_slice]
        else:
            cell_angles = None

        if cell_lengths is None and cell_angles is not None:
            warnings.warn('cell_lengths were found, but no cell_angles')
        if cell_lengths is not None and cell_angles is None:
            warnings.warn('cell_angles were found, but no cell_lengths')

        self._frame_index = self._frame_index + min(n_frames, total_n_frames)

        return coordinates, time, cell_lengths, cell_angles

    def write(self, coordinates, time=None, cell_lengths=None, cell_angles=None):
        """Write one or more frames of a molecular dynamics trajectory to disk
        in the AMBER NetCDF format.

        Parameters
        ----------
        coordinates : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of each atom, in units of angstroms.
        time : np.ndarray, dtype=np.float32, shape=(n_frames), optional
            The time index corresponding to each frame, in units of picoseconds.
        cell_lengths : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The lengths (a,b,c) of the unit cell for each frame.
        cell_angles : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame.

        Notes
        -----
        If the input arrays are of dimension deficient by one, for example
        if the coordinates array is two dimensional, the time is a single
        scalar or cell_lengths and cell_angles are a 1d array of length three,
        that is okay. You'll simply be saving a single frame.
        """
        self._validate_open()
        if self._mode not in ['w', 'ws', 'a', 'as']:
            raise IOError('The file was opened in mode=%s. Writing is not allowed.' % self._mode)

        coordinates = in_units_of(coordinates, 'angstroms')
        time = in_units_of(time, 'picoseconds')
        cell_lengths = in_units_of(cell_lengths, 'angstroms')
        cell_angles = in_units_of(cell_angles, 'degrees')

        # typecheck all of the input arguments rigorously
        coordinates = ensure_type(coordinates, np.float32, 3, 'coordinates', length=None,
            can_be_none=False, shape=(None, None, 3), warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        n_frames, n_atoms = coordinates.shape[0], coordinates.shape[1]

        time = ensure_type(time, np.float32, 1, 'time', length=n_frames,
            can_be_none=True, warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, np.float64, 2, 'cell_lengths', length=n_frames,
            can_be_none=True, shape=(n_frames, 3), warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        cell_angles = ensure_type(cell_angles, np.float64, 2, 'cell_angles', length=n_frames,
            can_be_none=True, shape=(n_frames, 3), warn_on_cast=False, add_newaxis_on_deficient_ndim=True)

        # are we dealing with a periodic system?
        if (cell_lengths is None and cell_angles is not None) or (cell_lengths is not None and cell_angles is None):
            provided, neglected = 'cell_lengths', 'cell_angles'
            if cell_lengths is None:
                provided, neglected = neglected, provided
            raise ValueError('You provided the variable "%s", but neglected to '
                             'provide "%s". They either BOTH must be provided, or '
                             'neither. Having one without the other is meaningless' % (
                                provided, neglected))

        if self._needs_initialization:
            self._initialize_headers(n_atoms)
            self._needs_initialization = False

        # this slice object says where we're going to put the data in the
        # arrays
        frame_slice = slice(self._frame_index, self._frame_index + n_frames)

        # deposit the data
        self._handle.variables['coordinates'][frame_slice, :, :] = coordinates
        if time is not None:
            self._handle.variables['time'][frame_slice] = time
        if cell_lengths is not None:
            self._handle.variables['cell_lengths'][frame_slice, :] = cell_lengths
        if cell_angles is not None:
            self._handle.variables['cell_angles'][frame_slice, :] = cell_angles

        # update the frame index pointers. this should be done at the
        # end so that if anything errors out, we don't actually get here
        self._frame_index += n_frames

    def flush(self):
        "Write all buffered data in the to the disk file."
        self._validate_open()
        self._handle.sync()

    def _initialize_headers(self, n_atoms):
        """Initialize the NetCDF file according to the AMBER NetCDF Convention,
        Version 1.0, revision B.

        The convention is defined here: http://ambermd.org/netcdf/nctraj.html
        """
        self._n_atoms = n_atoms

        # Set attributes.
        setattr(self._handle, 'title', '')
        setattr(self._handle, 'application', 'AMBER')
        setattr(self._handle, 'program', 'MDTraj')
        setattr(self._handle, 'programVersion', version.short_version)
        setattr(self._handle, 'Conventions', 'AMBER')
        setattr(self._handle, 'ConventionVersion', '1.0')

        # set the dimensions
        # unlimited number of frames in trajectory
        self._handle.createDimension('frame', 0)
        # number of spatial coordinates
        self._handle.createDimension('spatial', 3)
        # number of atoms
        self._handle.createDimension('atom', n_atoms)
        # three spatial coordinates for the length of the unit cell
        self._handle.createDimension('cell_spatial', 3)
        # three spatial coordinates for the angles that define the shape
        # of the unit cell
        self._handle.createDimension('cell_angular', 3)
        # length of the longest string used for a label
        self._handle.createDimension('label', 5)

        # Define variables to store unit cell data
        self._handle.createVariable('cell_spatial', 'c', ('cell_spatial',))
        cell_angles = self._handle.createVariable('cell_angular', 'c', ('cell_spatial', 'label'))
        cell_lengths = self._handle.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
        setattr(cell_lengths, 'units', 'angstrom')
        cell_angles = self._handle.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
        setattr(cell_angles, 'units', 'degree')

        self._handle.variables['cell_spatial'][0] = 'x'
        self._handle.variables['cell_spatial'][1] = 'y'
        self._handle.variables['cell_spatial'][2] = 'z'

        self._handle.variables['cell_angular'][0] = 'alpha'
        self._handle.variables['cell_angular'][1] = 'beta '
        self._handle.variables['cell_angular'][2] = 'gamma'

        # Define coordinates and snapshot times.
        frame_times = self._handle.createVariable('time', 'f', ('frame',))
        setattr(frame_times, 'units', 'picosecond')
        frame_coordinates = self._handle.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
        setattr(frame_coordinates, 'units', 'angstrom')

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
        if whence == 0 and offset >= 0:
            self._frame_index = offset
        elif whence == 1:
            self._frame_index = self._frame_index + offset
        elif whence == 2 and offset <= 0:
            self._frame_index = len(self._handle.dimensions['frame']) + offset
        else:
            raise IOError('Invalid argument')

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return int(self._frame_index)

    def close(self):
        """Close the NetCDF file handle"""
        if not self._closed and hasattr(self, '_handle'):
            self._handle.close()

        self._closed = True

    def __enter__(self):
        # supports the context manager protocol
        return self

    def __exit__(self, *exc_info):
        # supports the context manager protocol
        self.close()

    def __del__(self):
        self.close()

    def __len__(self):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        return len(self._handle.dimensions['frame'])
