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

import os
import socket
import warnings
from datetime import datetime

import numpy as np

import mdtraj
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import cast_indices, ensure_type, import_, in_units_of

__all__ = ["NetCDFTrajectoryFile", "load_netcdf"]

##############################################################################
# classes
##############################################################################


@FormatRegistry.register_loader(".nc")
@FormatRegistry.register_loader(".netcdf")
@FormatRegistry.register_loader(".ncdf")
def load_netcdf(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load an AMBER NetCDF file. Since the NetCDF format doesn't contain
    information to specify the topology, you need to supply a topology

    Parameters
    ----------
    filename : path-like
        filename of AMBER NetCDF file.
    top : {str, Trajectory, Topology}
        The NetCDF format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not None, then read only a subset of the atoms coordinates from the
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
    from mdtraj.core.trajectory import _parse_topology

    if top is None:
        raise ValueError('"top" argument is required for load_netcdf')

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with NetCDFTrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(
            topology,
            n_frames=n_frames,
            atom_indices=atom_indices,
            stride=stride,
        )


@FormatRegistry.register_fileobject(".nc")
@FormatRegistry.register_fileobject(".netcdf")
@FormatRegistry.register_fileobject(".ncdf")
class NetCDFTrajectoryFile:
    """Interface for reading and writing to AMBER NetCDF files. This is a
    file-like object, that supports both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    Parameters
    ----------
    filename : path-like
        The name of the file to open
    mode : {'r', 'w'}, default='r'
        The mode in which to open the file. Valid options are 'r' and 'w' for
        'read' and 'write', respectively.
    force_overwrite : bool, default=False
        In write mode, if a file named `filename` already exists, clobber
        it and overwrite it.
    """

    distance_unit = "angstroms"

    def __init__(self, filename, mode="r", force_overwrite=True):
        self._closed = True  # is the file currently closed?
        self._mode = mode  # what mode were we opened in

        try:
            # import netcdf4 if it's available
            import netCDF4

            netcdf = netCDF4.Dataset

            # set input args for netCDF4
            input_args = {"format": "NETCDF3_64BIT", "clobber": force_overwrite}

        except ImportError:
            import scipy.io

            netcdf = scipy.io.netcdf_file

            warning_message = (
                "Warning: The 'netCDF4' Python package is not installed. MDTraj is using the 'scipy' "
                "implementation to read and write netCDF files,which can be significantly slower.\n"
                "For improved performance, consider installing netCDF4. See installation instructions at:\n"
                "https://unidata.github.io/netcdf4-python/#quick-install"
            )

            warnings.warn(warning_message)

            # input args for scipy.io.netcdf_file
            # AMBER uses the NetCDF3 format, with 64 bit encodings, which
            # for scipy.io.netcdf_file is "version=2"
            input_args = {"version": 2}

        if mode not in ["r", "w"]:
            raise ValueError("mode must be one of ['r', 'w']")

        if mode == "w" and not force_overwrite and os.path.exists(filename):
            raise OSError('"%s" already exists' % filename)

        self._handle = netcdf(filename, mode=mode, **input_args)
        self._closed = False

        # self._frame_index is the current frame that we're at in the
        #     file
        # self._needs_initialization indicates whether we need to set the
        #     global properties of the file. This is required before the first
        #     write operation on a new file

        if mode == "w":
            self._frame_index = 0
            self._needs_initialization = True
        elif mode == "r":
            self._frame_index = 0
            self._needs_initialization = False
        else:
            raise RuntimeError()

    @property
    def n_atoms(self):
        self._validate_open()
        if self._needs_initialization:
            raise OSError("The file is uninitialized.")

        try:
            # netCDF4 requires use of .size on dimensions
            # to get size as integer.
            return self._handle.dimensions["atom"].size
        except AttributeError:
            # scipy case
            return self._handle.dimensions["atom"]

    @property
    def n_frames(self):
        self._validate_open()
        if not self._needs_initialization:
            return self._handle.variables["coordinates"].shape[0]
        return 0

    def _validate_open(self):
        if self._closed:
            raise OSError("The file is closed.")

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from a NetCDF file

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

        xyz, time, cell_lengths, cell_angles = self.read(
            n_frames=n_frames,
            stride=stride,
            atom_indices=atom_indices,
        )
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        xyz = in_units_of(
            xyz,
            self.distance_unit,
            Trajectory._distance_unit,
            inplace=True,
        )
        cell_lengths = in_units_of(
            cell_lengths,
            self.distance_unit,
            Trajectory._distance_unit,
            inplace=True,
        )

        return Trajectory(
            xyz=xyz,
            topology=topology,
            time=time,
            unitcell_lengths=cell_lengths,
            unitcell_angles=cell_angles,
        )

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
            The specific indices of the atoms you'd like to retrieve. If not
            supplied, all of the atoms will be retrieved.

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
            The angles (\alpha, \beta, \\gamma) defining the unit cell for
            each frame, or None if  the information is not present in the file.
        """
        self._validate_open()
        if self._mode != "r":
            raise OSError(
                "The file was opened in mode=%s. Reading is not allowed." % self._mode,
            )

        if n_frames is None:
            n_frames = np.inf
        elif stride is not None:
            # 'n_frames' frames should be read in total
            n_frames *= stride

        total_n_frames = self.n_frames
        frame_slice = slice(
            self._frame_index,
            self._frame_index + min(n_frames, total_n_frames),
            stride,
        )
        if self._frame_index >= total_n_frames:
            # just return something that'll look like len(xyz) == 0
            # this is basically just an alternative to throwing an indexerror
            return np.array([]), None, None, None

        if atom_indices is None:
            # get all of the atoms
            atom_slice = slice(None)
        else:
            atom_slice = ensure_type(
                atom_indices,
                dtype=int,
                ndim=1,
                name="atom_indices",
                warn_on_cast=False,
            )
            if not np.all(atom_slice < self.n_atoms):
                raise ValueError(
                    "As a zero-based index, the entries in "
                    "atom_indices must all be less than the number of atoms "
                    "in the trajectory, %d" % self.n_atoms,
                )
            if not np.all(atom_slice >= 0):
                raise ValueError(
                    "The entries in atom_indices must be greater " "than or equal to zero",
                )

        if "coordinates" in self._handle.variables:
            coordinates = self._handle.variables["coordinates"][
                frame_slice,
                atom_slice,
                :,
            ]
        else:
            raise ValueError(
                "No coordinates found in the NetCDF file. The only "
                "variables in the file were %s" % self._handle.variables.keys(),
            )

        if "time" in self._handle.variables:
            time = self._handle.variables["time"][frame_slice]
        else:
            time = None

        if "cell_lengths" in self._handle.variables:
            cell_lengths = self._handle.variables["cell_lengths"][frame_slice]
        else:
            cell_lengths = None

        if "cell_angles" in self._handle.variables:
            cell_angles = self._handle.variables["cell_angles"][frame_slice]
        else:
            cell_angles = None

        if cell_lengths is None and cell_angles is not None:
            warnings.warn("cell_lengths were found, but no cell_angles")
        if cell_lengths is not None and cell_angles is None:
            warnings.warn("cell_angles were found, but no cell_lengths")

        self._frame_index = self._frame_index + min(n_frames, total_n_frames)

        # scipy.io.netcdf variables are mem-mapped, and are only backed
        # by valid memory while the file handle is open. This is _bad_.
        # because we need to support the user opening the file, reading
        # the coordinates, and then closing it, and still having the
        # coordinates be a valid memory segment.
        # https://github.com/rmcgibbo/mdtraj/issues/440
        if coordinates is not None and not coordinates.flags["WRITEABLE"]:
            coordinates = np.array(coordinates, copy=True)
        if time is not None and not time.flags["WRITEABLE"]:
            time = np.array(time, copy=True)
        if cell_lengths is not None and not cell_lengths.flags["WRITEABLE"]:
            cell_lengths = np.array(cell_lengths, copy=True)
        if cell_angles is not None and not cell_angles.flags["WRITEABLE"]:
            cell_angles = np.array(cell_angles, copy=True)

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
            The angles (\alpha, \beta, \\gamma) defining the unit cell for
            each frame.

        Notes
        -----
        If the input arrays are of dimension deficient by one, for example
        if the coordinates array is two dimensional, the time is a single
        scalar or cell_lengths and cell_angles are a 1d array of length three,
        that is okay. You'll simply be saving a single frame.
        """
        self._validate_open()
        if self._mode not in ["w", "ws", "a", "as"]:
            raise OSError(
                "The file was opened in mode=%s. Writing is not allowed." % self._mode,
            )

        coordinates = in_units_of(coordinates, None, "angstroms")
        time = in_units_of(time, None, "picoseconds")
        cell_lengths = in_units_of(cell_lengths, None, "angstroms")
        cell_angles = in_units_of(cell_angles, None, "degrees")

        # typecheck all of the input arguments rigorously
        coordinates = ensure_type(
            coordinates,
            np.float32,
            3,
            "coordinates",
            length=None,
            can_be_none=False,
            shape=(None, None, 3),
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        n_frames, n_atoms = coordinates.shape[0], coordinates.shape[1]

        time = ensure_type(
            time,
            np.float32,
            1,
            "time",
            length=n_frames,
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        cell_lengths = ensure_type(
            cell_lengths,
            np.float64,
            2,
            "cell_lengths",
            length=n_frames,
            can_be_none=True,
            shape=(n_frames, 3),
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        cell_angles = ensure_type(
            cell_angles,
            np.float64,
            2,
            "cell_angles",
            length=n_frames,
            can_be_none=True,
            shape=(n_frames, 3),
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )

        # are we dealing with a periodic system?
        if (cell_lengths is None and cell_angles is not None) or (cell_lengths is not None and cell_angles is None):
            provided, neglected = "cell_lengths", "cell_angles"
            if cell_lengths is None:
                provided, neglected = neglected, provided
            raise ValueError(
                f'You provided the variable "{provided}", but neglected to '
                f'provide "{neglected}". They either BOTH must be provided, or '
                "neither. Having one without the other is meaningless",
            )

        if self._needs_initialization:
            self._initialize_headers(
                n_atoms=n_atoms,
                set_coordinates=True,
                set_time=(time is not None),
                set_cell=(cell_lengths is not None and cell_angles is not None),
            )
            self._needs_initialization = False

        # this slice object says where we're going to put the data in the
        # arrays
        frame_slice = slice(self._frame_index, self._frame_index + n_frames)

        # deposit the data
        try:
            self._handle.variables["coordinates"][frame_slice, :, :] = coordinates
            if time is not None:
                self._handle.variables["time"][frame_slice] = time
            if cell_lengths is not None:
                self._handle.variables["cell_lengths"][frame_slice, :] = cell_lengths
            if cell_angles is not None:
                self._handle.variables["cell_angles"][frame_slice, :] = cell_angles
        except KeyError as e:
            raise ValueError(
                "The file that you're trying to save to doesn't " "contain the field %s." % str(e),
            )

        # check for missing attributes
        missing = None
        if time is None and "time" in self._handle.variables:
            missing = "time"
        elif cell_angles is None and "cell_angles" in self._handle.variables:
            missing = "cell_angles"
        elif cell_lengths is None and "cell_lengths" in self._handle.variables:
            missing = "cell_lengths"
        if missing is not None:
            raise ValueError(
                "The file that you're saving to expects each frame "
                "to contain %s information, but you did not supply it."
                "I don't allow 'ragged' arrays." % missing,
            )

        # update the frame index pointers. this should be done at the
        # end so that if anything errors out, we don't actually get here
        self._frame_index += n_frames

    def flush(self):
        "Write all buffered data in the to the disk file."
        self._validate_open()
        self._handle.sync()

    def _initialize_headers(self, set_coordinates, n_atoms, set_time, set_cell):
        """Initialize the NetCDF file according to the AMBER NetCDF Convention,
        Version 1.0, revision B.

        The convention is defined here: http://ambermd.org/netcdf/nctraj.xhtml
        """
        # Set attributes.
        setattr(
            self._handle,
            "title",
            f"CREATED at {datetime.now()} on {socket.gethostname()}",
        )
        setattr(self._handle, "application", "Omnia")
        setattr(self._handle, "program", "MDTraj")
        setattr(self._handle, "programVersion", mdtraj.__version__)
        setattr(self._handle, "Conventions", "AMBER")
        setattr(self._handle, "ConventionVersion", "1.0")

        # set the dimensions
        # unlimited number of frames in trajectory
        self._handle.createDimension("frame", 0)
        # number of spatial coordinates
        self._handle.createDimension("spatial", 3)
        # number of atoms
        self._handle.createDimension("atom", n_atoms)

        if set_cell:
            # three spatial coordinates for the length of the unit cell
            self._handle.createDimension("cell_spatial", 3)
            # three spatial coordinates for the angles that define the shape
            # of the unit cell
            self._handle.createDimension("cell_angular", 3)
            # length of the longest string used for a label
            self._handle.createDimension("label", 5)

            # Define variables to store unit cell data
            cell_lengths = self._handle.createVariable(
                "cell_lengths",
                "d",
                ("frame", "cell_spatial"),
            )
            setattr(cell_lengths, "units", "angstrom")
            cell_angles = self._handle.createVariable(
                "cell_angles",
                "d",
                ("frame", "cell_angular"),
            )
            setattr(cell_angles, "units", "degree")

            self._handle.createVariable("cell_spatial", "c", ("cell_spatial",))
            self._handle.variables["cell_spatial"][0] = "a"
            self._handle.variables["cell_spatial"][1] = "b"
            self._handle.variables["cell_spatial"][2] = "c"

            self._handle.createVariable("cell_angular", "c", ("cell_spatial", "label"))
            self._handle.variables["cell_angular"][0] = "alpha"
            self._handle.variables["cell_angular"][1] = "beta "
            self._handle.variables["cell_angular"][2] = "gamma"

        if set_time:
            # Define coordinates and snapshot times.
            frame_times = self._handle.createVariable("time", "f", ("frame",))
            setattr(frame_times, "units", "picosecond")

        if set_coordinates:
            frame_coordinates = self._handle.createVariable(
                "coordinates",
                "f",
                ("frame", "atom", "spatial"),
            )
            setattr(frame_coordinates, "units", "angstrom")

            self._handle.createVariable("spatial", "c", ("spatial",))
            self._handle.variables["spatial"][0] = "x"
            self._handle.variables["spatial"][1] = "y"
            self._handle.variables["spatial"][2] = "z"

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
            self._frame_index = self.n_frames + offset
        else:
            raise OSError("Invalid argument")

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
        if not self._closed and hasattr(self, "_handle"):
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
            raise ValueError("I/O operation on closed file")
        return self.n_frames
