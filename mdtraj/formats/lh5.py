##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Robert McGibbon, Kyle A. Beauchamp
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

"""MSMBuilder2 "LH5" trajectory format."""

##############################################################################
# Imports
##############################################################################

import os
import sys
import warnings

import numpy as np

from mdtraj.core import element as elem
from mdtraj.formats.hdf5 import _check_mode
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import cast_indices, ensure_type, import_, in_units_of

MAXINT16 = np.iinfo(np.int16).max
MAXINT32 = np.iinfo(np.int32).max
DEFAULT_PRECISION = 1000

__all__ = ["LH5TrajectoryFile", "load_lh5"]

##############################################################################
# Utilities
##############################################################################


def _topology_from_arrays(AtomID, AtomNames, ChainID, ResidueID, ResidueNames):
    """Build topology object from the arrays stored in the lh5 file"""
    # Delayed import due to wacky recursive imports in compatibilty
    from mdtraj import Topology

    topology = Topology()

    # assert that the ChainID is just an array of empty strings, which appears
    # to be the case in our test systems for this legacy format
    if not np.all(chainid == "" for chainid in ChainID):
        raise NotImplementedError("Im not prepared to parse multiple chains")
    chain0 = topology.add_chain()

    # register the residues
    registered_residues = {}
    for i in np.argsort(ResidueID):
        residue_name = ResidueNames[i]
        if not isinstance(residue_name, str):
            residue_name = residue_name.decode()
        if ResidueID[i] not in registered_residues:
            res = topology.add_residue(residue_name, chain0)
            registered_residues[ResidueID[i]] = res

    # register the atoms
    for i in np.argsort(AtomID):
        atom_name = AtomNames[i]
        if not isinstance(atom_name, str):
            atom_name = atom_name.decode()
        element_symbol = atom_name.lstrip("0123456789")[0]

        try:
            element = elem.get_by_symbol(element_symbol)
        except KeyError:
            element = elem.virtual

        topology.add_atom(
            atom_name,
            element,
            registered_residues[ResidueID[i]],
        )

    topology.create_standard_bonds()
    return topology


def _convert_from_lossy_integers(X, precision=DEFAULT_PRECISION):
    """Implementation of the lossy compression used in Gromacs XTC using
    the pytables library.  Convert 16 bit integers into 32 bit floats."""
    X2 = X.astype("float32")
    X2 /= float(precision)
    return X2


def _convert_to_lossy_integers(X, precision=DEFAULT_PRECISION):
    """
    Implementation of the lossy compression used in Gromacs XTC using the pytables library.  Convert 32 bit floats
    into 16 bit integers.  These conversion functions have been optimized for memory use.  Further memory reduction
    would require an in-place astype() operation, which one could create using ctypes.
    """
    if np.max(X) * float(precision) < MAXINT16 and np.min(X) * float(precision) > -MAXINT16:
        X *= float(precision)
        Rounded = X.astype("int16")
        X /= float(precision)
    else:
        raise ValueError(
            "Data range too large for lh5. Try removing center of "
            "mass motion, check for 'blowing up, or use a different "
            "trajectory format",
        )
    return Rounded


##############################################################################
# Main code
##############################################################################


@FormatRegistry.register_loader(".lh5")
def load_lh5(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load an deprecated MSMBuilder2 LH5 trajectory file.

    Parameters
    ----------
    filename : path-like
        filename of AMBER NetCDF file.
    top : {path-like, Trajectory, Topology}
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

    See Also
    --------
    mdtraj.LH5TrajectoryFile :  Low level interface to LH5 files
    """
    atom_indices = cast_indices(atom_indices)
    with LH5TrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None
        return f.read_as_traj(
            n_frames=n_frames,
            stride=stride,
            atom_indices=atom_indices,
        )


@FormatRegistry.register_fileobject(".lh5")
class LH5TrajectoryFile:
    """Interface for reading and writing to a MSMBuilder2 "LH5" molecular
    dynamics trajectory file, a deprecated format.

    Parameters
    ----------
    filename : path-like
        Path to the file to open
    mode :  {'r, 'w'}
        Mode in which to open the file. 'r' is for reading and 'w' is for
        writing
    force_overwrite : bool
        In mode='w', how do you want to behave if a file by the name of `filename`
        already exists? if `force_overwrite=True`, it will be overwritten.
    """

    distance_unit = "nanometers"

    def __init__(self, filename, mode="r", force_overwrite=True):
        self._open = False
        self.filename = filename
        self.mode = mode
        if mode == "w" and not force_overwrite and os.path.exists(filename):
            raise OSError('"%s" already exists' % filename)
        # import tables
        self.tables = import_("tables")

        if mode == "w":
            print("Warning: The LH5 trajectory format is deprecated.", file=sys.stderr)
            # what frame are we currently reading or writing at?
            self._frame_index = 0
            # do we need to write the header information?
            self._needs_initialization = True
            if not os.fspath(filename).endswith(".lh5"):
                warnings.warn("The .lh5 extension is recommended.")
        elif mode == "r":
            self._frame_index = 0
            self._needs_initialization = False
        else:
            raise ValueError("mode must be one of ['r', 'w']")

        # Compression style of legacy MSMBuilder2 lh5 trajectory format
        compression = self.tables.Filters(
            complib="blosc",
            shuffle=True,
            complevel=1,
        )
        self._handle = self._open_file(
            filename,
            mode=mode,
            filters=compression,
        )
        self._open = True

    @property
    def topology(self):
        """Get the topology out from the file

        Returns
        -------
        topology : mdtraj.Topology
            A topology object
        """
        if np.all(self._handle.root.AtomID[:] == 0) and (
            np.all(self._handle.root.AtomNames[:] == b"") or np.all(self._handle.root.eAtomNames[:] == "")
        ):
            return None

        return _topology_from_arrays(
            self._handle.root.AtomID[:],
            self._handle.root.AtomNames[:],
            self._handle.root.ChainID[:],
            self._handle.root.ResidueID[:],
            self._handle.root.ResidueNames[:],
        )

    @topology.setter
    def topology(self, top):
        """Set the topology in the file

        Parameters
        ----------
        top : mdtraj.Topology
            A topology object
        """
        _check_mode(self.mode, ("w",))

        if self._needs_initialization:
            self._initialize_headers(top.n_atoms)
            self._needs_initialization = False

        top, bonds = top.to_dataframe()

        data = {
            "AtomID": top.index.values + 1,
            "AtomNames": top.name.values,
            "ResidueNames": top.resName.values,
            "ChainID": top.chainID.values,
            "ResidueID": top.resSeq.values + 1,
        }
        for key, val in data.items():
            node = self._get_node(where="/", name=key)[:] = val[:]
            node[:] = val[:]

    def read_as_traj(self, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from the LH5 file

        Parameters
        ----------
        n_frames : {int, None}
            The number of frames to read. If not supplied, all of the
            remaining frames will be read.
        stride : {int, None}
            By default all of the frames will be read, but you can pass this
            flag to read a subset of of the data by grabbing only every
            `stride`-th frame from disk.
        atom_indices : {int, None}
            By default all of the atom  will be read, but you can pass this
            flag to read only a subsets of the atoms for the `coordinates` and
            `velocities` fields. Note that you will have to carefully manage
            the indices and the offsets, since the `i`-th atom in the topology
            will not necessarily correspond to the `i`-th atom in your subset.

        Returns
        -------
        trajectory : Trajectory
            A trajectory object containing the loaded portion of the file.
        """
        _check_mode(self.mode, ("r",))

        from mdtraj.core.trajectory import Trajectory

        topology = self.topology
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        initial = int(self._frame_index)
        xyz = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        if stride is None:
            stride = 1
        time = (stride * np.arange(len(xyz))) + initial

        return Trajectory(xyz=xyz, topology=topology, time=time)

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read one or more frames of data from the file

        Parameters
        ----------
        n_frames : {int, None}
            The number of frames to read. If not supplied, all of the
            remaining frames will be read.
        stride : {int, None}
            By default all of the frames will be read, but you can pass this
            flag to read a subset of of the data by grabbing only every
            `stride`-th frame from disk.
        atom_indices : {int, None}
            By default all of the atom  will be read, but you can pass this
            flag to read only a subsets of the atoms for the `coordinates` and
            `velocities` fields. Note that you will have to carefully manage
            the indices and the offsets, since the `i`-th atom in the topology
            will not necessarily correspond to the `i`-th atom in your subset.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
            The cartesian coordinates, in nanometers
        """
        _check_mode(self.mode, ("r"))

        if n_frames is None:
            n_frames = np.inf
        if stride is not None:
            stride = int(stride)
        if atom_indices is None:
            atom_slice = slice(None)
        else:
            atom_slice = ensure_type(
                atom_indices,
                dtype=int,
                ndim=1,
                name="atom_indices",
                warn_on_cast=False,
            )

        total_n_frames = len(self._handle.root.XYZList)
        frame_slice = slice(
            self._frame_index,
            min(
                self._frame_index + n_frames,
                total_n_frames,
            ),
            stride,
        )
        if frame_slice.stop - frame_slice.start == 0:
            return np.array([], dtype=np.float32)

        xyz = self._handle.root.XYZList.__getitem__((frame_slice, atom_slice))
        if xyz.dtype == np.int16 or xyz.dtype == np.int32:
            xyz = _convert_from_lossy_integers(xyz)

        self._frame_index += frame_slice.stop - frame_slice.start
        return xyz

    def write(self, coordinates):
        """Write one or more frames of data to the file

        Parameters
        ----------
        coordinates : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms in every frame, in nanometers.
        """
        _check_mode(self.mode, ("w"))

        coordinates = ensure_type(
            coordinates,
            dtype=np.float32,
            ndim=3,
            name="coordinates",
            shape=(None, None, 3),
            can_be_none=False,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        if self._needs_initialization:
            self._initialize_headers(coordinates.shape[1])
            self._needs_initialization = False

        coordinates = _convert_to_lossy_integers(coordinates)
        self._get_node(where="/", name="XYZList").append(coordinates)

    def _initialize_headers(self, n_atoms):
        _check_mode(self.mode, ("w"))

        self._create_carray(
            where="/",
            name="AtomID",
            atom=self.tables.Int64Atom(),
            shape=(n_atoms,),
        )
        self._create_carray(
            where="/",
            name="AtomNames",
            atom=self.tables.StringAtom(itemsize=4),
            shape=(n_atoms,),
        )
        self._create_carray(
            where="/",
            name="ResidueNames",
            atom=self.tables.StringAtom(itemsize=4),
            shape=(n_atoms,),
        )
        self._create_carray(
            where="/",
            name="ChainID",
            atom=self.tables.StringAtom(itemsize=1),
            shape=(n_atoms,),
        )
        self._create_carray(
            where="/",
            name="ResidueID",
            atom=self.tables.Int64Atom(),
            shape=(n_atoms,),
        )
        self._create_earray(
            where="/",
            name="XYZList",
            atom=self.tables.Int16Atom(),
            shape=(0, n_atoms, 3),
        )

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
        _check_mode(self.mode, ("r",))

        if whence == 0 and offset >= 0:
            self._frame_index = offset
        elif whence == 1:
            self._frame_index = self._frame_index + offset
        elif whence == 2 and offset <= 0:
            self._frame_index = len(self._handle.root.XYZList) + offset
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
        "Close the HDF5 file handle"
        if self._open:
            self._handle.close()
            self._open = False

    def flush(self):
        "Write all buffered data in the to the disk file."
        if self._open:
            self._handle.flush()

    def __len__(self):
        "Number of frames in the file"
        if not self._open:
            raise ValueError("I/O operation on closed file")
        return len(self._handle.root.XYZList)

    def __del__(self):
        self.close()

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    # pytables 2/3 compatibility. pytables3 throws really annoying pending
    # deprecation warnings if you dont use the new method names
    @property
    def _open_file(self):
        return self.tables.open_file

    @property
    def _remove_node(self):
        return self._handle.remove_node

    @property
    def _create_carray(self):
        return self._handle.create_carray

    @property
    def _create_earray(self):
        return self._handle.create_earray

    @property
    def _get_node(self):
        return self._handle.get_node
