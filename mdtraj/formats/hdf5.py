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
This module implements the MDTraj HDF5 format described at
https://github.com/mdtraj/mdtraj/wiki/HDF5-Trajectory-Format
"""

##############################################################################
# Imports
##############################################################################

import operator
import os
import warnings
from collections import namedtuple

try:
    import simplejson as json
except ImportError:
    import json

# 3rd party
import numpy as np

import mdtraj.core.element as elem

# ours
import mdtraj
from mdtraj.core.topology import Topology
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import cast_indices, ensure_type, import_, in_units_of

__all__ = ["HDF5TrajectoryFile", "load_hdf5"]

Frames = namedtuple(
    "Frames",
    [
        "coordinates",
        "time",
        "cell_lengths",
        "cell_angles",
        "velocities",
        "kineticEnergy",
        "potentialEnergy",
        "temperature",
        "alchemicalLambda",
    ],
)

##############################################################################
# Code
##############################################################################


@FormatRegistry.register_loader(".h5")
@FormatRegistry.register_loader(".hdf5")
def load_hdf5(filename, stride=None, atom_indices=None, frame=None):
    """Load an MDTraj hdf5 trajectory file from disk.

    Parameters
    ----------
    filename : path-like
        Path of HDF Trajectory file.
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

    Examples
    --------
    >>> import mdtraj as md
    >>> traj = md.load_hdf5('output.h5')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    >>> traj2 = md.load_hdf5('output.h5', stride=2, top='topology.pdb')
    >>> print traj2
    <mdtraj.Trajectory with 250 frames, 423 atoms at 0x11136e410>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.HDF5TrajectoryFile :  Low level interface to HDF5 files
    """
    if not isinstance(filename, (str, os.PathLike)):
        raise TypeError(
            "filename must be of type path-like for load_lh5. " "you supplied %s" % type(filename),
        )

    atom_indices = cast_indices(atom_indices)

    with HDF5TrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None
        return f.read_as_traj(n_frames=n_frames, stride=stride, atom_indices=atom_indices)


@FormatRegistry.register_fileobject(".h5")
@FormatRegistry.register_fileobject(".hdf5")
class HDF5TrajectoryFile:
    """Interface for reading and writing to a MDTraj HDF5 molecular
    dynamics trajectory file, whose format is described
    `here <https://github.com/rmcgibbo/mdtraj/issues/36>`_.

    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The format is extremely flexible and high performance. It can hold a wide
    variety of information about a trajectory, including fields like the
    temperature and energies. Because it's built on the fantastic HDF5 library,
    it's easily extensible too.

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
    compression : {'zlib', None}
        Apply compression to the file? This will save space, and does not
        cost too many cpu cycles, so it's recommended.

    Attributes
    ----------
    root
    title
    application
    topology
    randomState
    forcefield
    reference
    constraints

    See Also
    --------
    mdtraj.load_hdf5 : High-level wrapper that returns a ``md.Trajectory``
    """

    distance_unit = "nanometers"

    def __init__(self, filename, mode="r", force_overwrite=True, compression="zlib"):
        self._open = False  # is the file handle currently open?
        self.mode = mode  # the mode in which the file was opened?

        if mode not in ["r", "w", "a"]:
            raise ValueError("mode must be one of ['r', 'w', 'a']")

        if mode == "w" and not force_overwrite and os.path.exists(filename):
            raise OSError('"%s" already exists' % filename)

        # import tables
        self.tables = import_("tables")

        if compression == "zlib":
            compression = self.tables.Filters(complib="zlib", shuffle=True, complevel=1)
        elif compression is None:
            compression = None
        else:
            raise ValueError('compression must be either "zlib" or None')

        self._handle = self._open_file(filename, mode=mode, filters=compression)
        self._open = True

        if mode == "w":
            # what frame are we currently reading or writing at?
            self._frame_index = 0
            # do we need to write the header information?
            self._needs_initialization = True
            if not os.fspath(filename).endswith(".h5"):
                warnings.warn("The .h5 extension is recommended.")

        elif mode == "a":
            try:
                self._frame_index = len(self._handle.root.coordinates)
                self._needs_initialization = False
            except self.tables.NoSuchNodeError:
                self._frame_index = 0
                self._needs_initialization = True
        elif mode == "r":
            self._frame_index = 0
            self._needs_initialization = False

    @property
    def root(self):
        """Direct access to the root group of the underlying Tables HDF5 file handle.

        This can be used for random or specific access to the underlying arrays
        on disk
        """
        _check_mode(self.mode, ("r", "a"))
        return self._handle.root

    #####################################################
    # title global attribute (optional)
    #####################################################

    @property
    def title(self):
        """User-defined title for the data represented in the file"""
        if hasattr(self._handle.root._v_attrs, "title"):
            return str(self._handle.root._v_attrs.title)
        return None

    @title.setter
    def title(self, value):
        """Set the user-defined title for the data represented in the file"""
        _check_mode(self.mode, ("w", "a"))
        self._handle.root._v_attrs.title = str(value)

    #####################################################
    # application global attribute (optional)
    #####################################################

    @property
    def application(self):
        "Suite of programs that created the file"
        if hasattr(self._handle.root._v_attrs, "application"):
            return str(self._handle.root._v_attrs.application)
        return None

    @application.setter
    def application(self, value):
        "Set the suite of programs that created the file"
        _check_mode(self.mode, ("w", "a"))
        self._handle.root._v_attrs.application = str(value)

    #####################################################
    # topology global attribute (optional, recommended)
    #####################################################

    @property
    def topology(self):
        """Get the topology out from the file

        Returns
        -------
        topology : mdtraj.Topology
            A topology object
        """
        try:
            raw = self._get_node("/", name="topology")[0]
            if not isinstance(raw, str):
                raw = raw.decode()
            topology_dict = json.loads(raw)
        except self.tables.NoSuchNodeError:
            return None

        topology = Topology()

        for chain_dict in sorted(topology_dict["chains"], key=operator.itemgetter("index")):
            chain = topology.add_chain()
            for residue_dict in sorted(chain_dict["residues"], key=operator.itemgetter("index")):
                try:
                    resSeq = residue_dict["resSeq"]
                except KeyError:
                    resSeq = None
                    warnings.warn("No resSeq information found in HDF file, defaulting to zero-based indices")
                try:
                    segment_id = residue_dict["segmentID"]
                except KeyError:
                    segment_id = ""
                residue = topology.add_residue(residue_dict["name"], chain, resSeq=resSeq, segment_id=segment_id)
                for atom_dict in sorted(residue_dict["atoms"], key=operator.itemgetter("index")):
                    try:
                        element = elem.get_by_symbol(atom_dict["element"])
                    except KeyError:
                        element = elem.virtual
                    topology.add_atom(atom_dict["name"], element, residue)

        atoms = list(topology.atoms)
        for index1, index2 in topology_dict["bonds"]:
            topology.add_bond(atoms[index1], atoms[index2])

        return topology

    @topology.setter
    def topology(self, topology_object):
        """Set the topology in the file

        Parameters
        ----------
        topology_object : mdtraj.Topology
            A topology object
        """
        _check_mode(self.mode, ("w", "a"))

        # we want to be able to handle the openmm Topology object
        # here too, so if it's not an mdtraj topology we'll just guess
        # that it's probably an openmm topology and convert
        if not isinstance(topology_object, Topology):
            topology_object = Topology.from_openmm(topology_object)

        try:
            topology_dict = {
                "chains": [],
                "bonds": [],
            }

            for chain in topology_object.chains:
                chain_dict = {
                    "residues": [],
                    "index": int(chain.index),
                }
                for residue in chain.residues:
                    residue_dict = {
                        "index": int(residue.index),
                        "name": str(residue.name),
                        "atoms": [],
                        "resSeq": int(residue.resSeq),
                        "segmentID": str(residue.segment_id),
                    }

                    for atom in residue.atoms:
                        try:
                            element_symbol_string = str(atom.element.symbol)
                        except AttributeError:
                            element_symbol_string = ""

                        residue_dict["atoms"].append(
                            {
                                "index": int(atom.index),
                                "name": str(atom.name),
                                "element": element_symbol_string,
                            },
                        )
                    chain_dict["residues"].append(residue_dict)
                topology_dict["chains"].append(chain_dict)

            for atom1, atom2 in topology_object.bonds:
                topology_dict["bonds"].append(
                    [
                        int(atom1.index),
                        int(atom2.index),
                    ],
                )

        except AttributeError as e:
            raise AttributeError(
                "topology_object fails to implement the"
                "chains() -> residue() -> atoms() and bond() protocol. "
                "Specifically, we encountered the following %s" % e,
            )

        # actually set the tables
        try:
            self._remove_node(where="/", name="topology")
        except self.tables.NoSuchNodeError:
            pass

        data = json.dumps(topology_dict)
        if not isinstance(data, bytes):
            data = data.encode("ascii")

        self._handle.create_array(where="/", name="topology", obj=[data])

    #####################################################
    # randomState global attribute (optional)
    #####################################################

    @property
    def randomState(self):
        "State of the creators internal random number generator at the start of the simulation"
        if hasattr(self._handle.root._v_attrs, "randomState"):
            return str(self._handle.root._v_attrs.randomState)
        return None

    @randomState.setter
    def randomState(self, value):
        "Set the state of the creators internal random number generator at the start of the simulation"
        _check_mode(self.mode, ("w", "a"))
        self._handle.root._v_attrs.randomState = str(value)

    #####################################################
    # forcefield global attribute (optional)
    #####################################################

    @property
    def forcefield(self):
        "Description of the hamiltonian used. A short, human readable string, like AMBER99sbildn."
        if hasattr(self._handle.root._v_attrs, "forcefield"):
            return str(self._handle.root._v_attrs.forcefield)
        return None

    @forcefield.setter
    def forcefield(self, value):
        "Set the description of the hamiltonian used. A short, human readable string, like AMBER99sbildn."
        _check_mode(self.mode, ("w", "a"))
        self._handle.root._v_attrs.forcefield = str(value)

    #####################################################
    # reference global attribute (optional)
    #####################################################

    @property
    def reference(self):
        "A published reference that documents the program or parameters used to generate the data"
        if hasattr(self._handle.root._v_attrs, "reference"):
            return str(self._handle.root._v_attrs.reference)
        return None

    @reference.setter
    def reference(self, value):
        "Set a published reference that documents the program or parameters used to generate the data"
        _check_mode(self.mode, ("w", "a"))
        self._handle.root._v_attrs.reference = str(value)

    #####################################################
    # constraints array (optional)
    #####################################################

    @property
    def constraints(self):
        """Constraints applied to the bond lengths

        Returns
        -------
        constraints : {None, np.array, dtype=[('atom1', '<i4'), ('atom2', '<i4'), ('distance', '<f4')])}
            A one dimensional array with the a int, int, float type giving
            the index of the two atoms involved in the constraints and the
            distance of the constraint. If no constraint information
            is in the file, the return value is None.
        """
        if hasattr(self._handle.root, "constraints"):
            return self._handle.root.constraints[:]
        return None

    @constraints.setter
    def constraints(self, value):
        """Set the constraints applied to bond lengths

        Returns
        -------
        valu : np.array, dtype=[('atom1', '<i4'), ('atom2', '<i4'), ('distance', '<f4')])
            A one dimensional array with the a int, int, float type giving
            the index of the two atoms involved in the constraints and the
            distance of the constraint.
        """
        _check_mode(self.mode, ("w", "a"))

        dtype = np.dtype(
            [
                ("atom1", np.int32),
                ("atom2", np.int32),
                ("distance", np.float32),
            ],
        )
        if not value.dtype == dtype:
            raise ValueError(
                "Constraints must be an array with dtype=%s. " "currently, I don't do any casting" % dtype,
            )

        if not hasattr(self._handle.root, "constraints"):
            self._create_table(
                where="/",
                name="constraints",
                description=dtype,
            )

        self._handle.root.constraints.truncate(0)
        self._handle.root.constraints.append(value)

    #####################################################
    # read/write methods for file-like behavior
    #####################################################

    def read_as_traj(self, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from the HDF5 file

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

        data = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(data) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(data.coordinates, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(data.cell_lengths, self.distance_unit, Trajectory._distance_unit, inplace=True)

        return Trajectory(
            xyz=data.coordinates,
            topology=topology,
            time=data.time,
            unitcell_lengths=data.cell_lengths,
            unitcell_angles=data.cell_angles,
        )

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

        Notes
        -----
        If you'd like more flexible access to the data, that is available by
        using the pytables group directly, which is accessible via the
        `root` property on this class.

        Returns
        -------
        frames : namedtuple
            The returned namedtuple will have the fields "coordinates", "time", "cell_lengths",
            "cell_angles", "velocities", "kineticEnergy", "potentialEnergy",
            "temperature" and "alchemicalLambda". Each of the fields in the
            returned namedtuple will either be a numpy array or None, dependening
            on if that data was saved in the trajectory. All of the data shall be
            n units of "nanometers", "picoseconds", "kelvin", "degrees" and
            "kilojoules_per_mole".
        """
        _check_mode(self.mode, ("r",))

        if n_frames is None:
            n_frames = np.inf
        if stride is not None:
            stride = int(stride)

        total_n_frames = len(self._handle.root.coordinates)
        frame_slice = slice(self._frame_index, min(self._frame_index + n_frames, total_n_frames), stride)
        if frame_slice.stop - frame_slice.start == 0:
            return []

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
            if not np.all(atom_slice < self._handle.root.coordinates.shape[1]):
                raise ValueError(
                    "As a zero-based index, the entries in "
                    "atom_indices must all be less than the number of atoms "
                    "in the trajectory, %d" % self._handle.root.coordinates.shape[1],
                )
            if not np.all(atom_slice >= 0):
                raise ValueError(
                    "The entries in atom_indices must be greater " "than or equal to zero",
                )

        def get_field(name, slice, out_units, can_be_none=True):
            try:
                node = self._get_node(where="/", name=name)
                data = node.__getitem__(slice)
                in_units = node.attrs.units
                if not isinstance(in_units, str):
                    in_units = in_units.decode()
                data = in_units_of(data, in_units, out_units)
                return data
            except self.tables.NoSuchNodeError:
                if can_be_none:
                    return None
                raise

        frames = Frames(
            coordinates=get_field(
                "coordinates",
                (frame_slice, atom_slice, slice(None)),
                out_units="nanometers",
                can_be_none=False,
            ),
            time=get_field("time", frame_slice, out_units="picoseconds"),
            cell_lengths=get_field("cell_lengths", (frame_slice, slice(None)), out_units="nanometers"),
            cell_angles=get_field("cell_angles", (frame_slice, slice(None)), out_units="degrees"),
            velocities=get_field(
                "velocities",
                (frame_slice, atom_slice, slice(None)),
                out_units="nanometers/picosecond",
            ),
            kineticEnergy=get_field("kineticEnergy", frame_slice, out_units="kilojoules_per_mole"),
            potentialEnergy=get_field("potentialEnergy", frame_slice, out_units="kilojoules_per_mole"),
            temperature=get_field("temperature", frame_slice, out_units="kelvin"),
            alchemicalLambda=get_field("lambda", frame_slice, out_units="dimensionless"),
        )

        self._frame_index += frame_slice.stop - frame_slice.start
        return frames

    def write(
        self,
        coordinates,
        time=None,
        cell_lengths=None,
        cell_angles=None,
        velocities=None,
        kineticEnergy=None,
        potentialEnergy=None,
        temperature=None,
        alchemicalLambda=None,
    ):
        """Write one or more frames of data to the file

        This method saves data that is associated with one or more simulation
        frames. Note that all of the arguments can either be raw numpy arrays
        or unitted arrays (with openmm.unit.Quantity). If the arrays are
        unittted, a unit conversion will be automatically done from the
        supplied units into the proper units for saving on disk. You won't have
        to worry about it.

        Furthermore, if you wish to save a single frame of simulation data, you
        can do so naturally, for instance by supplying a 2d array for the
        coordinates and a single float for the time. This "shape deficiency"
        will be recognized, and handled appropriately.

        Parameters
        ----------
        coordinates : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention, the
            lengths should be in units of nanometers.
        time : np.ndarray, shape=(n_frames,), optional
            You may optionally specify the simulation time, in picoseconds
            corresponding to each frame.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
            You may optionally specify the unitcell lengths.
            The length of the periodic box in each frame, in each direction,
            `a`, `b`, `c`. By convention the lengths should be in units
            of angstroms.
        cell_angles : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
            You may optionally specify the unitcell angles in each frame.
            Organized analogously to cell_lengths. Gives the alpha, beta and
            gamma angles respectively. By convention, the angles should be
            in units of degrees.
        velocities :  np.ndarray, shape=(n_frames, n_atoms, 3), optional
            You may optionally specify the cartesian components of the velocity
            for each atom in each frame. By convention, the velocities
            should be in units of nanometers / picosecond.
        kineticEnergy : np.ndarray, shape=(n_frames,), optional
            You may optionally specify the kinetic energy in each frame. By
            convention the kinetic energies should b in units of kilojoules per
            mole.
        potentialEnergy : np.ndarray, shape=(n_frames,), optional
            You may optionally specify the potential energy in each frame. By
            convention the kinetic energies should b in units of kilojoules per
            mole.
        temperature : np.ndarray, shape=(n_frames,), optional
            You may optionally specify the temperature in each frame. By
            convention the temperatures should b in units of Kelvin.
        alchemicalLambda : np.ndarray, shape=(n_frames,), optional
            You may optionally specify the alchemical lambda in each frame. These
            have no units, but are generally between zero and one.
        """
        _check_mode(self.mode, ("w", "a"))

        # these must be either both present or both absent. since
        # we're going to throw an error if one is present w/o the other,
        # lets do it now.
        if cell_lengths is None and cell_angles is not None:
            raise ValueError("cell_lengths were given, but no cell_angles")
        if cell_lengths is not None and cell_angles is None:
            raise ValueError("cell_angles were given, but no cell_lengths")

        # if the input arrays are openmm.unit.Quantities, convert them
        # into md units. Note that this acts as a no-op if the user doesn't
        # have openmm.unit installed (e.g. they didn't install OpenMM)
        coordinates = in_units_of(coordinates, None, "nanometers")
        time = in_units_of(time, None, "picoseconds")
        cell_lengths = in_units_of(cell_lengths, None, "nanometers")
        cell_angles = in_units_of(cell_angles, None, "degrees")
        velocities = in_units_of(velocities, None, "nanometers/picosecond")
        kineticEnergy = in_units_of(kineticEnergy, None, "kilojoules_per_mole")
        potentialEnergy = in_units_of(potentialEnergy, None, "kilojoules_per_mole")
        temperature = in_units_of(temperature, None, "kelvin")
        alchemicalLambda = in_units_of(alchemicalLambda, None, "dimensionless")

        # do typechecking and shapechecking on the arrays
        # this ensure_type method has a lot of options, but basically it lets
        # us validate most aspects of the array. Also, we can upconvert
        # on defficent ndim, which means that if the user sends in a single
        # frame of data (i.e. coordinates is shape=(n_atoms, 3)), we can
        # realize that. obviously the default mode is that they want to
        # write multiple frames at a time, so the coordinate shape is
        # (n_frames, n_atoms, 3)
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
        (
            n_frames,
            n_atoms,
        ) = coordinates.shape[0:2]
        time = ensure_type(
            time,
            dtype=np.float32,
            ndim=1,
            name="time",
            shape=(n_frames,),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        cell_lengths = ensure_type(
            cell_lengths,
            dtype=np.float32,
            ndim=2,
            name="cell_lengths",
            shape=(n_frames, 3),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        cell_angles = ensure_type(
            cell_angles,
            dtype=np.float32,
            ndim=2,
            name="cell_angles",
            shape=(n_frames, 3),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        velocities = ensure_type(
            velocities,
            dtype=np.float32,
            ndim=3,
            name="velocities",
            shape=(n_frames, n_atoms, 3),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        kineticEnergy = ensure_type(
            kineticEnergy,
            dtype=np.float32,
            ndim=1,
            name="kineticEnergy",
            shape=(n_frames,),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        potentialEnergy = ensure_type(
            potentialEnergy,
            dtype=np.float32,
            ndim=1,
            name="potentialEnergy",
            shape=(n_frames,),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        temperature = ensure_type(
            temperature,
            dtype=np.float32,
            ndim=1,
            name="temperature",
            shape=(n_frames,),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )
        alchemicalLambda = ensure_type(
            alchemicalLambda,
            dtype=np.float32,
            ndim=1,
            name="alchemicalLambda",
            shape=(n_frames,),
            can_be_none=True,
            warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True,
        )

        # if this is our first call to write(), we need to create the headers
        # and the arrays in the underlying HDF5 file
        if self._needs_initialization:
            self._initialize_headers(
                n_atoms=n_atoms,
                set_coordinates=True,
                set_time=(time is not None),
                set_cell=(cell_lengths is not None or cell_angles is not None),
                set_velocities=(velocities is not None),
                set_kineticEnergy=(kineticEnergy is not None),
                set_potentialEnergy=(potentialEnergy is not None),
                set_temperature=(temperature is not None),
                set_alchemicalLambda=(alchemicalLambda is not None),
            )
            self._needs_initialization = False

            # we need to check that that the entries that the user is trying
            # to save are actually fields in OUR file

        try:
            # try to get the nodes for all of the fields that we have
            # which are not None
            for name in [
                "coordinates",
                "time",
                "cell_angles",
                "cell_lengths",
                "velocities",
                "kineticEnergy",
                "potentialEnergy",
                "temperature",
            ]:
                contents = locals()[name]
                if contents is not None:
                    self._get_node(where="/", name=name).append(contents)
                if contents is None:
                    # for each attribute that they're not saving, we want
                    # to make sure the file doesn't explect it
                    try:
                        self._get_node(where="/", name=name)
                        raise AssertionError()
                    except self.tables.NoSuchNodeError:
                        pass

            # lambda is different, since the name in the file is lambda
            # but the name in this python function is alchemicalLambda
            name = "lambda"
            if alchemicalLambda is not None:
                self._get_node(where="/", name=name).append(alchemicalLambda)
            else:
                try:
                    self._get_node(where="/", name=name)
                    raise AssertionError()
                except self.tables.NoSuchNodeError:
                    pass

        except self.tables.NoSuchNodeError:
            raise ValueError(
                "The file that you're trying to save to doesn't "
                f"contain the field {name}. You can always save a new trajectory "
                "and have it contain this information, but I don't allow 'ragged' "
                f"arrays. If one frame is going to have {name} information, then I expect "
                "all of them to. So I can't save it for just these frames. Sorry "
                "about that :)",
            )
        except AssertionError:
            raise ValueError(
                "The file that you're saving to expects each frame "
                f"to contain {name} information, but you did not supply it."
                "I don't allow 'ragged' arrays. If one frame is going "
                f"to have {name} information, then I expect all of them to. ",
            )

        self._frame_index += n_frames
        self.flush()

    def _initialize_headers(
        self,
        n_atoms,
        set_coordinates,
        set_time,
        set_cell,
        set_velocities,
        set_kineticEnergy,
        set_potentialEnergy,
        set_temperature,
        set_alchemicalLambda,
    ):
        self._n_atoms = n_atoms

        self._handle.root._v_attrs.conventions = "Pande"
        self._handle.root._v_attrs.conventionVersion = "1.1"
        self._handle.root._v_attrs.program = "MDTraj"
        self._handle.root._v_attrs.programVersion = mdtraj.__version__
        self._handle.root._v_attrs.title = "title"

        # if the client has not the title attribute themselves, we'll
        # set it to MDTraj as a default option.
        if not hasattr(self._handle.root._v_attrs, "application"):
            self._handle.root._v_attrs.application = "MDTraj"

        # create arrays that store frame level informat
        if set_coordinates:
            self._create_earray(
                where="/",
                name="coordinates",
                atom=self.tables.Float32Atom(),
                shape=(0, self._n_atoms, 3),
            )
            self._handle.root.coordinates.attrs["units"] = "nanometers"

        if set_time:
            self._create_earray(
                where="/",
                name="time",
                atom=self.tables.Float32Atom(),
                shape=(0,),
            )
            self._handle.root.time.attrs["units"] = "picoseconds"

        if set_cell:
            self._create_earray(
                where="/",
                name="cell_lengths",
                atom=self.tables.Float32Atom(),
                shape=(0, 3),
            )
            self._create_earray(
                where="/",
                name="cell_angles",
                atom=self.tables.Float32Atom(),
                shape=(0, 3),
            )
            self._handle.root.cell_lengths.attrs["units"] = "nanometers"
            self._handle.root.cell_angles.attrs["units"] = "degrees"

        if set_velocities:
            self._create_earray(
                where="/",
                name="velocities",
                atom=self.tables.Float32Atom(),
                shape=(0, self._n_atoms, 3),
            )
            self._handle.root.velocities.attrs["units"] = "nanometers/picosecond"

        if set_kineticEnergy:
            self._create_earray(
                where="/",
                name="kineticEnergy",
                atom=self.tables.Float32Atom(),
                shape=(0,),
            )
            self._handle.root.kineticEnergy.attrs["units"] = "kilojoules_per_mole"

        if set_potentialEnergy:
            self._create_earray(
                where="/",
                name="potentialEnergy",
                atom=self.tables.Float32Atom(),
                shape=(0,),
            )
            self._handle.root.potentialEnergy.attrs["units"] = "kilojoules_per_mole"

        if set_temperature:
            self._create_earray(
                where="/",
                name="temperature",
                atom=self.tables.Float32Atom(),
                shape=(0,),
            )
            self._handle.root.temperature.attrs["units"] = "kelvin"

        if set_alchemicalLambda:
            self._create_earray(
                where="/",
                name="lambda",
                atom=self.tables.Float32Atom(),
                shape=(0,),
            )
            self._get_node("/", name="lambda").attrs["units"] = "dimensionless"

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
            self._frame_index = len(self._handle.root.coordinates) + offset
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

    def _validate(self):
        raise NotImplementedError()

        # check that all of the shapes are consistent
        # check that everything has units

    @property
    def _get_node(self):
        return self._handle.get_node

    @property
    def _create_earray(self):
        return self._handle.create_earray

    @property
    def _create_table(self):
        return self._handle.create_table

    @property
    def _remove_node(self):
        return self._handle.remove_node

    @property
    def _open_file(self):
        return self.tables.open_file

    def close(self):
        "Close the HDF5 file handle"
        if self._open:
            self._handle.close()
            self._open = False

    def flush(self):
        "Write all buffered data in the to the disk file."
        if self._open:
            self._handle.flush()

    def __del__(self):
        self.close()

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    def __len__(self):
        "Number of frames in the file"
        if not self._open:
            raise ValueError("I/O operation on closed file")
        return len(self._handle.root.coordinates)


def _check_mode(m, modes):
    if m not in modes:
        raise ValueError(
            "This operation is only available when a file " 'is open in mode="%s".' % m,
        )
