##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2022 Stanford University and the Authors
#
# Authors: Peter Eastman, Robert McGibbon
# Contributors: Carlos Hernandez, Jason Swails
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
# Portions copyright (c) 2012 Stanford University and the Authors.
#
# Portions of this code originate from the OpenMM molecular simulation
# toolkit. Those portions are Copyright 2008-2012 Stanford University
# and Peter Eastman, and distributed under the following license:
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
##############################################################################

import os
import warnings

import numpy as np

from mdtraj.core.topology import Topology
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import cast_indices, in_units_of, open_maybe_zipped
from mdtraj.utils.unitcell import (
    box_vectors_to_lengths_and_angles,
    lengths_and_angles_to_box_vectors,
)

__all__ = ["load_pdbx", "PDBxTrajectoryFile"]


##############################################################################
# Code
##############################################################################


@FormatRegistry.register_loader(".pdbx")
@FormatRegistry.register_loader(".cif")
def load_pdbx(
    filename,
    stride=None,
    atom_indices=None,
    frame=None,
    no_boxchk=False,
    top=None,
):
    """Load a PDBx/mmCIF file from disk.

    Parameters
    ----------
    filename : path-like
        Path to the PDBx/mmCIF file on disk
    stride : int, default=None
        Only read every stride-th model from the file
    atom_indices : array_like, default=None
        If not None, then read only a subset of the atoms coordinates from the
        file. These indices are zero-based (not 1 based, as used by the PDBx/mmCIF
        format). So if you want to load only the first atom in the file, you
        would supply ``atom_indices = np.array([0])``.
    frame : int, default=None
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.
    no_boxchk : bool, default=False
        By default, a heuristic check based on the particle density will be
        performed to determine if the unit cell dimensions are absurd. If the
        particle density is >1000 atoms per nm^3, the unit cell will be
        discarded. This is done because all PDBx/mmCIF files from RCSB contain a ``cell``
        record, even if there are no periodic boundaries, and dummy values are
        filled in instead. This check will filter out those false unit cells and
        avoid potential errors in geometry calculations. Set this variable to
        ``True`` in order to skip this heuristic check.
    top : mdtraj.core.Topology, default=None
        if you give a topology as input the topology won't be parsed from the file

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    Examples
    --------
    >>> import mdtraj as md
    >>> pdbx = md.load_pdbx('2EQQ.pdbx')
    >>> print(pdbx)
    <mdtraj.Trajectory with 20 frames, 423 atoms at 0x110740a90>

    See Also
    --------
    mdtraj.PDBxTrajectoryFile : Low level interface to PDBx/mmCIF files
    """
    from mdtraj import Trajectory

    if not isinstance(filename, (str, os.PathLike)):
        raise TypeError(
            "filename must be of type string or path-like for load_pdb. " "you supplied %s" % type(filename),
        )

    atom_indices = cast_indices(atom_indices)

    with PDBxTrajectoryFile(filename, top=top) as f:
        atom_slice = slice(None) if atom_indices is None else atom_indices
        if frame is not None:
            coords = f.positions[[frame], atom_slice, :]
        else:
            coords = f.positions[::stride, atom_slice, :]
        assert coords.ndim == 3, "internal shape error"
        n_frames = len(coords)

        topology = f.topology
        if atom_indices is not None:
            # The input topology shouldn't be modified because
            # subset makes a copy inside the function
            topology = topology.subset(atom_indices)

        if f.unitcell_angles is not None and f.unitcell_lengths is not None:
            unitcell_lengths = np.array([f.unitcell_lengths] * n_frames)
            unitcell_angles = np.array([f.unitcell_angles] * n_frames)
        else:
            unitcell_lengths = None
            unitcell_angles = None

        in_units_of(coords, f.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(
            unitcell_lengths,
            f.distance_unit,
            Trajectory._distance_unit,
            inplace=True,
        )

    time = np.arange(len(coords))
    if frame is not None:
        time *= frame
    elif stride is not None:
        time *= stride

    traj = Trajectory(
        xyz=coords,
        time=time,
        topology=topology,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
    )

    if not no_boxchk and traj.unitcell_lengths is not None:
        # Some PDBx/mmCIF files do not *really* have a unit cell, but still
        # have a cell record with a dummy definition. These boxes are usually
        # tiny (e.g., 1 A^3), so check that the particle density in the unit
        # cell is not absurdly high. Standard water density is ~55 M, which
        # yields a particle density ~100 atoms per cubic nm. It should be safe
        # to say that no particle density should exceed 10x that.
        particle_density = traj.top.n_atoms / traj.unitcell_volumes[0]
        if particle_density > 1000:
            warnings.warn(
                "Unlikely unit cell vectors detected in PDB file likely "
                "resulting from a dummy CRYST1 record. Discarding unit "
                "cell vectors.",
                category=UserWarning,
            )
            traj._unitcell_lengths = traj._unitcell_angles = None

    return traj


@FormatRegistry.register_fileobject(".pdbx")
@FormatRegistry.register_fileobject(".cif")
class PDBxTrajectoryFile:
    """Interface for reading and writing PDBx/mmCIF files

    Parameters
    ----------
    filename : path-like
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?
    top : mdtraj.core.Topology, default=None
        if you give a topology as input the topology won't be parsed from the file

    Attributes
    ----------
    positions : np.ndarray, shape=(n_frames, n_atoms, 3)
    topology : mdtraj.Topology
    closed : bool

    See Also
    --------
    mdtraj.load_pdbx : High-level wrapper that returns a ``md.Trajectory``
    """

    distance_unit = "nanometers"

    def __init__(self, filename, mode="r", force_overwrite=True, top=None):
        self._open = False
        self._mode = mode
        from openmm.app import PDBxFile
        from openmm.unit import nanometers

        if mode == "r":
            self._open = True
            pdbx = PDBxFile(filename)
            if top is None:
                self._topology = Topology.from_openmm(pdbx.topology)
            else:
                self._topology = top
            positions = [
                pdbx.getPositions(asNumpy=True, frame=i).value_in_unit(nanometers) for i in range(pdbx.getNumFrames())
            ]
            self._positions = np.array(positions)
            vectors = pdbx.topology.getPeriodicBoxVectors()
            if vectors is not None:
                vectors = [np.array(v.value_in_unit(nanometers)) for v in vectors]
                l1, l2, l3, alpha, beta, gamma = box_vectors_to_lengths_and_angles(
                    *vectors,
                )
                self._unitcell_lengths = (l1, l2, l3)
                self._unitcell_angles = (alpha, beta, gamma)
            else:
                self._unitcell_lengths = None
                self._unitcell_angles = None
        elif mode == "w":
            self._open = True
            self._next_model = 0
            self._file = open_maybe_zipped(filename, "w", force_overwrite)
        else:
            raise ValueError("invalid mode: %s" % mode)

    def write(self, positions, topology, unitcell_lengths=None, unitcell_angles=None):
        """Write one frame of a molecular dynamics trajectory to disk in PDBx/mmCIF format.

        Parameters
        ----------
        positions : array_like
            The list of atomic positions to write.
        topology : mdtraj.Topology
            The Topology defining the model to write.
        unitcell_lengths : {tuple, None}
            Lengths of the three unit cell vectors, or None for a non-periodic system
        unitcell_angles : {tuple, None}
            Angles between the three unit cell vectors, or None for a non-periodic system
        """
        if not self._mode == "w":
            raise ValueError("file not opened for writing")
        from openmm.app import PDBxFile
        from openmm.unit import nanometers

        if self._next_model == 0:
            self._openmm_topology = topology.to_openmm()
            if unitcell_lengths is None:
                self._openmm_topology.setPeriodicBoxVectors(None)
            else:
                vectors = lengths_and_angles_to_box_vectors(
                    *unitcell_lengths[0],
                    *unitcell_angles[0],
                )
                self._openmm_topology.setPeriodicBoxVectors(vectors * nanometers)
            PDBxFile.writeHeader(self._openmm_topology, self._file)
            self._next_model = 1
        if len(positions.shape) == 3:
            positions = positions[0]
        PDBxFile.writeModel(
            self._openmm_topology,
            positions * nanometers,
            self._file,
            self._next_model,
        )
        self._next_model += 1

    @property
    def positions(self):
        """The cartesian coordinates of all of the atoms in each frame. Available when a file is opened in mode='r'"""
        return self._positions

    @property
    def topology(self):
        """The topology from this PDBx/mmCIF file. Available when a file is opened in mode='r'"""
        return self._topology

    @property
    def unitcell_lengths(self):
        """The unitcell lengths (3-tuple) in this PDBx/mmCIF file. May be None"""
        return self._unitcell_lengths

    @property
    def unitcell_angles(self):
        """The unitcell angles (3-tuple) in this PDBx/mmCIF file. May be None"""
        return self._unitcell_angles

    @property
    def closed(self):
        """Whether the file is closed"""
        return not self._open

    def close(self):
        """Close the PDBx/mmCIF file"""
        if self._mode == "w" and self._open:
            self._file.close()
        self._open = False

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __len__(self):
        "Number of frames in the file"
        if str(self._mode) != "r":
            raise NotImplementedError('len() only available in mode="r" currently')
        return len(self._positions)
