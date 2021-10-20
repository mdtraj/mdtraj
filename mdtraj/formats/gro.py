"""gro.py: Used for loading Gromacs GRO files.
"""
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Robert McGibbon, Lee-Ping Wang, Peter Eastman
# Contributors: Jason Swails
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
#
# Portions of this code originate from the OpenMM molecular simulation
# toolkit, copyright (c) 2012 Stanford University and the Authors. Those
# portions are distributed under the following terms:
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

##############################################################################
# Imports
##############################################################################


import os
import sys
import itertools
from re import sub, match, findall
# import element as elem
import numpy as np
import warnings

import mdtraj as md
from mdtraj.utils import in_units_of, cast_indices, ensure_type
from mdtraj.formats import pdb
from mdtraj.core import element as elem
from mdtraj.formats.registry import FormatRegistry


##############################################################################
# Code
##############################################################################

@FormatRegistry.register_loader('.gro')
def load_gro(filename, stride=None, atom_indices=None, frame=None, top=None):
    """Load a GROMACS GRO file.

    Parameters
    ----------
    filename : path-like
        Path to the GRO file on disk.
    stride : int, default=None
        Only read every stride-th model from the file
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. These indices are zero-based.
    frame : int, optional
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.
    top : mdtraj.core.Topology, default=None
        if you give a topology as input the topology won't be parsed from the gro file
        it saves time if you have to parse a big number of files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    with GroTrajectoryFile(filename, 'r', top=top) as f:
        topology = f.topology
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None
        return f.read_as_traj(n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)


@FormatRegistry.register_fileobject('.gro')
class GroTrajectoryFile(object):
    """Interface for reading and writing to GROMACS GRO files.

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
        if you give a topology as input the topology won't be parsed from the gro file
        it saves time if you have to parse a big number of files

    Attributes
    ----------
    n_atoms : int
        The number of atoms in the file
    topology : md.Topology
        The topology. TODO(rmcgibbo) note about chain

    See Also
    --------
    load_gro : High-level wrapper that returns a ``md.Trajectory``
    """
    distance_unit = 'nanometers'

    def __init__(self, filename, mode='r', force_overwrite=True, top=None):
        self._open = False
        self._file = None
        self._mode = mode

        self.topology = top

        if mode == 'r':
            self._open = True
            self._frame_index = 0
            self._file = open(filename, 'r')
            try:
                if self.topology is None:
                    self.n_atoms, self.topology = self._read_topology()
                else:
                    self.n_atoms = self.topology.n_atoms
            finally:
                self._file.seek(0)
        elif mode == 'w':
            self._open = True
            if os.path.exists(filename) and not force_overwrite:
                raise IOError('"%s" already exists' % filename)
            self._frame_index = 0
            self._file = open(filename, 'w')
        else:
            raise ValueError("invalid mode: %s" % mode)


    def write(self, coordinates, topology, time=None, unitcell_vectors=None,
              precision=3):
        """Write one or more frames of a molecular dynamics trajectory to disk
        in the GROMACS GRO format.

        Parameters
        ----------
        coordinates : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of each atom, in units of nanometers.
        topology : mdtraj.Topology
            The Topology defining the model to write.
        time : np.ndarray, dtype=float32, shape=(n_frames), optional
            The simulation time corresponding to each frame, in picoseconds.
            If not supplied, the numbers 0..n_frames will be written.
        unitcell_vectors : np.ndarray, dtype=float32, shape=(n_frames, 3, 3), optional
            The periodic box vectors of the simulation in each frame, in nanometers.
        precision : int, optional
            The number of decimal places to print for coordinates. Default is 3
        """
        if not self._open:
            raise ValueError('I/O operation on closed file')
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')

        coordinates = ensure_type(coordinates, dtype=np.float32, ndim=3, name='coordinates', can_be_none=False, warn_on_cast=False)
        time = ensure_type(time, dtype=float, ndim=1, name='time', can_be_none=True, shape=(len(coordinates),), warn_on_cast=False)
        unitcell_vectors = ensure_type(unitcell_vectors, dtype=float, ndim=3, name='unitcell_vectors',
            can_be_none=True, shape=(len(coordinates), 3, 3), warn_on_cast=False)

        for i in range(coordinates.shape[0]):
            frame_time = None if time is None else time[i]
            frame_box = None if unitcell_vectors is None else unitcell_vectors[i]
            self._write_frame(coordinates[i], topology, frame_time, frame_box, precision)

    def read_as_traj(self, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from a gro file

        Parameters
        ----------
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
        topology = self.topology
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        coordinates, time, unitcell_vectors = self.read(stride=stride, atom_indices=atom_indices)
        if len(coordinates) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        coordinates = in_units_of(coordinates, self.distance_unit, Trajectory._distance_unit, inplace=True)
        unitcell_vectors = in_units_of(unitcell_vectors, self.distance_unit, Trajectory._distance_unit, inplace=True)

        traj = Trajectory(xyz=coordinates, topology=topology, time=time)
        traj.unitcell_vectors = unitcell_vectors
        return traj


    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a molecular dynamics trajectory in the GROMACS GRO
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
            The cartesian coordinates of the atoms, in units of nanometers.
        time : np.ndarray, None
            The time corresponding to each frame, in units of picoseconds, or
            None if no time information is present in the trajectory.
        unitcell_vectors : np.ndarray, shape=(n_frames, 3, 3)
            The box vectors in each frame, in units of nanometers
        """
        if not self._open:
            raise ValueError('I/O operation on closed file')
        if not self._mode == 'r':
            raise ValueError('file not opened for reading')

        coordinates = []
        unitcell_vectors = []
        time = []
        contains_time = True

        atom_indices = cast_indices(atom_indices)
        atom_slice = slice(None) if atom_indices is None else atom_indices

        if n_frames is None:
            frameiter = itertools.count()
        else:
            frameiter = range(n_frames)

        for i in frameiter:
            try:
                frame_xyz, frame_box, frame_time = self._read_frame()
                contains_time = contains_time and (frame_time is not None)
                coordinates.append(frame_xyz[atom_slice])
                unitcell_vectors.append(frame_box)
                time.append(frame_time)
            except StopIteration:
                break

        coordinates, unitcell_vectors, time = map(np.array, (coordinates, unitcell_vectors, time))

        if not contains_time:
            time = None
        else:
            time = time[::stride]

        return coordinates[::stride], time, unitcell_vectors[::stride]

    def _read_topology(self):
        if not self._open:
            raise ValueError('I/O operation on closed file')
        if not self._mode == 'r':
            raise ValueError('file not opened for reading')
        pdb.PDBTrajectoryFile._loadNameReplacementTables()

        n_atoms = None
        topology = md.Topology()
        chain = topology.add_chain()
        residue = None
        atomReplacements = {}

        # This is needed because sometimes
        # residue names get replaced and this
        # brings to wrong residue parsing
        old_resname = None
        for ln, line in enumerate(self._file):
            if ln == 1:
                n_atoms = int(line.strip())
            elif ln > 1 and ln < n_atoms + 2:
                (thisresnum, thisresname, thisatomname, thisatomnum) = \
                    [line[i*5:i*5+5].strip() for i in range(4)]
                thisresnum, thisatomnum = map(int, (thisresnum, thisatomnum))
                if residue is None or residue.resSeq != thisresnum or old_resname != thisresname:
                    if residue is not None and thisresnum == residue.resSeq:
                        warnings.warn("WARNING: two consecutive residues with same number (%s, %s)" % (thisresname, old_resname))
                    old_resname = thisresname
                    if thisresname in pdb.PDBTrajectoryFile._residueNameReplacements:
                        thisresname = pdb.PDBTrajectoryFile._residueNameReplacements[thisresname]
                    residue = topology.add_residue(thisresname, chain, resSeq=thisresnum)
                    if thisresname in pdb.PDBTrajectoryFile._atomNameReplacements:
                        atomReplacements = pdb.PDBTrajectoryFile._atomNameReplacements[thisresname]
                    else:
                        atomReplacements = {}

                thiselem = thisatomname
                if len(thiselem) > 1:
                    thiselem = thiselem[0] + sub('[A-Z0-9]','',thiselem[1:])
                try:
                    element = elem.get_by_symbol(thiselem)
                except KeyError:
                    element = elem.virtual
                if thisatomname in atomReplacements:
                    thisatomname = atomReplacements[thisatomname]

                topology.add_atom(thisatomname, element=element, residue=residue,
                                  serial=thisatomnum)
        topology.create_standard_bonds()
        return n_atoms, topology

    def _read_frame(self):
        if not self._open:
            raise ValueError('I/O operation on closed file')
        if not self._mode == 'r':
            raise ValueError('file not opened for reading')

        atomcounter = itertools.count()
        comment = None
        boxvectors = None
        topology = None
        xyz = np.zeros((self.n_atoms, 3), dtype=np.float32)

        got_line = False
        firstDecimalPos = None
        atomindex = -1
        for ln, line in enumerate(self._file):
            got_line = True
            if ln == 0:
                comment = line.strip()
                continue
            elif ln == 1:
                assert self.n_atoms == int(line.strip())
                continue
            if firstDecimalPos is None:
                try:
                    firstDecimalPos = line.index('.', 20)
                    secondDecimalPos = line.index('.', firstDecimalPos+1)
                except ValueError:
                    firstDecimalPos = secondDecimalPos = None
            crd = _parse_gro_coord(line, firstDecimalPos, secondDecimalPos)
            if crd is not None and atomindex < self.n_atoms - 1:
                atomindex = next(atomcounter)
                xyz[atomindex, :] = (crd[0], crd[1], crd[2])
            elif _is_gro_box(line) and ln == self.n_atoms + 2:
                sline = line.split()
                boxvectors = tuple([float(i) for i in sline])
                # the gro_box line comes at the end of the record
                break
            else:
                raise Exception("Unexpected line in .gro file: "+line)

        if not got_line:
            raise StopIteration()

        time = None
        if 't=' in comment:
            # title string (free format string, optional time in ps after 't=')
            time = float(findall('t= *(\d+\.\d+)',comment)[-1])

        # box vectors (free format, space separated reals), values: v1(x) v2(y)
        # v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be
        # omitted (they will be set to zero).
        box = [boxvectors[i] if i < len(boxvectors) else 0 for i in range(9)]
        unitcell_vectors = np.array([
            [box[0], box[3], box[4]],
            [box[5], box[1], box[6]],
            [box[7], box[8], box[2]]])

        return xyz, unitcell_vectors, time

    def _write_frame(self, coordinates, topology, time, box, precision):
        comment = 'Generated with MDTraj'
        if time is not None:
            comment += ', t= %s' % time

        varwidth = precision + 5
        fmt = '%%5d%%-5s%%5s%%5d%%%d.%df%%%d.%df%%%d.%df' % (
                varwidth, precision, varwidth, precision, varwidth, precision)
        assert topology.n_atoms == coordinates.shape[0]
        lines = [comment, ' %d' % topology.n_atoms]
        if box is None:
            box = np.zeros((3,3))

        for i in range(topology.n_atoms):
            atom = topology.atom(i)
            residue = atom.residue
            serial = atom.serial
            if serial is None:
                serial = atom.index
            if serial >= 100000:
                serial %= 100000
            lines.append(fmt % (residue.resSeq, residue.name, atom.name, serial,
                                coordinates[i, 0], coordinates[i, 1], coordinates[i, 2]))

        lines.append('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f' % (
            box[0,0], box[1,1], box[2,2],
            box[0,1], box[0,2], box[1,0],
            box[1,2], box[2,0], box[2,1]))

        self._file.write('\n'.join(lines))
        self._file.write('\n')

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
        raise NotImplementedError()

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return self._frame_index

    def close(self):
        "Close the file"
        if self._open:
            self._file.close()
            self._open = False

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

##############################################################################
# Utilities
##############################################################################


def _isint(word):
    """ONLY matches integers! If you have a decimal point? None shall pass!

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is an integer (only +/- sign followed by digits)

    """
    return match('^[-+]?[0-9]+$',word)

def _isfloat(word):
    """Matches ANY number; it can be a decimal, scientific notation, what have you
    CAUTION - this will also match an integer.

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is any number

    """
    return match('^[-+]?[0-9]*\.?[0-9]*([eEdD][-+]?[0-9]+)?$',word)

def _parse_gro_coord(line, firstDecimal, secondDecimal):
    """ Determines whether a line contains GROMACS data or not

    @param[in] line The line to be tested

    """
    if firstDecimal is None or secondDecimal is None:
        return None
    digits = secondDecimal - firstDecimal
    try:
        return tuple(float(line[20+i*digits:20+(i+1)*digits]) for i in range(3))
    except ValueError:
        return None

def _is_gro_box(line):
    """ Determines whether a line contains a GROMACS box vector or not

    @param[in] line The line to be tested

    """
    sline = line.split()
    if len(sline) == 9 and all([_isfloat(i) for i in sline]):
        return 1
    elif len(sline) == 3 and all([_isfloat(i) for i in sline]):
        return 1
    else:
        return 0
