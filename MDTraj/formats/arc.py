##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Lee-Ping Wang
# Contributors: Robert McGibbon and Jason Swails
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
from mdtraj.utils import cast_indices, in_units_of
from mdtraj.formats.registry import _FormatRegistry
from mdtraj.utils.six import string_types
from mdtraj.utils.six.moves import xrange

__all__ = ['ArcTrajectoryFile', 'load_arc']


##############################################################################
# Classes
##############################################################################

class _EOF(IOError):
    pass


@_FormatRegistry.register_loader('.arc')
@_FormatRegistry.register_loader('.xyz')
def load_arc(filename, top=None, stride=None, atom_indices=None):
    """Load a TINKER .arc file from disk.

    Parameters
    ----------
    filename : str
        String filename of TINKER .arc file.
    top : {str, Trajectory, Topology}
        The .arc format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
        information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.ArcTrajectoryFile :  Low level interface to TINKER .arc files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    # we make it not required in the signature, but required here. although this
    # is a little weird, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
#   if top is None:
#       raise ValueError('"top" argument is required for load_arc')

    if not isinstance(filename, string_types):
        raise TypeError('filename must be of type string for load_arc. '
            'you supplied %s' % type(filename))

    atom_indices = cast_indices(atom_indices)

    with ArcTrajectoryFile(filename) as f:
        xyz, abc, ang = f.read(stride=stride, atom_indices=atom_indices)
        in_units_of(xyz, f.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(abc, f.distance_unit, Trajectory._distance_unit, inplace=True)
        if top is None:
            topology = f.topology
        else:
            topology = _parse_topology(top)

    if atom_indices is not None:
        topology = topology.subset(atom_indices)

    time = np.arange(len(xyz))
    if stride is not None:
        # if we loaded with a stride, the Trajectories's time field should
        # respect that
        time *= stride

    t = Trajectory(xyz=xyz, topology=topology, time=time,
                   unitcell_lengths=abc,
                   unitcell_angles=ang)
    return t


@_FormatRegistry.register_fileobject('.arc')
@_FormatRegistry.register_fileobject('.xyz')
class ArcTrajectoryFile(object):
    """Interface for reading and writing to an TINKER archive files.
    (Note that the TINKER .xyz format is identical to this.)  This is
    a file-like object, that both reading or writing depending on the
    `mode` flag. It implements the context manager protocol, so you
    can also use it with the python 'with' statement.

    The conventional units in the arc file is angstrom. The format only
    supports storing the cartesian coordinates and box lengths.

    Attributes
    ----------
    topology : Topology
        A single-chain, single-residue topology generated from the atom and bond
        information found in the TINKER archive/xyz file. It is only generated
        from the first member of the archive

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r'}
        The mode in which to open the file, only 'r' for read is supported.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?
    """

    distance_unit = 'angstroms'

    def __init__(self, filename, mode='r', force_overwrite=True):
        """Open an TINKER.arc file for reading/writing.
        """
        self._is_open = False
        self._filename = filename
        self._mode = mode
        self.topology = None

        if mode == 'w':
            raise ValueError('Writing TINKER .arc files is not supported at this time')

        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0

        if mode == 'r':
            # if n_atoms is None:
            #     raise ValueError('To open a mdcrd file in mode="r", you must '
            #                      'supply the number of atoms, "n_atoms"')
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self._fh = open(filename, 'r')
            self._is_open = True
        else:
            raise ValueError('mode must be "r". '
                             'you supplied "%s"' % mode)

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
        raise NotImplementedError()

    def close(self):
        """Close the .arc file"""
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

    def __len__(self):
        "Number of frames in the file"
        raise NotImplementedError()

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a TINKER .arc file.  Note that only the
        Cartesian coordinates are read in.  The .arc file also
        contains TINKER-specific numeric atom types and some bonding
        information, which we do not read in.

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

        coords = []
        lengths = []
        angles = []
        for i in frame_counter:
            try:
                coord, length, angle = self._read()
                if atom_indices is not None:
                    coord = coord[atom_indices, :]
            except _EOF:
                break

            coords.append(coord)
            lengths.append(length)
            angles.append(angle)

            for j in range(stride - 1):
                # throw away these frames
                self._read()

        coords = np.array(coords)
        if None in lengths:
            lengths = angles = None
        else:
            lengths = np.array(lengths)
            angles = np.array(angles)

        return coords, lengths, angles

    def _read(self):
        "Read a single frame"
        from mdtraj.core.topology import Topology
        from mdtraj.core.element import Element
        # Read in the number of atoms.
        line = self._fh.readline()
        if line == '':
            raise _EOF()

        self._n_atoms = int(line.split()[0])
        self._line_counter += 1

        coords = np.empty((self._n_atoms, 3), dtype=np.float32)
        bond_partners = [[] for i in xrange(self._n_atoms)]
        atom_names = ['' for i in xrange(self._n_atoms)]
        line = self._fh.readline()
        s = line.split()
        self._line_counter += 1
        # See if we have box info on this line or not
        cell_lengths = cell_angles = None
        if len(s) == 6:
            try:
                cell_lengths = np.asarray(
                                [float(s[0]), float(s[1]), float(s[2])]
                )
                cell_angles = np.asarray(
                                [float(s[3]), float(s[4]), float(s[5])]
                )
                line = self._fh.readline()
                self._line_counter += 1
            except ValueError:
                pass
        i = 0
        while i < self._n_atoms - 1:
            atom_names[i] = s[1]
            bond_partners[i] = [int(x) for x in s[6:]]
            coords[i,:] = [float(s[pos]) for pos in [2, 3, 4]]
            i += 1
            line = self._fh.readline()
            s = line.split()
            self._line_counter += 1
        # Now do the last atom
        atom_names[i] = s[1]
        bond_partners[i] = [int(x) for x in s[6:]]
        coords[i,:] = [float(s[pos]) for pos in [2, 3, 4]]
        # Now see if we have to build a topology
        if self.topology is None:
            self.topology = top = Topology()
            chain = top.add_chain()                # only 1 chain
            res = top.add_residue('RES', chain, 1) # only 1 residue
            for at in atom_names:
                # First get the element. Try for common 2-letter elements, then
                # use the first letter only (default to None if I can't find it)
                if at[:2].upper() in ('NA', 'CL', 'MG'):
                    elem = Element.getBySymbol(at[:2])
                else:
                    try:
                        elem = Element.getBySymbol(at[0])
                    except KeyError:
                        elem = None
                top.add_atom(at, elem, res)
            # Now add the bonds
            atoms = list(top.atoms)
            for i, bonds in enumerate(bond_partners):
                me = atoms[i]
                for b in bonds:
                    b -= 1
                    if b < i: continue
                    it = atoms[b]
                    top.add_bond(me, it)

        return coords, cell_lengths, cell_angles

    def write(self, xyz):
        """ The ArcTrajectoryFile does not have a write method,
        because TINKER .arc files have special numerical atom types
        which are not shared by any other trajectory file format.

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write.

        """

        raise RuntimeError('write() is not available for .arc files')

