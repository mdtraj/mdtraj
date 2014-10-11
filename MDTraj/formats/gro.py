"""gro.py: Used for loading Gromacs GRO files.
"""
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Robert McGibbon, Lee-Ping Wang, Peter Eastman
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
from re import sub, match
# import element as elem
import numpy as np

from mdtraj import Topology
from mdtraj.formats import pdb
from mdtraj.core import element as elem
from mdtraj.formats.registry import _FormatRegistry


##############################################################################
# Code
##############################################################################

@_FormatRegistry.register_loader('.gro')
def load_gro(filename, stride=None, atom_indices=None, frame=None):
    """Load a GROMACS GRO file.
    
    Parameters
    ----------
    filename : str
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
    """
    raise NotImplementedError()


@_FormatRegistry.register_fileobject('.gro')
class GroTrajectoryFile(object):
    """Interface for reading and writing to GROMACS GRO files.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?

    See Also
    --------
    load_gro : High-level wrapper that returns a ``md.Trajectory``
    """

    def __init__(self, filename, mode='r', force_overwrite=True):
        self._open = False
        self._file = None
        self.n_atoms = 0
        
        if mode == 'r':
            pdb.PDBTrajectoryFile._loadNameReplacementTables()
            self._frame_index = 0
            self._file = open(filename, 'r')
            self._initialize_read()
        elif mode == 'w':
            if os.path.exists(filename) and not force_overwrite:
                raise IOError('"%s" already exists' % filename)
            self._frame_index = 0
            self._file = open(filename, 'w')
        else:
            raise ValueError("invalid mode: %s" % mode)

        self._open = True
    
    def _initialize_read(self):
        line0, line1 = [self._file.readline() for i in range(2)]
        self.n_atoms = int(line1.strip())
        self._file.seek(0)

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
            The cartesian coordinates of the atoms, in units of angstroms.
        time : np.ndarray, None
            The time corresponding to each frame, in units of picoseconds, or
            None if no time information is present in the trajectory.
        cell_vectors : np.ndarray
        """
        if self._frame_index == 0:
            xyz, unitcell_vectors, time, topology = self._read_frame(True)
            print xyz

    def _read_frame(self, parse_topology):
        atomcounter = itertools.count()
        comment = None
        boxvectors = None
        topology = None
        xyz = np.zeros((self.n_atoms, 3), dtype=np.float32)
        if parse_topology:
            topology = Topology()
            chain = topology.add_chain()
            residue = None

        for ln, line in enumerate(self._file):
            if ln == 0:
                comment = line.strip()
            elif ln == 1:
                assert self.n_atoms == int(line.strip())
            elif _is_gro_coord(line):
                atomindex = next(atomcounter)
                if parse_topology:
                    (thisresnum, thisresname, thisatomname, thisatomnum) = \
                        [line[i*5:i*5+5].strip() for i in range(4)]
                    thisresnum, thisatomnum = map(int, (thisresnum, thisatomnum))
                    if residue is None or residue.resSeq != thisresnum:
                        residue = topology.add_residue(thisresname, chain, resSeq=thisresnum)

                    thiselem = thisatomname
                    if len(thiselem) > 1:
                        thiselem = thiselem[0] + sub('[A-Z0-9]','',thiselem[1:])
                        try:
                            element = elem.get_by_symbol(thiselem)
                        except KeyError:
                            element = None

                    topology.add_atom(thisatomname, element=element, residue=residue,
                                      serial=thisatomnum)

                firstDecimalPos = line.index('.', 20)
                secondDecimalPos = line.index('.', firstDecimalPos+1)
                digits = secondDecimalPos-firstDecimalPos
                pos = [float(line[20+i*digits:20+(i+1)*digits]) for i in range(3)]
                xyz[atomindex, :] = (pos[0], pos[1], pos[2])
            elif _is_gro_box(line) and ln == self.n_atoms + 2:
                sline = line.split()
                boxvectors = tuple([float(i) for i in sline])
                # the gro_box line comes at the end of the record
                break
            else:
                raise Exception("Unexpected line in .gro file: "+line)
        
        time = None
        if 't=' in comment:
            # title string (free format string, optional time in ps after 't=')
            time = float(comment[comment.index('t=')+2:].strip())

        # box vectors (free format, space separated reals), values: v1(x) v2(y)
        # v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be
        # omitted (they will be set to zero).
        box = [boxvectors[i] if i < len(boxvectors) else 0 for i in range(9)]
        unitcell_vectors = np.array([
            [box[0], box[3], box[4]],
            [box[1], box[5], box[6]],
            [box[2], box[7], box[8]]])

        return xyz, unitcell_vectors, time, topology


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

def _is_gro_coord(line):
    """ Determines whether a line contains GROMACS data or not

    @param[in] line The line to be tested

    """
    sline = line.split()
    if len(sline) == 6 or len(sline) == 9:
        return all([_isint(sline[2]), _isfloat(sline[3]), _isfloat(sline[4]), _isfloat(sline[5])])
    elif len(sline) == 5 or len(sline) == 8:
        return all([_isint(line[15:20]), _isfloat(sline[2]), _isfloat(sline[3]), _isfloat(sline[4])])
    else:
        return 0

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
