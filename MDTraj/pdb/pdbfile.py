# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

"""
pdbfile.py: Used for loading PDB files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import print_function, division

__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import numpy as np
import xml.etree.ElementTree as etree
from copy import copy
from .pdbstructure import PdbStructure
from mdtraj.topology import Topology
from mdtraj.utils import ilen
from . import element as elem


class PDBTrajectoryFile(object):
    """Interface for reading and writing Protein Data Bank (PDB) files

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?

    Attributes
    ----------
    positions : np.ndarray, shape=(n_frames, n_atoms, 3)
    topology : mdtraj.Topology
    closed : bool

    Notes
    -----
    When writing pdb files, mdtraj follows the PDB3.0 standard as closely as
    possible. During *reading* however, we try to be more lenient. For instance,
    we will parse common nonstandard atom names during reading, and convert them
    into the standard names. The replacement table used by mdtraj is at
    {mdtraj_source}/pdb/data/pdbNames.xml.

    See Also
    --------
    mdtraj.load_pdb : High-level wrapper that returns a ``md.Trajectory``
    """
    distance_unit = 'angstroms'
    _residueNameReplacements = {}
    _atomNameReplacements = {}
    _chain_names = [chr(ord('A') + i) for i in range(26)]

    def __init__(self, filename, mode='r', force_overwrite=True):
        self._open = False
        self._file = None
        self._topology = None
        self._positions = None
        self._mode = mode

        if mode == 'r':
            PDBTrajectoryFile._loadNameReplacementTables()
            self._file = open(filename, 'r')
            self._read_models()
        elif mode == 'w':
            self._header_written = False
            self._footer_written = False
            if os.path.exists(filename) and not force_overwrite:
                raise IOError('"%s" already exists' % filename)
            self._file = open(filename, 'w')
        else:
            raise ValueError("invalid mode: %s" % mode)

        self._open = True

    def write(self, positions, topology, modelIndex=None, unitcell_lengths=None, unitcell_angles=None):
        """Write a PDB file to disk

        Parameters
        ----------
        positions : array_like
            The list of atomic positions to write.
        topology : mdtraj.Topology
            The Topology defining the model to write.
        modelIndex : {int, None}
            If not None, the model will be surrounded by MODEL/ENDMDL records
            with this index
        unitcell_lengths : {tuple, None}
            Lengths of the three unit cell vectors, or None for a non-periodic system
        unitcell_angles : {tuple, None}
            Angles between the three unit cell vectors, or None for a non-periodic system
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        if not self._header_written:
            self._write_header(unitcell_lengths, unitcell_angles)
            self._header_written = True

        if ilen(topology.atoms) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if np.any(np.isnan(positions)):
            raise ValueError('Particle position is NaN')
        if np.any(np.isinf(positions)):
            raise ValueError('Particle position is infinite')

        atomIndex = 1
        posIndex = 0
        if modelIndex is not None:
            print("MODEL     %4d" % modelIndex, file=self._file)
        for (chainIndex, chain) in enumerate(topology.chains):
            chainName = self._chain_names[chainIndex % len(self._chain_names)]
            residues = list(chain.residues)
            for (resIndex, res) in enumerate(residues):
                if len(res.name) > 3:
                    resName = res.name[:3]
                else:
                    resName = res.name
                for atom in res.atoms:
                    if len(atom.name) < 4 and atom.name[:1].isalpha() and (atom.element is None or len(atom.element.symbol) < 2):
                        atomName = ' '+atom.name
                    elif len(atom.name) > 4:
                        atomName = atom.name[:4]
                    else:
                        atomName = atom.name
                    coords = positions[posIndex]
                    if atom.element is not None:
                        symbol = atom.element.symbol
                    else:
                        symbol = ' '
                    line = "ATOM  %5d %-4s %3s %s%4d    %s%s%s  1.00  0.00          %2s  " % (
                        atomIndex % 100000, atomName, resName, chainName,
                        (resIndex + 1) % 10000, _format_83(coords[0]),
                        _format_83(coords[1]), _format_83(coords[2]), symbol)
                    assert len(line) == 80, 'Fixed width overflow detected'
                    print(line, file=self._file)
                    posIndex += 1
                    atomIndex += 1
                if resIndex == len(residues)-1:
                    print("TER   %5d      %3s %s%4d" % (atomIndex, resName, chainName, resIndex+1), file=self._file)
                    atomIndex += 1

        if modelIndex is not None:
            print("ENDMDL", file=self._file)

    def _write_header(self, unitcell_lengths, unitcell_angles):
        """Write out the header for a PDB file.

        Parameters
        ----------
        unitcell_lengths : {tuple, None}
            The lengths of the three unitcell vectors, ``a``, ``b``, ``c``
        unitcell_angles : {tuple, None}
            The angles between the three unitcell vectors, ``alpha``,
            ``beta``, ``gamma``
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')

        if unitcell_lengths is None and unitcell_angles is None:
            return
        if unitcell_lengths is not None and unitcell_angles is not None:
            if not len(unitcell_lengths) == 3:
                raise ValueError('unitcell_lengths must be length 3')
            if not len(unitcell_angles) == 3:
                raise ValueError('unitcell_angles must be length 3')
        else:
            raise ValueError('either unitcell_lengths and unitcell_angles'
                             'should both be spefied, or neither')

        box = list(unitcell_lengths) + list(unitcell_angles)
        assert len(box) == 6

        print("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % tuple(box), file=self._file)

    def _write_footer(self):
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        print("END", file=self._file)
        self._footer_written = True

    @classmethod
    def set_chain_names(cls, values):
        """Set the cycle of chain names used when writing PDB files

        When writing PDB files, PDBTrajectoryFile translates each chain's
        index into a name -- the name is what's written in the file. By
        default, chains are named with the letters A-Z.

        Parameters
        ----------
        values : list
            A list of chacters (strings of length 1) that the PDB writer will
            cycle through to choose chain names.
        """
        for item in values:
            if not isinstance(item, basestring) and len(item) == 1:
                raise TypeError('Names must be a single character string')
        cls._chain_names = values

    @property
    def positions(self):
        """The cartesian coordinates of all of the atoms in each frame. Available when a file is opened in mode='r'
        """
        return self._positions

    @property
    def topology(self):
        """The topology from this PDB file. Available when a file is opened in mode='r'
        """
        return self._topology

    @property
    def unitcell_lengths(self):
        "The unitcell lengths (3-tuple) in this PDB file. May be None"
        return self._unitcell_lengths

    @property
    def unitcell_angles(self):
        "The unitcell angles (3-tuple) in this PDB file. May be None"
        return self._unitcell_angles

    @property
    def closed(self):
        "Whether the file is closed"
        return not self._open

    def close(self):
        "Close the PDB file"
        if self._mode == 'w' and not self._footer_written:
            self._write_footer()
        if self._open:
            self._file.close()
        self._open = False

    def _read_models(self):
        if not self._mode == 'r':
            raise ValueError('file not opened for reading')

        self._topology = Topology()

        pdb = PdbStructure(self._file, load_all_models=True)

        atomByNumber = {}
        for chain in pdb.iter_chains():
            c = self._topology.add_chain()
            for residue in chain.iter_residues():
                resName = residue.get_name()
                if resName in PDBTrajectoryFile._residueNameReplacements:
                    resName = PDBTrajectoryFile._residueNameReplacements[resName]
                r = self._topology.add_residue(resName, c)
                if resName in PDBTrajectoryFile._atomNameReplacements:
                    atomReplacements = PDBTrajectoryFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}
                for atom in residue.atoms:
                    atomName = atom.get_name()
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]
                    atomName = atomName.strip()
                    element = atom.element
                    if element is None:
                        element = self._guess_element(atomName, residue)

                    newAtom = self._topology.add_atom(atomName, element, r)
                    atomByNumber[atom.serial_number] = newAtom

        # load all of the positions (from every model)
        _positions = []
        for model in pdb.iter_models(use_all_models=True):
            coords = []
            for chain in model.iter_chains():
                for residue in chain.iter_residues():
                    for atom in residue.atoms:
                        coords.append(atom.get_position())
            _positions.append(coords)
        self._positions = np.array(_positions)

        ## The atom positions read from the PDB file
        self._unitcell_lengths = pdb.get_unit_cell_lengths()
        self._unitcell_angles = pdb.get_unit_cell_angles()
        self._topology.create_standard_bonds()
        self._topology.create_disulfide_bonds(self.positions[0])

        # Add bonds based on CONECT records.
        connectBonds = []
        for connect in pdb.models[0].connects:
            i = connect[0]
            for j in connect[1:]:
                connectBonds.append((atomByNumber[i], atomByNumber[j]))
        if len(connectBonds) > 0:
            # Only add bonds that don't already exist.
            existingBonds = set(self._topology.bonds)
            for bond in connectBonds:
                if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                    self._topology.add_bond(bond[0], bond[1])
                    existingBonds.add(bond)

    @staticmethod
    def _loadNameReplacementTables():
        """Load the list of atom and residue name replacements."""
        if len(PDBTrajectoryFile._residueNameReplacements) == 0:
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', 'pdbNames.xml'))
            allResidues = {}
            proteinResidues = {}
            nucleicAcidResidues = {}
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                if name == 'All':
                    PDBTrajectoryFile._parseResidueAtoms(residue, allResidues)
                elif name == 'Protein':
                    PDBTrajectoryFile._parseResidueAtoms(residue, proteinResidues)
                elif name == 'Nucleic':
                    PDBTrajectoryFile._parseResidueAtoms(residue, nucleicAcidResidues)
            for atom in allResidues:
                proteinResidues[atom] = allResidues[atom]
                nucleicAcidResidues[atom] = allResidues[atom]
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                for id in residue.attrib:
                    if id == 'name' or id.startswith('alt'):
                        PDBTrajectoryFile._residueNameReplacements[residue.attrib[id]] = name
                if 'type' not in residue.attrib:
                    atoms = copy(allResidues)
                elif residue.attrib['type'] == 'Protein':
                    atoms = copy(proteinResidues)
                elif residue.attrib['type'] == 'Nucleic':
                    atoms = copy(nucleicAcidResidues)
                else:
                    atoms = copy(allResidues)
                PDBTrajectoryFile._parseResidueAtoms(residue, atoms)
                PDBTrajectoryFile._atomNameReplacements[name] = atoms

    def _guess_element(self, atom_name, residue):
        "Try to guess the element name"

        upper = atom_name.upper()
        if upper.startswith('CL'):
            element = elem.chlorine
        elif upper.startswith('NA'):
            element = elem.sodium
        elif upper.startswith('MG'):
            element = elem.magnesium
        elif upper.startswith('BE'):
            element = elem.beryllium
        elif upper.startswith('LI'):
            element = elem.lithium
        elif upper.startswith('K'):
            element = elem.potassium
        elif len(residue) == 1 and upper.startswith('CA'):
            element = elem.calcium

        # TJL has edited this. There are a few issues here. First,
        # parsing for the element is non-trivial, so I do my best
        # below. Second, there is additional parsing code in
        # pdbstructure.py, and I am unsure why it doesn't get used
        # here...
        elif len(residue) > 1 and upper.startswith('CE'):
            element = elem.carbon  # (probably) not Celenium...
        elif len(residue) > 1 and upper.startswith('CD'):
            element = elem.carbon  # (probably) not Cadmium...
        elif residue.name in ['TRP', 'ARG', 'GLN', 'HIS'] and upper.startswith('NE'):
            element = elem.nitrogen  # (probably) not Neon...
        elif residue.name in ['ASN'] and upper.startswith('ND'):
            element = elem.nitrogen  # (probably) not ND...
        elif residue.name == 'CYS' and upper.startswith('SG'):
            element = elem.sulfur  # (probably) not SG...
        else:
            try:
                element = elem.get_by_symbol(atom_name[0])
            except KeyError:
                try:
                    symbol = atom_name[0:2].strip().rstrip("AB0123456789").lstrip("0123456789")
                    element = elem.get_by_symbol(symbol)                                    
                except KeyError:
                    element = None

        return element

    @staticmethod
    def _parseResidueAtoms(residue, map):
        for atom in residue.findall('Atom'):
            name = atom.attrib['name']
            for id in atom.attrib:
                map[atom.attrib[id]] = name

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError"""
    if -999.999 < f < 9999.999:
        return '%8.3f' % f
    if -9999999 < f < 99999999:
        return ('%8.3f' % f)[:8]
    raise ValueError('coordinate "%s" could not be represnted '
                     'in a width-8 field' % f)
