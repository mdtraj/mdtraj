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
__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import sys
import math
import numpy as np
import xml.etree.ElementTree as etree
from copy import copy
from pdbstructure import PdbStructure
from mdtraj.topology import Topology
import element as elem


class PDBTrajectoryFile(object):
    """Interface for reading and writing Protein Data Bank (PDB) files
    """
    distance_unit = 'angstroms'
    _residueNameReplacements = {}
    _atomNameReplacements = {}

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

    def write(self, positions, topology, modelIndex):
        """Write a PDB file

        Parameters
        ----------
        positions : list
            The list of atomic positions to write.
        topology : Topology, optional
            The Topology defining the model to write.
        modelIndex : {int, None}
            If not None, the model will be surrounded by MODEL/ENDMDL records
            with this index
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        if not self._header_written:
            self._write_header(topology)
            self._header_written = True

        if len(list(topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if np.any(np.isnan(positions)):
            raise ValueError('Particle position is NaN')
        if np.any(np.isinf(positions)):
            raise ValueError('Particle position is infinite')

        atomIndex = 1
        posIndex = 0
        if modelIndex is not None:
            print >> self._file, "MODEL     %4d" % modelIndex
        for (chainIndex, chain) in enumerate(topology.chains()):
            chainName = chr(ord('A')+chainIndex%26)
            residues = list(chain.residues())
            for (resIndex, res) in enumerate(residues):
                if len(res.name) > 3:
                    resName = res.name[:3]
                else:
                    resName = res.name
                for atom in res.atoms():
                    if len(atom.name) < 4 and atom.name[:1].isalpha() and (atom.element is None or len(atom.element.symbol) < 2):
                        atomName = ' '+atom.name
                    elif len(atom.name) > 4:
                        atomName = atom.name[:4]
                    else:
                        atomName = atom.name
                    coords = positions[posIndex]
                    line = "ATOM  %5d %-4s %3s %s%4d    %s%s%s  1.00  0.00" % (
                        atomIndex%100000, atomName, resName, chainName,
                        (resIndex+1)%10000, _format_83(coords[0]),
                        _format_83(coords[1]), _format_83(coords[2]))
                    assert len(line) == 66, 'Fixed width overflow detected'
                    print >> self._file, line
                    posIndex += 1
                    atomIndex += 1
                if resIndex == len(residues)-1:
                    print >> self._file, "TER   %5d      %3s %s%4d" % (atomIndex, resName, chainName, resIndex+1)
                    atomIndex += 1

        if modelIndex is not None:
            print >> self._file, "ENDMDL"

    def _write_header(self, topology):
        """Write out the header for a PDB file.

        Parameters
        ----------
         topology : Topology
            The Topology defining the molecular system being written
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        boxSize = topology.getUnitCellDimensions()
        if boxSize is not None:
            print >>self._file, "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 " % boxSize

    def _write_footer(self):
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        print >>self._file, "END"
        self._footer_written = True

    @property
    def positions(self):
        return self._positions

    @property
    def topology(self):
        "Get the topology from this PDB file"
        return self._topology

    @property
    def closed(self):
        "True if the file is closed"
        return not self._open

    def close(self):
        "Close the file"
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
            c = self._topology.addChain()
            for residue in chain.iter_residues():
                resName = residue.get_name()
                if resName in PDBTrajectoryFile._residueNameReplacements:
                    resName = PDBTrajectoryFile._residueNameReplacements[resName]
                r = self._topology.addResidue(resName, c)
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
                        # Try to guess the element.

                        upper = atomName.upper()
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
                        elif( len( residue ) == 1 and upper.startswith('CA') ):
                            element = elem.calcium

                        # TJL has edited this. There are a few issues here. First,
                        # parsing for the element is non-trivial, so I do my best
                        # below. Second, there is additional parsing code in
                        # pdbstructure.py, and I am unsure why it doesn't get used
                        # here...
                        elif ( len( residue ) > 1 and upper.startswith('CE') ):
                            element = elem.carbon # (probably) not Celenium...
                        elif ( len( residue ) > 1 and upper.startswith('CD') ):
                            element = elem.carbon # (probably) not Cadmium...
                        elif ( residue.name in ['TRP', 'ARG', 'GLN', 'HIS'] and upper.startswith('NE') ):
                            element = elem.nitrogen # (probably) not Neon...
                        elif ( residue.name in ['ASN'] and upper.startswith('ND') ):
                            element = elem.nitrogen # (probably) not ND...
                        elif ( residue.name == 'CYS' and upper.startswith('SG') ):
                            element = elem.sulfur # (probably) not SG...
                        else:
                            try:
                                symbol = atomName[0:2].strip().rstrip("AB0123456789").lstrip("0123456789")
                                element = elem.get_by_symbol(symbol)
                            except KeyError:
                                try:
                                    element = elem.get_by_symbol(atomName[0])
                                except KeyError:
                                    pass

                    newAtom = self._topology.addAtom(atomName, element, r)
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
        #self.positions = np.array(coords)
        #print self.positions.shape
        self._topology.setUnitCellDimensions(pdb.get_unit_cell_dimensions())
        self._topology.createStandardBonds()
        self._topology.createDisulfideBonds(self.positions[0])

        # Add bonds based on CONECT records.
        connectBonds = []
        for connect in pdb.models[0].connects:
            i = connect[0]
            for j in connect[1:]:
                connectBonds.append((atomByNumber[i], atomByNumber[j]))
        if len(connectBonds) > 0:
            # Only add bonds that don't already exist.
            existingBonds = set(self._topology.bonds())
            for bond in connectBonds:
                if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                    self._topology.addBond(bond[0], bond[1])
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
