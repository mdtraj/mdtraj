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


class PDBFile(object):
    """PDBFile parses a Protein Data Bank (PDB) file and constructs a Topology and a set of atom positions from it.

    This class also provides methods for creating PDB files.  To write a file containing a single model, call
    writeFile().  You also can create files that contain multiple models.  To do this, first call writeHeader(),
    then writeModel() once for each model in the file, and finally writeFooter() to complete the file."""

    _residueNameReplacements = {}
    _atomNameReplacements = {}

    def __init__(self, file):
        """Load a PDB file.

        The atom positions and Topology can be retrieved by calling getPositions() and getTopology().

        Parameters:
         - file (string) the name of the file to load
        """
        top = Topology()
        coords = [];
        ## The Topology read from the PDB file
        self.topology = top

        # Load the PDB file

        pdb = PdbStructure(open(file))
        PDBFile._loadNameReplacementTables()

        # Build the topology

        # TJL QUESTION:
        # Isn't the following parsing supported in pdbstructure.py, which we just
        # called??

        atomByNumber = {}
        for chain in pdb.iter_chains():
            c = top.addChain()
            for residue in chain.iter_residues():
                resName = residue.get_name()
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}
                for atom in residue.atoms:
                    atomName = atom.get_name()
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]
                    atomName = atomName.strip()
                    element = None

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
                            except KeyError: pass

                    newAtom = top.addAtom(atomName, element, r)
                    atomByNumber[atom.serial_number] = newAtom
                    pos = atom.get_position()
                    coords.append(pos)
        ## The atom positions read from the PDB file
        self.positions = np.array(coords)
        self.topology.setUnitCellDimensions(pdb.get_unit_cell_dimensions())
        self.topology.createStandardBonds()
        self.topology.createDisulfideBonds(self.positions)

        # Add bonds based on CONECT records.

        connectBonds = []
        for connect in pdb.models[0].connects:
            i = connect[0]
            for j in connect[1:]:
                connectBonds.append((atomByNumber[i], atomByNumber[j]))
        if len(connectBonds) > 0:
            # Only add bonds that don't already exist.
            existingBonds = set(top.bonds())
            for bond in connectBonds:
                if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                    top.addBond(bond[0], bond[1])
                    existingBonds.add(bond)

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology

    def getPositions(self):
        """Get the atomic positions."""
        return self.positions

    @staticmethod
    def _loadNameReplacementTables():
        """Load the list of atom and residue name replacements."""
        if len(PDBFile._residueNameReplacements) == 0:
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', 'pdbNames.xml'))
            allResidues = {}
            proteinResidues = {}
            nucleicAcidResidues = {}
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                if name == 'All':
                    PDBFile._parseResidueAtoms(residue, allResidues)
                elif name == 'Protein':
                    PDBFile._parseResidueAtoms(residue, proteinResidues)
                elif name == 'Nucleic':
                    PDBFile._parseResidueAtoms(residue, nucleicAcidResidues)
            for atom in allResidues:
                proteinResidues[atom] = allResidues[atom]
                nucleicAcidResidues[atom] = allResidues[atom]
            for residue in tree.getroot().findall('Residue'):
                name = residue.attrib['name']
                for id in residue.attrib:
                    if id == 'name' or id.startswith('alt'):
                        PDBFile._residueNameReplacements[residue.attrib[id]] = name
                if 'type' not in residue.attrib:
                    atoms = copy(allResidues)
                elif residue.attrib['type'] == 'Protein':
                    atoms = copy(proteinResidues)
                elif residue.attrib['type'] == 'Nucleic':
                    atoms = copy(nucleicAcidResidues)
                else:
                    atoms = copy(allResidues)
                PDBFile._parseResidueAtoms(residue, atoms)
                PDBFile._atomNameReplacements[name] = atoms

    @staticmethod
    def _parseResidueAtoms(residue, map):
        for atom in residue.findall('Atom'):
            name = atom.attrib['name']
            for id in atom.attrib:
                map[atom.attrib[id]] = name

    @staticmethod
    def writeFile(topology, positions, file=sys.stdout, modelIndex=None):
        """Write a PDB file containing a single model.

        Parameters:
         - topology (Topology) The Topology defining the model to write
         - positions (list) The list of atomic positions to write
         - file (file=stdout) A file to write to
        """
        PDBFile.writeHeader(topology, file)
        PDBFile.writeModel(topology, positions, file)
        PDBFile.writeFooter(topology, file)

    @staticmethod
    def writeHeader(topology, file=sys.stdout):
        """Write out the header for a PDB file.

        Parameters:
         - topology (Topology) The Topology defining the molecular system being written
         - file (file=stdout) A file to write the file to
        """
        boxSize = topology.getUnitCellDimensions()
        if boxSize is not None:
            print >>file, "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 " % boxSize

    @staticmethod
    def writeModel(topology, positions, file=sys.stdout, modelIndex=None):
        """Write out a model to a PDB file.

        Parameters:
         - topology (Topology) The Topology defining the model to write
         - positions (list) The list of atomic positions to write
         - file (file=stdout) A file to write the model to
         - modelIndex (int=None) If not None, the model will be surrounded by MODEL/ENDMDL records with this index
        """
        if len(list(topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if np.any(np.isnan(positions)):
            raise ValueError('Particle position is NaN')
        if np.any(np.isinf(positions)):
            raise ValueError('Particle position is infinite')
        atomIndex = 1
        posIndex = 0
        if modelIndex is not None:
            print >>file, "MODEL     %4d" % modelIndex
        for (chainIndex, chain) in enumerate(topology.chains()):
            chainName = chr(ord('A')+chainIndex)
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
                    print >>file, "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00" % (atomIndex%100000, atomName, resName, chainName, (resIndex+1)%10000, coords[0], coords[1], coords[2])
                    posIndex += 1
                    atomIndex += 1
                if resIndex == len(residues)-1:
                    print >>file, "TER   %5d      %3s %s%4d" % (atomIndex, resName, chainName, resIndex+1)
                    atomIndex += 1
        if modelIndex is not None:
            print >>file, "ENDMDL"

    @staticmethod
    def writeFooter(topology, file=sys.stdout):
        """Write out the footer for a PDB file.

        Parameters:
         - topology (Topology) The Topology defining the molecular system being written
         - file (file=stdout) A file to write the file to
        """
        print >>file, "END"
