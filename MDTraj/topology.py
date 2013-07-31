# Copyright 2012 the authors
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
#
# topology.py: Used for storing topological information about a system.
#
# This is part of the OpenMM molecular simulation toolkit originating from
# Simbios, the NIH National Center for Physics-Based Simulation of
# Biological Structures at Stanford, funded under the NIH Roadmap for
# Medical Research, grant U54 GM072970. See https://simtk.org.
#
# Portions copyright (c) 2012 Stanford University and the Authors.
# Authors: Peter Eastman
# Contributors: mdtraj developers
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
# Imports
##############################################################################

import cPickle as pickle
import os
import numpy as np
import xml.etree.ElementTree as etree

from mdtraj.utils import ilen, import_

##############################################################################
# Utilities
##############################################################################


def _topology_from_subset(topology, atom_indices):
    """Create a new topology that only contains the supplied indices

    Note
    ----
    This really should be a copy constructor (class method) on Topology,
    but I want it to work on either the mdtraj topology OR the OpenMM
    topology. An inplace version for the topology object we have here
    is also available.

    Parameters
    ----------
    topology : topology
        The base topology
    atom_indices : list([int])
        The indices of the atoms to keep
    """
    newTopology = Topology()
    old_atom_to_new_atom = {}

    for chain in topology._chains:
        newChain = newTopology.add_chain()
        for residue in chain._residues:
            newResidue = newTopology.add_residue(residue.name, newChain)
            for atom in residue._atoms:
                if atom.index in atom_indices:
                    newAtom = newTopology.add_atom(atom.name, atom.element, newResidue)
                    old_atom_to_new_atom[atom] = newAtom

    bondsiter = topology.bonds
    if not hasattr(bondsiter, '__iter__'):
        bondsiter = bondsiter()

    for atom1, atom2 in bondsiter:
        try:
            newTopology.add_bond(old_atom_to_new_atom[atom1],
                                old_atom_to_new_atom[atom2])
        except KeyError:
            pass
            # we only put bonds into the new topology if both of their partners
            # were indexed and thus HAVE a new atom

    # Delete empty residues
    for chain in newTopology._chains:
        chain._residues = [r for r in chain._residues if len(r._atoms) > 0]
    # Delete empty chains
    newTopology._chains = [c for c in newTopology._chains if len(c._residues) > 0]
    # Re-set the numAtoms and numResidues
    newTopology._numAtoms = ilen(newTopology.atoms)
    newTopology._numResidues = ilen(newTopology.residues)

    return newTopology


##############################################################################
# Classes
##############################################################################


class Topology(object):
    """Topology stores the topological information about a system.

    The structure of a Topology object is similar to that of a PDB file.  It consists of a set of Chains
    (often but not always corresponding to polymer chains).  Each Chain contains a set of Residues,
    and each Residue contains a set of Atoms.  In addition, the Topology stores a list of which atom
    pairs are bonded to each other, and the dimensions of the crystallographic unit cell.

    Atom and residue names should follow the PDB 3.0 nomenclature for all molecules for which one exists.

    Attributes
    ----------
    chains : generator
        Iterator over all Chains in the Topology.
    residues : genetator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.
    """

    _standardBonds = {}

    def __init__(self):
        """Create a new Topology object"""
        self._chains = []
        self._numResidues = 0
        self._numAtoms = 0
        self._bonds = []

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_openmm(self):
        """Convert this topology into OpenMM topology

        Returns
        -------
        topology : simtk.openmm.app.Topology
           This topology, as an OpenMM topology
        """
        app = import_('simtk.openmm.app')

        out = app.Topology()
        atom_mapping = {}

        for chain in self.chains:
            c = out.addChain()
            for residue in chain.residues:
                r = out.addResidue(residue.name, c)
                for atom in residue.atoms:
                    a = out.addAtom(atom.name, app.Element.getBySymbol(atom.element.symbol), r)
                    atom_mapping[atom] = a

        for a1, a2 in self.bonds:
            out.addBond(atom_mapping[a1], atom_mapping[a2])

        return out


    @classmethod
    def from_openmm(cls, value):
        """Create a mdtraj topology from an OpenMM topology

        Parameters
        ----------
        value : simtk.openmm.app.Topology
            An OpenMM topology that you wish to convert to a
            mdtraj topology.
        """
        from mdtraj import pdb
        app = import_('simtk.openmm.app')

        if not isinstance(value, app.Topology):
            raise TypeError('value must be an OpenMM Topology. '
                            'You supplied a %s' % type(value))

        out = cls()
        atom_mapping = {}

        for chain in value.chains():
            c = out.add_chain()
            for residue in chain.residues():
                r = out.add_residue(residue.name, c)
                for atom in residue.atoms():
                    a = out.add_atom(atom.name, pdb.element.get_by_symbol(atom.element.symbol), r)
                    atom_mapping[atom] = a

        for a1, a2 in value.bonds():
            out.add_bond(atom_mapping[a1], atom_mapping[a2])

        return out


    def to_dataframe(self):
        """Convert this topology into a pandas dataframe

        Returns
        ------
        df : pandas.DataFrame
            This topology, represented as a data frame.
        """
        pd = import_('pandas')
        data = []
        for atom in self.atoms:
            data.append((atom.index, atom.name, atom.element.symbol,
                         atom.residue.index, atom.residue.name,
                         atom.residue.chain.index))

        out = pd.DataFrame(data, columns = ["index", "atom", "element",
                                            "rindex" , "residue", "chain"])
        return out
        

    @classmethod
    def from_dataframe(cls, value):
        """Create a mdtraj topology from a pandas data frame

        Parameters
        ----------
        value : pandas.DataFrame
            A pandas dataframe topology that you wish to convert to a
            mdtraj topology.
        """
        pd = import_('pandas')
        out = cls()
        if not isinstance(value, pd.DataFrame):
            raise TypeError('value must be an instance of pandas.DataFrame. '
                            'You supplied a %s' % type(value))

        for ci in np.unique(value['chain']):
            chain_atoms = value[value['chain'] == ci]
            c = out.add_chain()

            for ri in np.unique(chain_atoms['rindex']):
                rnames = chain_atoms[chain_atoms['rindex'] == ri]['residue']
                print len(rnames)
                print rnames
                print np.array(rnames)[0]

                
                if not np.all(rnames == np.array(rnames)[0]):
                    raise ValueError('All of the atoms with residue index %d do not share the same residue name' % ri)
                r = out.add_residue(rnames[0], c)



        import IPython as ip; ip.embed()
        #for i, row in value.iterrows():
        #    c = out._chains[row['chain']]
            


        return out

    def __eq__(self, other):
        """Are two topologies equal?

        Parameters
        ----------
        other : object
            The object to compare to

        Returns
        -------
        equality : bool
            Are the two topologies identical?
        """
        if not isinstance(other, Topology):
            return False
        if self is other:
            return True

        if len(self._chains) != len(other._chains):
            return False

        for c1, c2 in zip(self.chains, other.chains):
            if c1.index != c2.index:
                return False
            if len(c1._residues) != len(c2._residues):
                return False

            for r1, r2 in zip(c1.residues, c2.residues):
                if (r1.index != r1.index) or (r1.name != r2.name):
                    return False
                if len(r1._atoms) != len(r2._atoms):
                    return False

                for a1, a2 in zip(r1.atoms, r2.atoms):
                    if (a1.index != a2.index)  or (a1.name != a2.name):
                        return False
                    for attr in ['atomic_number', 'name', 'symbol']:
                        if getattr(a1.element, attr) != getattr(a2.element, attr):
                            return False
        return True


    def add_chain(self):
        """Create a new Chain and add it to the Topology.

        Returns
        -------
        chain : mdtraj.topology.Chain
            the newly created Chain
        """
        chain = Chain(len(self._chains), self)
        self._chains.append(chain)
        return chain

    def add_residue(self, name, chain):
        """Create a new Residue and add it to the Topology.

        Parameters
        ----------
        name : str
            The name of the residue to add
        chain : mdtraj.topology.Chain
            The Chain to add it to

        Returns
        -------
        residue : mdtraj.topology.Residue
            The newly created Residue
        """
        residue = Residue(name, self._numResidues, chain)
        self._numResidues += 1
        chain._residues.append(residue)
        return residue

    def add_atom(self, name, element, residue):
        """Create a new Atom and add it to the Topology.

        Parameters
        ----------
        name : str
            The name of the atom to add
        element : mdtraj.pdb.element.Element
            The element of the atom to add
        residue : mdtraj.topology.Residue
            The Residue to add it to

        Returns
        -------
        atom : mdtraj.topology.Atom
            the newly created Atom
        """
        atom = Atom(name, element, self._numAtoms, residue)
        self._numAtoms += 1
        residue._atoms.append(atom)
        return atom

    def add_bond(self, atom1, atom2):
        """Create a new bond and add it to the Topology.

        Parameters
        ----------
        atom1 : mdtraj.topology.Atom
            The first Atom connected by the bond
        atom2 : mdtraj.topology.Atom
            The second Atom connected by the bond
        """
        self._bonds.append((atom1, atom2))

    @property
    def chains(self):
        """Iterator over all Chains in the Topology.

        Returns
        -------
        chainiter : listiterator
            Iterator over all Chains in the Topology.
        """
        return iter(self._chains)

    @property
    def residues(self):
        """Iterator over all Residues in the Topology.

        Returns
        -------
        residueiter : generator
            Iterator over all Residues in the Topology.
        """
        for chain in self._chains:
            for residue in chain._residues:
                yield residue

    @property
    def atoms(self):
        """Iterator over all Atoms in the Topology.

        Returns
        -------
        atomiter : generator
            Iterator over all Atoms in the Topology.
        """
        for chain in self._chains:
            for residue in chain._residues:
                for atom in residue._atoms:
                    yield atom

    @property
    def bonds(self):
        """Iterator over all bonds (each represented as a tuple of two Atoms) in the Topology.

        Returns
        -------
        atomiter : generator
            Iterator over all tuple of Atoms in the Trajectory involved in a bond.
        """
        return iter(self._bonds)


    def create_standard_bonds(self):
        """Create bonds based on the atom and residue names for all standard residue types.
        """
        if len(Topology._standardBonds) == 0:
            # Load the standard bond defitions.

            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'pdb', 'data', 'residues.xml'))
            for residue in tree.getroot().findall('Residue'):
                bonds = []
                Topology._standardBonds[residue.attrib['name']] = bonds
                for bond in residue.findall('Bond'):
                    bonds.append((bond.attrib['from'], bond.attrib['to']))
        for chain in self._chains:
            # First build a map of atom names to atoms.

            atomMaps = []
            for residue in chain._residues:
                atomMap = {}
                atomMaps.append(atomMap)
                for atom in residue._atoms:
                    atomMap[atom.name] = atom

            # Loop over residues and construct bonds.

            for i in range(len(chain._residues)):
                name = chain._residues[i].name
                if name in Topology._standardBonds:
                    for bond in Topology._standardBonds[name]:
                        if bond[0].startswith('-') and i > 0:
                            fromResidue = i-1
                            fromAtom = bond[0][1:]
                        elif bond[0].startswith('+') and i <len(chain._residues):
                            fromResidue = i+1
                            fromAtom = bond[0][1:]
                        else:
                            fromResidue = i
                            fromAtom = bond[0]
                        if bond[1].startswith('-') and i > 0:
                            toResidue = i-1
                            toAtom = bond[1][1:]
                        elif bond[1].startswith('+') and i <len(chain._residues):
                            toResidue = i+1
                            toAtom = bond[1][1:]
                        else:
                            toResidue = i
                            toAtom = bond[1]
                        if fromAtom in atomMaps[fromResidue] and toAtom in atomMaps[toResidue]:
                            self.add_bond(atomMaps[fromResidue][fromAtom], atomMaps[toResidue][toAtom])

    def create_disulfide_bonds(self, positions):
        """Identify disulfide bonds based on proximity and add them to the Topology.

        Parameters
        ----------
        positions : list
            The list of atomic positions based on which to identify bonded atoms
        """
        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names

        cyx = [res for res in self.residues if res.name == 'CYS' and isCyx(res)]
        atomNames = [[atom.name for atom in res._atoms] for res in cyx]
        for i in range(len(cyx)):
            sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
            pos1 = positions[sg1.index]
            for j in range(i):
                sg2 = cyx[j]._atoms[atomNames[j].index('SG')]
                pos2 = positions[sg2.index]
                delta = [x-y for (x,y) in zip(pos1, pos2)]
                distance = np.sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
                if distance < 0.3: # this is supposed to be nm. I think we're good
                    self.add_bond(sg1, sg2)

    def subset(self, atom_indices):
        """Create a new Topology from a subset of the atoms in an existing topology.

        Notes
        -----
        The existing topology will not be altered.

        Parameters
        ----------
        atom_indices array_like
            A list of the indices corresponding to the atoms in that you'd
            like to retain.
        """
        return _topology_from_subset(self, atom_indices)


class Chain(object):
    """A Chain object represents a chain within a Topology.

    Attributes
    ----------
    index : int
        The index of the Chain within its Topology
    topology : mdtraj.Topology
        The Topology this Chain belongs to
    residues : genetator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.
    """
    def __init__(self, index, topology):
        """Construct a new Chain.  You should call add_chain() on the Topology instead of calling this directly."""
        ## The index of the Chain within its Topology
        self.index = index
        ## The Topology this Chain belongs to
        self.topology = topology
        self._residues = []

    @property
    def residues(self):
        """Iterator over all Residues in the Chain.

        Returns
        -------
        residueiter : listiterator
            Iterator over all Residues in the Topology.
        """
        return iter(self._residues)

    @property
    def atoms(self):
        """Iterator over all Atoms in the Chain.

        Returns
        -------
        atomiter : generator
            Iterator over all Atoms in the Chain.
        """
        for residue in self._residues:
            for atom in residue._atoms:
                yield atom

class Residue(object):
    """A Residue object represents a residue within a Topology.

    Attributes
    ----------
    atoms : genetator
    """
    def __init__(self, name, index, chain):
        """Construct a new Residue.  You should call add_residue()
        on the Topology instead of calling this directly."""
        ## The name of the Residue
        self.name = name
        ## The index of the Residue within its Topology
        self.index = index
        ## The Chain this Residue belongs to
        self.chain = chain
        self._atoms = []

    @property
    def atoms(self):
        """Iterator over all Atoms in the Residue.

        Returns
        -------
        atomiter : listiterator
            Iterator over all Atoms in the Residue.
        """
        return iter(self._atoms)

class Atom(object):
    """An Atom object represents a residue within a Topology.

    Attributes
    ----------
    name : str
        The name of the Atom
    element : mdtraj.pdb.element.Element
        The element of the Atoms
    index : int
        The index of the Atom within its Topology
    residue : mdtraj.topology.Residue
        The Residue this Atom belongs to
    """

    def __init__(self, name, element, index, residue):
        """Construct a new Atom.  You should call add_atom() on the Topology instead of calling this directly."""
        ## The name of the Atom
        self.name = name
        ## That Atom's element
        self.element = element
        ## The index of the Atom within its Topology
        self.index = index
        ## The Residue this Atom belongs to
        self.residue = residue
