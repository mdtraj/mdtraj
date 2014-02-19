##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Peter Eastman, Robert McGibbon
# Contributors: Kyle A. Beauchamp
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
# toolkit, copyright (c) 2012 Stanford University and Peter Eastman. Those
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

from __future__ import print_function, division
import os
import numpy as np
import itertools
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
            resSeq = getattr(residue, 'resSeq', None) or residue.index
            newResidue = newTopology.add_residue(residue.name, newChain, resSeq)
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

    The structure of a Topology object is similar to that of a PDB file.
    It consists of a set of Chains (often but not always corresponding to
    polymer chains).  Each Chain contains a set of Residues, and each Residue
    contains a set of Atoms.  In addition, the Topology stores a list of which
    atom pairs are bonded to each other.

    Atom and residue names should follow the PDB 3.0 nomenclature for all
    molecules for which one exists.

    Attributes
    ----------
    chains : generator
        Iterator over all Chains in the Topology.
    residues : genetator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.

    Examples
    --------
    >>> topology = md.load('example.pdb').topology
    >>> print(topology)
    <mdtraj.Topology with 1 chains, 3 residues, 22 atoms, 21 bonds at 0x105a98e90>
    >>> table, bonds = topology.to_dataframe()
    >>> print(table.head())
       serial name element  resSeq resName  chainID
    0       0   H1       H       0     CYS        0
    1       1  CH3       C       0     CYS        0
    2       2   H2       H       0     CYS        0
    3       3   H3       H       0     CYS        0
    4       4    C       C       0     CYS        0
    >>> # rename residue "CYS" to "CYSS"
    >>> table[table['residue'] == 'CYS']['residue'] = 'CYSS'
    >>> print(table.head())
       serial name element  resSeq resName   chainID
    0       0   H1       H       0     CYSS        0
    1       1  CH3       C       0     CYSS        0
    2       2   H2       H       0     CYSS        0
    3       3   H3       H       0     CYSS        0
    4       4    C       C       0     CYSS        0
    >>> t2 = md.Topology.from_dataframe(table, bonds)
    """

    _standardBonds = {}

    def __init__(self):
        """Create a new Topology object"""
        self._chains = []
        self._numResidues = 0
        self._numAtoms = 0
        self._bonds = []
        self._atoms = []
        self._residues = []

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "<mdtraj.Topology with %d chains, %d residues, %d atoms, %d bonds>" % (self.n_chains, self.n_residues, self.n_atoms, len(self._bonds))

    def __repr__(self):
        return "<mdtraj.Topology with %d chains, %d residues, %d atoms, %d bonds at 0x%02x>" % (self.n_chains, self.n_residues, self.n_atoms, len(self._bonds), id(self))

    def copy(self):
        """Return a copy of the topology

        Returns
        -------
        out : Topology
            A copy of this topology
        """
        out = Topology()
        for chain in self.chains:
            c = out.add_chain()
            for residue in chain.residues:
                r = out.add_residue(residue.name, c, residue.resSeq)
                for atom in residue.atoms:
                    out.add_atom(atom.name, atom.element, r)

        for a1, a2 in self.bonds:
            out.add_bond(a1, a2)

        return out

    def __copy__(self, *args):
        return self.copy()

    def __deepcopy__(self, *args):
        return self.copy()

    def join(self, other):
        """Join two topologies together

        Parameters
        ----------
        other : Topology
            Another topology object

        Returns
        -------
        out : Topology
            A joint topology, with all of the atoms/residues/chains/bonds
            in each of the individual topologies
        """
        if not isinstance(other, Topology):
            raise ValueError('other must be an instance of Topology to join')
        out = self.copy()

        atom_mapping = {}
        for chain in other.chains:
            c = out.add_chain()
            for residue in chain.residues:
                r = out.add_residue(residue.name, c, residue.resSeq)
                for atom in residue.atoms:
                    a = out.add_atom(atom.name, atom.element, r)
                    atom_mapping[atom] = a

        for a1, a2 in other.bonds:
            out.add_bond(atom_mapping[a1], atom_mapping[a2])

        return out

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
        -------
        atoms : pandas.DataFrame
            The atoms in the topology, represented as a data frame.
        bonds : np.ndarray
            The bonds in this topology, represented as an n_bonds x 2 array
            of the indices of the atoms involved in each bond.
        """
        pd = import_('pandas')
        data = []
        for atom in self.atoms:
            if atom.element is None:
                element_symbol = ""
            else:
                element_symbol = atom.element.symbol
            data.append((atom.index, atom.name, element_symbol,
                         atom.residue.resSeq, atom.residue.name,
                         atom.residue.chain.index))

        atoms = pd.DataFrame(data, columns=["serial", "name", "element",
                                            "resSeq", "resName", "chainID"])
        atoms = atoms.set_index("serial")
        bonds = np.array([(a.index, b.index) for (a, b) in self.bonds])
        return atoms, bonds

    @classmethod
    def from_dataframe(cls, atoms, bonds=None):
        """Create a mdtraj topology from a pandas data frame

        Parameters
        ----------
        atoms : pandas.DataFrame
            The atoms in the topology, represented as a data frame. This data
            frame should have columns "serial" (atom index), "name" (atom name),
            "element" (atom's element), "resSeq" (index of the residue)
            "resName" (name of the residue), "chainID" (index of the chain),
            following the same conventions as wwPDB 3.0 format.
        bonds : np.ndarray, shape=(n_bonds, 2), dtype=int, optional
            The bonds in the topology, represented as an n_bonds x 2 array
            of the indices of the atoms involved in each bond. Specifiying
            bonds here is optional. To create standard protein bonds, you can
            use `create_standard_bonds` to "fill in" the bonds on your newly
            created Topology object

        See Also
        --------
        create_standard_bonds
        """
        pd = import_('pandas')
        from mdtraj import pdb

        for col in ["name", "element", "resSeq" , "resName", "chainID"]:
            if col not in atoms.columns:
                raise ValueError('dataframe must have column %s' % col)

        out = cls()
        if not isinstance(atoms, pd.DataFrame):
            raise TypeError('atoms must be an instance of pandas.DataFrame. '
                            'You supplied a %s' % type(atoms))
        if not isinstance(bonds, np.ndarray):
            raise TypeError('bonds must be an instance of numpy.ndarray. '
                            'You supplied a %s' % type(bonds))

        if not np.all(np.arange(len(atoms)) == atoms.index):
            raise ValueError('atoms must be uniquely numbered starting from zero.')
        out._atoms = [None for i in range(len(atoms))]
        for ci in np.unique(atoms['chainID']):
            chain_atoms = atoms[atoms['chainID'] == ci]
            c = out.add_chain()

            for ri in np.unique(chain_atoms['resSeq']):
                residue_atoms = chain_atoms[chain_atoms['resSeq'] == ri]
                rnames = residue_atoms['resName']
                residue_name = np.array(rnames)[0]
                if not np.all(rnames == residue_name):
                    raise ValueError('All of the atoms with residue index %d do not share the same residue name' % ri)
                r = out.add_residue(residue_name, c, ri)

                for ai, atom in residue_atoms.iterrows():
                    if atom['element'] == "":
                        element = None
                    else:
                        element = pdb.element.get_by_symbol(atom['element'])
                    a = Atom(atom['name'], element, ai, r)
                    out._atoms[ai] = a
                    r._atoms.append(a)

        if bonds is not None:
            for ai1, ai2 in bonds:
                out.add_bond(out.atom(ai1), out.atom(ai2))

        out._numAtoms = out.n_atoms
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
                if (r1.index != r1.index) or (r1.name != r2.name): # or (r1.resSeq != r2.resSeq):
                    return False
                if len(r1._atoms) != len(r2._atoms):
                    return False

                for a1, a2 in zip(r1.atoms, r2.atoms):
                    if (a1.index != a2.index)  or (a1.name != a2.name):
                        return False
                    if a1.element is not None and a2.element is not None:
                        if a1.element != a2.element:
                            return False
                        #for attr in ['atomic_number', 'name', 'symbol']:
                        #    if getattr(a1.element, attr) != getattr(a2.element, attr):
                        #        return False

        if len(self._bonds) != len(other._bonds):
            return False
        for (a1, b1), (a2, b2) in zip(self.bonds, other.bonds):
            if (a1.index != a2.index) or (b1.index != b2.index):
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

    def add_residue(self, name, chain, resSeq=None):
        """Create a new Residue and add it to the Topology.

        Parameters
        ----------
        name : str
            The name of the residue to add
        chain : mdtraj.topology.Chain
            The Chain to add it to
        resSeq : int, optional
            Residue sequence number, such as from a PDB record. These sequence
            numbers are arbitrary, and do not necessarily start at 0 (or 1).
            If not supplied, the resSeq attribute will be set to the
            residue's sequential (0 based) index.

        Returns
        -------
        residue : mdtraj.topology.Residue
            The newly created Residue
        """
        if resSeq is None:
            resSeq = self._numResidues
        residue = Residue(name, self._numResidues, chain, resSeq)
        self._residues.append(residue)
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
        self._atoms.append(atom)
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

    def chain(self, index):
        """Get a specific chain by index.  These indices
        start from zero.

        Returns
        -------
        chain : Chain
            The `index`-th chain in the topology.
        """
        return self._chains[index]

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
    def n_chains(self):
        """Get the number of chains in the Topology"""
        return len(self._chains)

    def residue(self, index):
        """Get a specific residue by index.  These indices
        start from zero.

        Returns
        -------
        residue : Residue
            The `index`-th residue in the topology.
        """
        return self._residues[index]

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
    def n_residues(self):
        """Get the number of residues in the Topology"""
        return len(self._residues)

    def atom(self, index):
        """Get a specific atom by index. These indices
        start from zero.

        Returns
        -------
        atom : Atom
            The `index`-th atom in the topology.
        """
        return self._atoms[index]

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
    def n_atoms(self):
        """Get the number of atoms in the Topology"""
        return len(self._atoms)

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

    def residue(self, index):
        """Get a specific residue in this Chain

        Returns
        -------
        residue : Residue
        """
        return self._residue[index]

    @property
    def n_residues(self):
        "Get the number of residues in this Chain"
        return len(self._residues)

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

    def atom(self, index):
        """Get a specific atom in this Chain

        Returns
        -------
        atom : Atom
        """
        # this could be made faster by caching the list
        # of atoms internally if necessary
        return next(itertools.islice(self.atoms, index, index+1))

    @property
    def n_atoms(self):
        """Get the number of atoms in this Chain"""
        return sum(r.n_atoms for r in self._residues)


class Residue(object):
    """A Residue object represents a residue within a Topology.

    Attributes
    ----------
    name : str
        The name of the Residue
    index : int
        The index of the Residue within its Topology
    chain : int
        The residue sequence number
    """
    def __init__(self, name, index, chain, resSeq):
        """Construct a new Residue.  You should call add_residue()
        on the Topology instead of calling this directly."""
        self.name = name
        self.index = index
        self.chain = chain
        self.resSeq = resSeq
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

    def atom(self, index):
        """Get a specific atom in this Residue.

        Returns
        -------
        atom : Atom
        """
        return self._atoms[index]

    @property
    def n_atoms(self):
        """Get the number of atoms in this Residue"""
        return len(self._atoms)

    def  __str__(self):
        return '%s%s' % (self.name, self.resSeq)

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

    def __eq__(self, other):
        """ Check whether two Atom objects are equal. """
        if self.name != other.name:
            return False
        if self.index != other.index:
            return False
        if self.element.name != other.element.name:
            return False
        if self.residue.name != other.residue.name:
            return False
        if self.residue.index != other.residue.index:
            return False
        if self.residue.chain.index != other.residue.chain.index:
            return False
        return True

    def __hash__(self):
        """ A quick comparison. """
        return self.index

    def __str__(self):
        return '%s-%s' % (self.residue, self.name)
