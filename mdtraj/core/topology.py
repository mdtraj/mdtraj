##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Peter Eastman, Robert McGibbon
# Contributors: Kyle A. Beauchamp, Matthew Harrigan, Carlos Xavier Hernandez
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

import itertools
import numpy as np
import os
import xml.etree.ElementTree as etree
from collections import namedtuple
import warnings

from mdtraj.core import element as elem
from mdtraj.core.residue_names import (_PROTEIN_RESIDUES, _WATER_RESIDUES,
                                       _AMINO_ACID_CODES)
from mdtraj.core.selection import parse_selection
from mdtraj.utils import ilen, import_, ensure_type
from mdtraj.utils.six import string_types
from mdtraj.utils.singleton import Singleton

##############################################################################
# Utilities
##############################################################################


def _topology_from_subset(topology, atom_indices):
    """Create a new topology that only contains the supplied indices

    Note
    ----
    This really should be a copy constructor (class method) on Topology.
    It used to work on OpenMM topologies, but we've diverged where that no
    longer works.

    Parameters
    ----------
    topology : mdtraj.Topology
        The base topology
    atom_indices : array_like, dtype=int
        The indices of the atoms to keep
    """
    np_indices = np.array(atom_indices)
    if len(np_indices) > 1:
        if (np_indices[1:] - np_indices[:-1]).min() < 1:
            warnings.warn('atom_indices are not monotonically increasing')
        if len(np.unique(np_indices)) < len(np_indices):
            warnings.warn('atom_indices are not unique')

    newTopology = Topology()
    old_atom_to_new_atom = {}

    for chain in topology._chains:
        newChain = newTopology.add_chain()
        for residue in chain._residues:
            resSeq = getattr(residue, 'resSeq', None) or residue.index
            newResidue = None
            for atom in residue._atoms:
                if atom.index in atom_indices:
                    try:  # OpenMM Topology objects don't have serial attributes, so we have to check first.
                        serial = atom.serial
                    except AttributeError:
                        serial = None
                    if newResidue is None:
                        newResidue = newTopology.add_residue(residue.name, newChain,
                                                             resSeq, residue.segment_id)
                    newAtom = newTopology.add_atom(atom.name, atom.element,
                                                   newResidue, serial=serial)
                    old_atom_to_new_atom[atom] = newAtom

    bondsiter = topology.bonds
    if not hasattr(bondsiter, '__iter__'):
        bondsiter = bondsiter()

    for bond in bondsiter:
        try:
            atom1, atom2 = bond
            newTopology.add_bond(old_atom_to_new_atom[atom1],
                                 old_atom_to_new_atom[atom2],
                                 type=bond.type,
                                 order=bond.order)
        except KeyError:
            pass
            # we only put bonds into the new topology if both of their partners
            # were indexed and thus HAVE a new atom

    # Delete empty residues
    newTopology._residues = [r for r in newTopology._residues if len(r._atoms) > 0]
    for chain in newTopology._chains:
        chain._residues = [r for r in chain._residues if len(r._atoms) > 0]

    # Delete empty chains
    newTopology._chains = [c for c in newTopology._chains
                           if len(c._residues) > 0]
    # Re-set the numAtoms and numResidues
    newTopology._numAtoms = ilen(newTopology.atoms)
    newTopology._numResidues = ilen(newTopology.residues)

    # Reset the chain indices
    for i, chain in enumerate(newTopology.chains):
        chain.index = i
    # Reset the residue indices
    for i, res in enumerate(newTopology.residues):
        res.index = i

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
    residues : generator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.
    bonds : generator
        Iterator over all Bonds in the Topology

    Examples
    --------
    >>> topology = md.load('example.pdb').topology
    >>> print(topology)
    <mdtraj.Topology with 1 chains, 3 residues, 22 atoms, 21 bonds at 0x105a98e90>
    >>> table, bonds = topology.to_dataframe()
    >>> print(table.head())
       serial name element  resSeq resName  chainID
    0       0   H1       H       1     CYS        0
    1       1  CH3       C       1     CYS        0
    2       2   H2       H       1     CYS        0
    3       3   H3       H       1     CYS        0
    4       4    C       C       1     CYS        0
    >>> # rename residue "CYS" to "CYSS"
    >>> table[table['residue'] == 'CYS']['residue'] = 'CYSS'
    >>> print(table.head())
       serial name element  resSeq resName   chainID
    0       0   H1       H       1     CYSS        0
    1       1  CH3       C       1     CYSS        0
    2       2   H2       H       1     CYSS        0
    3       3   H3       H       1     CYSS        0
    4       4    C       C       1     CYSS        0
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
        return "<%s>" % (self._string_summary_basic())

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def _string_summary_basic(self):
        return ("mdtraj.Topology with %d chains, %d residues, "
                "%d atoms, %d bonds" % (self.n_chains, self.n_residues,
                                        self.n_atoms, len(self._bonds)))

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
                r = out.add_residue(residue.name, c, residue.resSeq, residue.segment_id)
                for atom in residue.atoms:
                    out.add_atom(atom.name, atom.element, r, serial=atom.serial)

        for bond in self.bonds:
            a1, a2 = bond
            out.add_bond(a1, a2, type=bond.type, order=bond.order)

        return out

    def __copy__(self, *args):
        return self.copy()

    def __deepcopy__(self, *args):
        return self.copy()

    def __hash__(self):
        hash_value = hash(tuple(self._chains))
        hash_value ^= hash(tuple(self._atoms))
        hash_value ^= hash(tuple(self._bonds))
        hash_value ^= hash(tuple(self._residues))

        return hash_value

    def join(self, other, keep_resSeq=True):
        """Join two topologies together

        Parameters
        ----------
        other : Topology
            Another topology object
        keep_resSeq : bool, optional, default=True
            if False the residue numbers (resSeq) of the
            topology that is joined are updated in order to
            continue from the last resSeq of this topology

        Returns
        -------
        out : Topology
            A joint topology, with all of the atoms/residues/chains/bonds
            in each of the individual topologies
        """
        if not isinstance(other, Topology):
            raise ValueError('other must be an instance of Topology to join')
        out = self.copy()

        #I need this in order to have the resSeq of the
        #new residues to continue from the one of the
        #las residue of this topology
        if not keep_resSeq:
            out_resSeq = out.atom(-1).residue.resSeq

        atom_mapping = {}
        for chain in other.chains:
            c = out.add_chain()
            for residue in chain.residues:
                if keep_resSeq:
                    out_resSeq = residue.resSeq
                else:
                    out_resSeq += 1
                r = out.add_residue(residue.name, c, out_resSeq, residue.segment_id)
                for atom in residue.atoms:
                    a = out.add_atom(atom.name, atom.element, r,
                                     serial=atom.serial)
                    atom_mapping[atom] = a

        for bond in other.bonds:
            a1, a2 = bond
            out.add_bond(atom_mapping[a1], atom_mapping[a2], type=bond.type, order=bond.order)

        return out

    def to_fasta(self, chain=None):
        """Convert this topology into FASTA string

        Parameters
        ----------
        chain : Integer, optional, default=None
            If specified, will return the FASTA string for this chain in the
            Topology.

        Returns
        -------
        fasta : String or list of Strings
           A FASTA string for each chain specified.
        """
        fasta = lambda c: "".join([res.code for res in c.residues
                                   if res.is_protein and res.code is not None])
        if chain is not None:
            if not isinstance(chain, int):
                raise ValueError('chain must be an Integer.')
            return fasta(self._chains[chain])
        else:
            return [fasta(c) for c in self._chains]

    def to_openmm(self, traj=None):
        """Convert this topology into OpenMM topology

        Parameters
        ----------
        traj : MDTraj.Trajectory, optional, default=None
            If specified, use the first frame from this trajectory to
            set the unitcell information in the openmm topology.

        Returns
        -------
        topology : openmm.app.Topology
           This topology, as an OpenMM topology
        """
        app = import_('openmm.app')
        mm = import_('openmm')
        u = import_('openmm.unit')

        out = app.Topology()
        atom_mapping = {}
        bond_mapping = {Single: app.Single,
                        Double: app.Double,
                        Triple: app.Triple,
                        Amide: app.Amide,
                        Aromatic: app.Aromatic,
                        None: None}

        for chain in self.chains:
            c = out.addChain()
            for residue in chain.residues:
                r = out.addResidue(residue.name, c, id=str(residue.resSeq))
                for atom in residue.atoms:
                    if atom.element is elem.virtual:
                        element = None
                    else:
                        element = app.Element.getBySymbol(atom.element.symbol)
                    a = out.addAtom(atom.name, element, r)
                    atom_mapping[atom] = a

        for bond in self.bonds:
            a1, a2 = bond
            out.addBond(atom_mapping[a1], atom_mapping[a2], type=bond_mapping[bond.type], order=bond.order)

        if traj is not None:
            angles = traj.unitcell_angles[0]

            if np.linalg.norm(angles - 90.0) > 1E-4:
                raise(ValueError("Unitcell angles must be 90.0 to use "
                                 "in OpenMM topology."))

            box_vectors = mm.Vec3(*traj.unitcell_lengths[0]) * u.nanometer
            out.setUnitCellDimensions(box_vectors)

        return out

    @classmethod
    def from_openmm(cls, value):
        """Create a mdtraj topology from an OpenMM topology

        Parameters
        ----------
        value : openmm.app.Topology
            An OpenMM topology that you wish to convert to a
            mdtraj topology.
        """
        app = import_('openmm.app')
        bond_mapping = {app.Single: Single,
                        app.Double: Double,
                        app.Triple: Triple,
                        app.Amide: Amide,
                        app.Aromatic: Aromatic,
                        None: None}

        if not isinstance(value, app.Topology):
            raise TypeError('value must be an OpenMM Topology. '
                            'You supplied a %s' % type(value))

        out = cls()
        atom_mapping = {}

        for chain in value.chains():
            c = out.add_chain()
            for residue in chain.residues():
                try:
                    r = out.add_residue(residue.name, c, resSeq=int(residue.id))
                except ValueError:
                    r = out.add_residue(residue.name, c)
                for atom in residue.atoms():
                    if atom.element is None:
                        element = elem.virtual
                    else:
                        element = elem.get_by_symbol(atom.element.symbol)
                    a = out.add_atom(atom.name, element, r)
                    atom_mapping[atom] = a

        for bond in value.bonds():
            a1, a2 = bond
            out.add_bond(atom_mapping[a1], atom_mapping[a2], type=bond_mapping[bond.type], order=bond.order)

        return out

    def to_dataframe(self):
        """Convert this topology into a pandas dataframe

        Returns
        -------
        atoms : pandas.DataFrame
            The atoms in the topology, represented as a data frame.
        bonds : np.ndarray, shape=(n_bonds, 4), dtype=float, Optional
            The bonds in this topology, represented as an n_bonds x 4 array indicating the two atom indices of the bond,
            the bond type, and bond order, cast to floats as the type is mapped from the classes to fractional.
            The atom indices and order are integers cast to float
        """
        pd = import_('pandas')
        data = [(atom.serial, atom.name, atom.element.symbol,
                 atom.residue.resSeq, atom.residue.name,
                 atom.residue.chain.index,atom.segment_id) for atom in self.atoms]

        atoms = pd.DataFrame(data, columns=["serial", "name", "element",
                                            "resSeq", "resName", "chainID","segmentID"])

        bonds = np.zeros([len(self._bonds), 4], dtype=float)
        for index, bond in enumerate(self.bonds):
            if bond.order is None:
                order = 0.0
            else:
                order = bond.order
            try:
                bond_type = float(bond.type)
            except TypeError:
                # Trap the None case
                bond_type = 0.0
            bonds[index] = bond.atom1.index, bond.atom2.index, bond_type, order
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
            and optionally "segmentID", following the same conventions
            as wwPDB 3.0 format.
        bonds : np.ndarray, shape=(n_bonds, 4) or (n_bonds, 2), dtype=float, Optional
            The bonds in the topology, represented as a n_bonds x 4 or n_bonds x 2
            size array of the indices of the atoms involved, type, and order of each
            bond, represented as floats. Type and order are optional. Specifying bonds
            here is optional. To create standard protein bonds, you can use
            `create_standard_bonds` to "fill in" the bonds on your newly created
            Topology object, although type and order of bond are not computed if
            that method is used.

        See Also
        --------
        create_standard_bonds
        """
        pd = import_('pandas')

        if bonds is None:
            bonds = np.zeros([0, 4], dtype=float)

        for col in ["name", "element", "resSeq",
                    "resName", "chainID", "serial"]:
            if col not in atoms.columns:
                raise ValueError('dataframe must have column %s' % col)

        if "segmentID" not in atoms.columns:
            atoms["segmentID"] = ""

        out = cls()
        if not isinstance(atoms, pd.DataFrame):
            raise TypeError('atoms must be an instance of pandas.DataFrame. '
                            'You supplied a %s' % type(atoms))
        if not isinstance(bonds, np.ndarray):
            raise TypeError('bonds must be an instance of numpy.ndarray. '
                            'You supplied a %s' % type(bonds))

        if not np.all(np.arange(len(atoms)) == atoms.index):
            raise ValueError('atoms must be uniquely numbered '
                             'starting from zero.')
        out._atoms = [None for i in range(len(atoms))]

        c = None
        r = None
        previous_chainID = None
        previous_resName = None
        previous_resSeq = None
        for atom_index, atom in atoms.iterrows():

            int(atom_index) # Fixes bizarre hashing issue on Py3K.  See #545

            if atom['chainID'] != previous_chainID:
                previous_chainID = atom['chainID']

                c = out.add_chain()

            if atom['resSeq'] != previous_resSeq or atom['resName'] != previous_resName or c.n_atoms==0:
                previous_resSeq = atom['resSeq']
                previous_resName = atom['resName']

                r = out.add_residue(atom['resName'], c, atom['resSeq'], atom['segmentID'])

            a = Atom(atom['name'], elem.get_by_symbol(atom['element']),
                                atom_index, r, serial=atom['serial'])
            out._atoms[atom_index] = a
            r._atoms.append(a)

        for bond in bonds:
            ai1 = int(bond[0])
            ai2 = int(bond[1])
            try:
                bond_type = float_to_bond_type(bond[2])
                bond_order = int(bond[3])
                if bond_order == 0:
                    bond_order = None
            except IndexError:
                # Does not exist
                bond_type = None
                bond_order = None
            out.add_bond(out.atom(ai1), out.atom(ai2), bond_type, bond_order)

        out._numAtoms = out.n_atoms
        return out

    def to_bondgraph(self):
        """Create a NetworkX graph from the atoms and bonds in this topology

        Returns
        -------
        g : nx.Graph
            A graph whose nodes are the Atoms in this topology, and
            whose edges are the bonds

        See Also
        --------
        atoms
        bonds

        Notes
        -----
        This method requires the NetworkX python package.
        """
        nx = import_('networkx')
        g = nx.Graph()
        g.add_nodes_from(self.atoms)
        g.add_edges_from(self.bonds)
        return g

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
                if (r1.index != r1.index) or (r1.name != r2.name):  # or (r1.resSeq != r2.resSeq):
                    return False
                if len(r1._atoms) != len(r2._atoms):
                    return False

                for a1, a2 in zip(r1.atoms, r2.atoms):
                    if (a1.index != a2.index) or (a1.name != a2.name):
                        return False
                    if a1.element is not a2.element:
                        return False
                    # for attr in ['atomic_number', 'name', 'symbol']:
                    #     if getattr(a1.element, attr) != getattr(a2.element, attr):
                    #         return False

        if len(self._bonds) != len(other._bonds):
            return False

        # the bond ordering is somewhat ambiguous, so try and fix it for comparison
        self_sorted_bonds = sorted([bond for bond in self.bonds])
        other_sorted_bonds = sorted([bond for bond in other.bonds])

        for i, bond in enumerate(self_sorted_bonds):
            if bond != other_sorted_bonds[i]:
                return False

        return True

    def add_chain(self, chain_id=None):
        """Create a new Chain and add it to the Topology.

        Parameters
        ----------
        chain_id : str
            The PDB chainID of the Chain.

        Returns
        -------
        chain : mdtraj.topology.Chain
            the newly created Chain
        """
        chain = Chain(len(self._chains), self, chain_id)
        self._chains.append(chain)
        return chain

    def add_residue(self, name, chain, resSeq=None, segment_id=""):
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
        segment_id : str, optional
            A label for the segment to which this residue belongs

        Returns
        -------
        residue : mdtraj.topology.Residue
            The newly created Residue
        """
        if resSeq is None:
            resSeq = self._numResidues
        residue = Residue(name, self._numResidues, chain, resSeq, segment_id)
        self._residues.append(residue)
        self._numResidues += 1
        chain._residues.append(residue)
        return residue

    def insert_atom(self, name, element, residue, index=None, rindex=None, serial=None):
        """Create a new Atom and insert it into the Topology at a specific position.

        Parameters
        ----------
        name : str
            The name of the atom to add
        element : mdtraj.element.Element
            The element of the atom to add
        residue : mdtraj.topology.Residue
            The Residue to add it to
        index : int
            If provided, the desired index for this atom
            within the topology. Existing atoms with
            indices >= index will be pushed back.
        rindex : int
            If provided, the desired position for this atom
            within the residue
        serial : int
            Serial number associated with the atom.
            This has nothing to do with the actual ordering
            and is solely for PDB labeling purposes.

        Returns
        -------
        atom : mdtraj.topology.Atom
            the newly created Atom
        """
        if element is None:
            element = elem.virtual
        if index is None:
            atom = Atom(name, element, self._numAtoms, residue, serial=serial)
            self._atoms.append(atom)
        else:
            atom = Atom(name, element, index, residue, serial=serial)
            for i in range(index, len(self._atoms)):
                self._atoms[i].index += 1
            self._atoms.insert(index, atom)
        self._numAtoms += 1
        if rindex is None:
            residue._atoms.append(atom)
        else:
            residue._atoms.insert(rindex, atom)
        return atom

    def delete_atom_by_index(self, index):
        """Delete an Atom from the topology.

        Parameters
        ----------
        index : int
            The index of the atom to be removed.
        """
        a = self._atoms[index]
        if a.index != index:
            raise RuntimeError("Index of selected atom does not match order in topology.")
        for i in range(index+1, len(self._atoms)):
            self._atoms[i].index -= 1
        a.residue._atoms.remove(a)
        self._atoms.remove(a)
        self._numAtoms -= 1

    def add_atom(self, name, element, residue, serial=None):
        """Create a new Atom and add it to the Topology.

        Parameters
        ----------
        name : str
            The name of the atom to add
        element : mdtraj.element.Element
            The element of the atom to add
        residue : mdtraj.topology.Residue
            The Residue to add it to
        serial : int
            Serial number associated with the atom.

        Returns
        -------
        atom : mdtraj.topology.Atom
            the newly created Atom
        """
        if element is None:
            element = elem.virtual
        atom = Atom(name, element, self._numAtoms, residue, serial=serial)
        self._atoms.append(atom)
        self._numAtoms += 1
        residue._atoms.append(atom)
        return atom

    def add_bond(self, atom1, atom2, type=None, order=None):
        """Create a new bond and add it to the Topology.

        Parameters
        ----------
        atom1 : mdtraj.topology.Atom
            The first Atom connected by the bond
        atom2 : mdtraj.topology.Atom
            The second Atom connected by the bond
        type : mdtraj.topology.Singleton or None, Default: None, Optional
            Bond type of the bond, or None if not known/provided
        order : 1, 2, 3 or None, Default: None, Optional
            Characteristic order of the bond if known, defaults None
        """
        if atom1.index < atom2.index:
            self._bonds.append(Bond(atom1, atom2, type=type, order=order))
        else:
            self._bonds.append(Bond(atom2, atom1, type=type, order=order))

    def chain(self, index):
        """Get a specific chain by index.  These indices
        start from zero.

        Parameters
        ----------
        index : int
            The index of the chain to select.

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

        Parameters
        ----------
        index : int
            The index of the residue to select.

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
        """Get the number of residues in the Topology. """
        return len(self._residues)

    def atom(self, index):
        """Get a specific atom by index. These indices
        start from zero.

        Parameters
        ----------
        index : int
            The index of the atom to select.

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

    def atoms_by_name(self, name):
        """Iterator over all Atoms in the Topology with a specified name

        Parameters
        ----------
        name : str
            The particular atom name of interest.

        Examples
        --------
        >>> for atom in topology.atoms_by_name('CA'):
        ...     print(atom)

        Returns
        -------
        atomiter : generator
        """
        for atom in self.atoms:
            if atom.name == name:
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
        bonditer : generator
            Iterator over all bonds between Atoms in the Trajectory.
            Each bond can then be iterated to get the atoms.
            I.e.: `for atom1, atom2 in Topology.bonds` yields iterator over pairs of Atoms
                  `for bond in Topology.bonds` yields iterator over bonds
        """
        return iter(self._bonds)

    @property
    def n_bonds(self):
        """Get the number of bonds in the Topology"""
        return len(self._bonds)

    def create_standard_bonds(self):
        """Create bonds based on the atom and residue names for all standard residue types.
        """
        if len(Topology._standardBonds) == 0:
            # Load the standard bond defitions.

            tree = etree.parse(os.path.join(os.path.dirname(__file__), '..',
                               'formats', 'pdb', 'data', 'residues.xml'))
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
                        elif (bond[0].startswith('+')
                              and i < len(chain._residues)):
                            fromResidue = i+1
                            fromAtom = bond[0][1:]
                        else:
                            fromResidue = i
                            fromAtom = bond[0]
                        if bond[1].startswith('-') and i > 0:
                            toResidue = i-1
                            toAtom = bond[1][1:]
                        elif (bond[1].startswith('+')
                              and i < len(chain._residues)):
                            toResidue = i+1
                            toAtom = bond[1][1:]
                        else:
                            toResidue = i
                            toAtom = bond[1]
                        if (fromAtom in atomMaps[fromResidue]
                                and toAtom in atomMaps[toResidue]):
                            self.add_bond(atomMaps[fromResidue][fromAtom],
                                          atomMaps[toResidue][toAtom])

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

        cyx = [res for res in self.residues
               if res.name == 'CYS' and isCyx(res)]
        atomNames = [[atom.name for atom in res._atoms] for res in cyx]
        for i in range(len(cyx)):
            sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
            pos1 = positions[sg1.index]
            for j in range(i):
                sg2 = cyx[j]._atoms[atomNames[j].index('SG')]
                pos2 = positions[sg2.index]
                delta = [x-y for (x, y) in zip(pos1, pos2)]
                distance = np.sqrt(
                    delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
                if distance < 0.3:  # this is supposed to be nm. I think we're good
                    self.add_bond(sg1, sg2)

    def subset(self, atom_indices):
        """Create a new Topology from a subset of the atoms in an existing topology.

        Notes
        -----
        The existing topology will not be altered.

        Parameters
        ----------
        atom_indices : array_like
            A list of the indices corresponding to the atoms in that you'd
            like to retain.
        """
        return _topology_from_subset(self, atom_indices)

    def select_expression(self, selection_string):
        """Translate a atom selection expression into a pure python expression.

        Parameters
        ----------
        selection_string : str
            An expression in the MDTraj atom selection DSL

        Examples
        --------
        >>> topology.select_expression('name O and water')
        "[atom.index for atom in topology.atoms if ((atom.name == 'O') and atom.residue.is_water)]")

        Returns
        -------
        python_string : str
            A string containing a pure python expression, equivalent to the
            selection expression.
        """
        condition = parse_selection(selection_string).source
        fmt_string = "[atom.index for atom in topology.atoms if {condition}]"
        return fmt_string.format(condition=condition)

    def select(self, selection_string):
        """Execute a selection against the topology

        Parameters
        ----------
        selection_string : str
            An expression in the MDTraj atom selection DSL

        Examples
        --------
        >>> topology.select('name O and water')
        array([1, 3, 5, 10, ...])

        Returns
        -------
        indices : np.ndarray, dtype=int, ndim=1
            Array of the indices of the atoms matching the selection expression.

        See Also
        --------
        select_expression, mdtraj.core.selection.parse_selection
        """

        filter_func = parse_selection(selection_string).expr
        indices = np.array([a.index for a in self.atoms if filter_func(a)])
        return indices

    def select_atom_indices(self, selection='minimal'):
        """Get the indices of biologically-relevant groups by name.

        Parameters
        ----------
        selection : {'all', 'alpha', 'minimal', 'heavy', 'water'}
            What types of atoms to select.

            ``all``
                All atoms
            ``alpha``
                Protein residue alpha carbons
            ``minimal``
                Keep the atoms in protein residues with names in {CA, CB, C, N, O}
            ``heavy``
                All non-hydrogen protein atoms.
            ``water``
                Water oxygen atoms

        Returns
        ----------
        indices : np.ndarray (N,)
            An array of the indices of the selected atoms.
        """
        selection = selection.lower()
        options = ['all', 'alpha', 'minimal', 'heavy', 'water']
        if selection == 'all':
            atom_indices = np.arange(self.n_atoms)
        elif selection == 'alpha':
            atom_indices = [a.index for a in self.atoms if
                            a.name == 'CA'
                            and a.residue.is_protein]
        elif selection == 'minimal':
            atom_indices = [a.index for a in self.atoms if
                            a.name in ['CA', 'CB', 'C', 'N', 'O']
                            and a.residue.is_protein]
        elif selection == 'heavy':
            atom_indices = [a.index for a in self.atoms if
                            a.element != elem.hydrogen
                            and a.residue.is_protein]
        elif selection == 'water':
            atom_indices = [a.index for a in self.atoms if
                            a.name in ['O', 'OW']
                            and a.residue.is_water]
        else:
            raise ValueError(
                '%s is not a valid option. Selection must be one of %s' % (
                    selection, ', '.join(options)))

        indices = np.array(atom_indices)
        return indices

    def select_pairs(self, selection1=None, selection2=None):
        """Generate unique pairs of atom indices.

        If a selecton is a string, it will be resolved using the atom selection
        DSL, otherwise it is expected to be an array of atom indices.

        Parameters
        ----------
        selection1 : str or array-like, shape=(n_indices, ), dtype=int
            A selection for `select()` or an array of atom indices.
        selection2 : str or array-like, shape=(n_indices, ), dtype=int
            A selection for `select()` or an array of atom indices.

        Returns
        -------
        pairs : array-like, shape=(n_pairs, 2), dtype=int
            Each row gives the indices of two atoms.

        """
        # Resolve selections using the atom selection DSL...
        if isinstance(selection1, string_types):
            a_indices = self.select(selection1)
        else:  # ...or use a provided array of indices.
            a_indices = ensure_type(selection1, dtype=np.int32, ndim=1,
                                    name='a_indices', warn_on_cast=False)
        if isinstance(selection2, string_types):
            b_indices = self.select(selection2)
        else:
            b_indices = ensure_type(selection2, dtype=np.int32, ndim=1,
                                    name='b_indices', warn_on_cast=False)
        a_indices.sort()
        b_indices.sort()

        # Create unique pairs from the indices.
        # In the cases where a_indices and b_indices are identical or mutually
        # exclusive, we can utilize a more efficient and memory friendly
        # approach by removing the intermediate set creation required in
        # the general case.
        if np.array_equal(a_indices, b_indices):
            pairs = self._unique_pairs_equal(a_indices)
        elif len(np.intersect1d(a_indices, b_indices)) == 0:
            pairs = self._unique_pairs_mutually_exclusive(a_indices, b_indices)
        else:
            pairs = self._unique_pairs(a_indices, b_indices)
        return pairs

    @classmethod
    def _unique_pairs(cls, a_indices, b_indices):
        return np.array(list(set(
            (a, b) if a > b else (b, a)
            for a, b in itertools.product(a_indices, b_indices)
            if a != b)), dtype=np.int32)

    @classmethod
    def _unique_pairs_mutually_exclusive(cls, a_indices, b_indices):
        pairs = np.fromiter(itertools.chain.from_iterable(
            itertools.product(a_indices, b_indices)),
            dtype=np.int32, count=len(a_indices) * len(b_indices) * 2)
        return np.vstack((pairs[::2], pairs[1::2])).T

    @classmethod
    def _unique_pairs_equal(cls, a_indices):
        pairs = np.fromiter(itertools.chain.from_iterable(
            itertools.combinations(a_indices, 2)),
            dtype=np.int32, count=len(a_indices) * (len(a_indices) - 1))
        return np.vstack((pairs[::2], pairs[1::2])).T

    def find_molecules(self):
        """Identify molecules based on bonds.

        A molecule is defined as a set of atoms that are connected to each other by bonds.
        This method uses the list of bonds to divide up the Topology's atoms into molecules.

        Returns
        -------
        molecules : list of sets
            Each entry represents one molecule, and is the set of all Atoms in that molecule
        """
        if len(self._bonds) == 0 and any(res.n_atoms > 1 for res in self._residues):
            raise ValueError('Cannot identify molecules because this Topology does not include bonds')

        # Make a list of every other atom to which each atom is connected.

        num_atoms = self.n_atoms
        atom_bonds = [[] for i in range(num_atoms)]
        for atom1, atom2 in self.bonds:
            atom_bonds[atom1.index].append(atom2.index)
            atom_bonds[atom2.index].append(atom1.index)

        # This is essentially a recursive algorithm, but it is reformulated as a loop to avoid
        # stack overflows.  It selects an atom, marks it as a new molecule, then recursively
        # marks every atom bonded to it as also being in that molecule.

        atom_molecule = [-1]*num_atoms
        num_molecules = 0
        for i in range(num_atoms):
            if atom_molecule[i] == -1:
                # Start a new molecule.

                atom_stack = [i]
                neighbor_stack = [0]
                molecule = num_molecules
                num_molecules += 1

                # Recursively tag all the bonded atoms.

                while len(atom_stack) > 0:
                    atom = atom_stack[-1]
                    atom_molecule[atom] = molecule
                    while neighbor_stack[-1] < len(atom_bonds[atom]) and atom_molecule[atom_bonds[atom][neighbor_stack[-1]]] != -1:
                        neighbor_stack[-1] += 1
                    if neighbor_stack[-1] < len(atom_bonds[atom]):
                        atom_stack.append(atom_bonds[atom][neighbor_stack[-1]])
                        neighbor_stack.append(0)
                    else:
                        del atom_stack[-1]
                        del neighbor_stack[-1]

        # Build the final output.

        molecules = [set() for i in range(num_molecules)]
        for atom in self.atoms:
            molecules[atom_molecule[atom.index]].add(atom)
        return molecules

    def guess_anchor_molecules(self):
        """Guess anchor molecules for imaging

        Returns
        -------
        anchor_molecules : list of atom sets
            List of anchor molecules

        See Also
        --------
        Trajectory.image_molecules
        """
        molecules = self.find_molecules()

        # Select the anchor molecules.
        molecules.sort(key=lambda x: -len(x))
        atoms_cutoff = max(len(molecules[int(0.1*len(molecules))]),
                           int(0.1*len(molecules[0])))
        anchor_molecules = [mol for mol in molecules if len(mol) > atoms_cutoff]
        num_anchors = len(anchor_molecules)
        if num_anchors == 0:
            raise ValueError("Could not find any anchor molecules. Based on "
                             "our heuristic, those should be molecules with "
                             "more than {} atoms. Perhaps your topology "
                             "doesn't give an accurate bond graph?"
                             .format(atoms_cutoff))
        return anchor_molecules


class Chain(object):
    """A Chain object represents a chain within a Topology.

    Attributes
    ----------
    index : int
        The index of the Chain within its Topology
    topology : mdtraj.Topology
        The Topology this Chain belongs to
    chain_id : str
        The PDB chainID of the Chain.
    residues : generator
        Iterator over all Residues in the Chain.
    atoms : generator
        Iterator over all Atoms in the Chain.
    """

    def __init__(self, index, topology, chain_id=None):
        """Construct a new Chain.  You should call add_chain() on the Topology instead of calling this directly."""
        # The index of the Chain within its Topology
        self.index = index
        # The Topology this Chain belongs to
        self.topology = topology
        self._residues = []
        # PDB format chainID
        self.chain_id = chain_id

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
        """Get a specific residue in this Chain.

        Parameters
        ----------
        index : int
            The index of the residue to select.

        Returns
        -------
        residue : Residue
        """
        return self._residues[index]

    @property
    def n_residues(self):
        """Get the number of residues in this Chain. """
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

    def atoms_by_name(self, name):
        """Iterator over all Atoms in the Chain with a specified name.

        Parameters
        ----------
        name : str
            The particular atom name of interest.

        Examples
        --------
        >>> for atom in chain.atoms_by_name('CA'):
        ...     print(atom)

        Returns
        -------
        atomiter : generator
        """
        for atom in self.atoms:
            if atom.name == name:
                yield atom

    def atom(self, index):
        """Get a specific atom in this Chain.

        Parameters
        ----------
        index : int
            The index of the atom to select.

        Returns
        -------
        atom : Atom
        """
        # this could be made faster by caching the list
        # of atoms internally if necessary
        return next(itertools.islice(self.atoms, index, index + 1))

    @property
    def n_atoms(self):
        """Get the number of atoms in this Chain"""
        return sum(r.n_atoms for r in self._residues)

    def __hash__(self):
        return self.index


class Residue(object):
    """A Residue object represents a residue within a Topology.

    Attributes
    ----------
    name : str
        The name of the Residue
    index : int
        The index of the Residue within its Topology
    chain : mdtraj.topology.Chain
        The chain within which this residue belongs
    resSeq : int
        The residue sequence number
    segment_id : str, optional
        A label for the segment to which this residue belongs
    """

    def __init__(self, name, index, chain, resSeq, segment_id=''):
        """Construct a new Residue.  You should call add_residue()
        on the Topology instead of calling this directly."""
        self.name = name
        self.index = index
        self.chain = chain
        self.resSeq = resSeq
        self.segment_id = segment_id
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

    def atoms_by_name(self, name):
        """Iterator over all Atoms in the Residue with a specified name

        Parameters
        ----------
        name : str
            The particular atom name of interest.

        Examples
        --------
        >>> for atom in residue.atoms_by_name('CA'):
        ...     print(atom)

        Returns
        -------
        atomiter : generator
        """
        for atom in self.atoms:
            if atom.name == name:
                yield atom

    def atom(self, index_or_name):
        """Get a specific atom in this Residue.

        Parameters
        ----------
        index_or_name : {int, str}
            Either a (zero-based) index, or the name of the atom. If a string
            is passed in, the first atom -- in index order -- with a matching
            name wil be returned.

        Returns
        -------
        atom : Atom
        """
        try:
            return self._atoms[index_or_name]
        except TypeError:
            try:
                return next(self.atoms_by_name(index_or_name))
            except StopIteration:
                raise KeyError('no matching atom found')

    @property
    def n_atoms(self):
        """Get the number of atoms in this Residue"""
        return len(self._atoms)

    @property
    def is_protein(self):
        """Whether the residue is one found in proteins."""
        return self.name in _PROTEIN_RESIDUES

    @property
    def code(self):
        """Get the one letter code for this Residue"""
        if self.is_protein:
            return _AMINO_ACID_CODES[self.name]
        else:
            return None

    @property
    def is_water(self):
        """Whether the residue is water.

        Residue names according to VMD

        References
        ----------
        http://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node133.html
        """
        return self.name in _WATER_RESIDUES

    @property
    def is_nucleic(self):
        """Whether the residue is one found in nucleic acids."""
        raise NotImplementedError

    def __str__(self):
        return '%s%s' % (self.name, self.resSeq)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash((self.name, self.index, self.resSeq, self.segment_id))


class Atom(object):
    """An Atom object represents a residue within a Topology.

    Attributes
    ----------
    name : str
        The name of the Atom
    element : mdtraj.element.Element
        The element of the Atoms
    index : int
        The index of the Atom within its Topology
    residue : mdtraj.topology.Residue
        The Residue this Atom belongs to
    serial : int
        The serial number from the PDB specification. Unlike index,
        this may not be contiguous or 0-indexed.
    """

    def __init__(self, name, element, index, residue, serial=None):
        """Construct a new Atom.  You should call add_atom() on the Topology instead of calling this directly."""
        # The name of the Atom
        self.name = name
        # That Atom's element
        self.element = element
        # The index of the Atom within its Topology
        self.index = index
        # The Residue this Atom belongs to
        self.residue = residue
        # The not-necessarily-contiguous "serial" number from the PDB spec
        self.serial = serial

    @property
    def n_bonds(self):
        """Number of bonds in which the atom participates."""
        # TODO: this info could be cached.
        return ilen(bond for bond in self.residue.chain.topology.bonds
                    if self in bond)

    @property
    def is_backbone(self):
        """Whether the atom is in the backbone of a protein residue"""
        return (self.name in set(['C', 'CA', 'N', 'O'])
                and self.residue.is_protein)

    @property
    def is_sidechain(self):
        """Whether the atom is in the sidechain of a protein residue"""
        return (self.name not in set(['C', 'CA', 'N', 'O', 'HA', 'H'])
                and self.residue.is_protein)

    @property
    def segment_id(self):
        """User specified segment_id of the residue to which this atom belongs"""
        return self.residue.segment_id

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
        """A quick comparison. """
        return self.index

    def __str__(self):
        return '%s-%s' % (self.residue, self.name)

    def __repr__(self):
        return str(self)


# Enumerated values for bond type

class Single(Singleton):
    def __repr__(self):
        return 'Single'

    def __float__(self):
        return 1.0
Single = Single()


class Double(Singleton):
    def __repr__(self):
        return 'Double'

    def __float__(self):
        return 2.0
Double = Double()


class Triple(Singleton):
    def __repr__(self):
        return 'Triple'

    def __float__(self):
        return 3.0
Triple = Triple()


class Aromatic(Singleton):
    def __repr__(self):
        return 'Aromatic'

    def __float__(self):
        return 1.5
Aromatic = Aromatic()


class Amide(Singleton):
    def __repr__(self):
        return 'Amide'

    def __float__(self):
        return 1.25
Amide = Amide()


def float_to_bond_type(bond_float):
    """
    Convert a float to known bond type class, or None if no matched class if found

    Parameters
    ----------
    bond_float : float
        Representation of built in bond types as a float,
        Maps this float to specific bond class, if the float has no map, None is returned instead

    Returns
    -------
    bond_type : mdtraj.topology.Singleton subclass or None
        Bond type matched to the float if known
        If no match is found, returns None (which is also a valid type for the Bond class)
    """
    all_bond_types = [Single, Double, Triple, Aromatic, Amide]
    for bond_type in all_bond_types:
        if float(bond_type) == float(bond_float):
            return bond_type
    return None


class Bond(namedtuple('Bond', ['atom1', 'atom2'])):
    """A Bond representation of a bond between two Atoms within a Topology

    Attributes
    ----------
    atom1 : mdtraj.topology.Atom
        The first atom in the bond
    atom2 : mdtraj.topology.Atom
        The second atom in the bond
    order : instance of mdtraj.topology.Singleton or None
    type : int on [1,3] domain or None
    """

    def __new__(cls, atom1, atom2, type=None, order=None):
        """Construct a new Bond.  You should call add_bond()
        on the Topology instead of calling this directly.

        Must use __new__ constructor since this is an immutable class
        """
        bond = super(Bond, cls).__new__(cls, atom1, atom2)
        assert isinstance(type, Singleton) or type is None, "Type must be None or a Singleton"
        assert order is None or 1 <= order <= 3, "Order must be int between 1 to 3 or None"
        bond.type = type
        bond.order = order
        return bond

    def __getnewargs__(self):
        """
        Support for pickle protocol 2:
        http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances
        """
        return self[0], self[1], self.type, self.order

    @property
    def _equality_tuple(self):
        # Hierarchy of parameters: Atom1 index -> Atom2 index -> type -> order
        return (self[0].index, self[1].index,
                float(self.type) if self.type is not None else 0.0,
                self.order if self.order is not None else 0)

    def __deepcopy__(self, memo):
        return Bond(self[0], self[1], self.type, self.order)

    def __eq__(self, other):
        if not isinstance(other, Bond):
            return False
        return self._equality_tuple == other._equality_tuple

    def __repr__(self):
        s = "Bond(%s, %s" % (self[0], self[1])
        if self.type is not None:
            s = "%s, type=%s" % (s, self.type)
        if self.order is not None:
            s = "%s, order=%d" % (s, self.order)
        s += ")"
        return s

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        # Set of atoms making up bonds, the type, and the order
        return hash((self[0], self[1], self.type, self.order))

    @staticmethod
    def _other_is_bond(other):
        # Ensure other type for inequalities is a bond
        if not isinstance(other, Bond):
            raise TypeError("Bond inequalities can only be compared with other bonds")

    def __gt__(self, other):
        # Cannot use total_ordering because namedtuple
        # has its own __gt__, __lt__, etc. methods, which
        # supersede total_ordering
        self._other_is_bond(other)
        return self._equality_tuple > other._equality_tuple

    def __ge__(self, other):
        self._other_is_bond(other)
        return self._equality_tuple >= other._equality_tuple

    def __lt__(self, other):
        self._other_is_bond(other)
        return self._equality_tuple < other._equality_tuple

    def __le__(self, other):
        self._other_is_bond(other)
        return self._equality_tuple <= other._equality_tuple

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getstate__(self):
        # This is required for pickle because the parent class
        # does not properly return a state
        return self.__dict__
