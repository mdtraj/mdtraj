##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Tim Moore
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
##############################################################################
import networkx as nx
from xml.etree import cElementTree

import numpy as np
from mdtraj.formats.registry import _FormatRegistry

__all__ = ['load_hoomdxml']


@_FormatRegistry.register_loader('.hoomdxml')
def load_hoomdxml(filename, top=None):
    """Load a single conformation from an HOOMD-Blue XML file.

    For more information on this file format, see:
    http://codeblue.umich.edu/hoomd-blue/doc/page_xml_file_format.html
    Notably, all node names and attributes are in all lower case.
    HOOMD-Blue does not contain residue and chain information explicitly. 
    For this reason, chains will be found by looping over all the bonds and 
    finding what is bonded to what. 
    Each chain consisists of exactly one residue. 

    Parameters
    ----------
    filename : string
        The path on disk to the XML file

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object, with corresponding 
        Topology.
    """
    from mdtraj.core.trajectory import Trajectory
    from mdtraj.core.topology import Topology
    topology = Topology()
    tree = cElementTree.parse(filename)
    config = tree.getroot().find('configuration')
    position = config.find('position')
    bond = config.find('bond')
    atom_type = config.find('type')  # MDTraj calls this "name"

    box = config.find('box')
    # be generous for case of box attributes
    for L in ['lx', 'Lx', 'lX', 'LX']:
        try:
            lx = float(box.attrib[L])
        except KeyError:
            pass
    for L in ['ly', 'Ly', 'lY', 'LY']:
        try:
            ly = float(box.attrib[L])
        except KeyError:
            pass
    for L in ['lz', 'Lz', 'lZ', 'LZ']:
        try:
            lz = float(box.attrib[L])
        except KeyError:
            pass
    try:
        xy = float(box.attrib['xy'])
        xz = float(box.attrib['xz'])
        yz = float(box.attrib['yz'])
    except:
        xy = 0.0
        xz = 0.0
        yz = 0.0
    unitcell_vectors = np.array([[[lx,  xy*ly, xz*lz],
                                  [0.0, ly,    yz*lz],
                                  [0.0, 0.0,   lz   ]]])

    positions, types = [], {}
    for pos in position.text.splitlines()[1:]:
        positions.append((float(pos.split()[0]),
                          float(pos.split()[1]),
                          float(pos.split()[2])))

    for idx, atom_name in enumerate(atom_type.text.splitlines()[1:]):
        types[idx] = str(atom_name.split()[0])
    if len(types) != len(positions):
        raise ValueError('Different number of types and positions in xml file')

    # ignore the bond type
    bonds = []
    for b in bond.text.splitlines()[1:]:
        bonds.append((int(b.split()[1]),  # b.split()[0] == bond type
                      int(b.split()[2])))
    chains = _find_chains(bonds)
    ions = _find_ions(chains, len(types))

    # add chains, bonds and ions (each chain = 1 residue)
    for chain in chains:
        t_chain = topology.add_chain()
        t_residue = topology.add_residue('A', t_chain)
        for atom in chain:
            topology.add_atom(types[atom], 'U', t_residue)
    for ion in ions:
        t_chain = topology.add_chain()
        t_residue = topology.add_residue('A', t_chain)
        topology.add_atom(types[atom], 'U', t_residue)
    for bond in bonds:
        atom1, atom2 = bond[0], bond[1]
        topology.add_bond(topology.atom(atom1), topology.atom(atom2))

    traj = Trajectory(xyz=np.array(positions), topology=topology)
    traj.unitcell_vectors = unitcell_vectors

    return traj

def _find_chains(bond_list):
    """Given a set of bonds, find unique molecules, with the assumption that
    there are no bonds between separate chains (i.e., only INTRAmolecular
    bonds), which also implies that each atom can be in exactly one chain.
    
    Parameters
    ----------
    bond_list : list of (int, int)
        The list of bonds

    Returns
    _______
    chains : list of list of int
        List of atoms in each chain
    """
    chains = []
    bond_list = np.asarray(bond_list)
    molecules = nx.Graph()
    molecules.add_nodes_from(set(bond_list.flatten()))
    molecules.add_edges_from(bond_list)
    return list(nx.connected_components(molecules))

def _find_ions(chains, n):
    """Find atoms that are not included in any chains.

    Parameters
    __________
    chains : list of list of int 
        The list of atoms in chains
    n : int
        Total number of atoms
    """
    ions = []
    for i in range(n):
        if not _in_chain(chains, i):
            ions.append(i)
    return ions

def _in_chain(lists, n):
    """Check if an item is in a list of lists"""
    for l in lists:
        if n in l:
            return True
    return False

if __name__ == '__main__':
    t = load_hoomdxml('start.xml')
