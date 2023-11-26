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
from xml.etree import cElementTree

import numpy as np

from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils import ilen, import_, ensure_type
from mdtraj.core.element import virtual_site


__all__ = ['load_hoomdxml']


@FormatRegistry.register_loader('.hoomdxml')
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
    filename : path-like
        The path on disk to the XML file
    top : None
        This argumet is ignored

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object, with corresponding 
        Topology.

    Notes
    -----
    This function requires the NetworkX python package.
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
    box.attrib = dict((key.lower(), val) for key, val in box.attrib.items())
    # be generous for case of box attributes
    lx = float(box.attrib['lx'])
    ly = float(box.attrib['ly'])
    lz = float(box.attrib['lz'])
    try:
        xy = float(box.attrib['xy'])
        xz = float(box.attrib['xz'])
        yz = float(box.attrib['yz'])
    except (ValueError, KeyError):
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
    if hasattr(bond, 'text'):
        bonds = [(int(b.split()[1]), int(b.split()[2])) for b in bond.text.splitlines()[1:]]
        chains = _find_chains(bonds)
    else:
        chains = []
        bonds = []

    # Relate the first index in the bonded-group to mdtraj.Residue
    bonded_to_residue = {}
    for i, _ in enumerate(types):
        bonded_group = _in_chain(chains, i)
        if bonded_group is not None:
            if bonded_group[0] not in bonded_to_residue:
                t_chain = topology.add_chain()
                t_residue = topology.add_residue('A', t_chain)
                bonded_to_residue[bonded_group[0]] = t_residue
            topology.add_atom(types[i], virtual_site, 
                    bonded_to_residue[bonded_group[0]])
        if bonded_group is None:
            t_chain = topology.add_chain()
            t_residue = topology.add_residue('A', t_chain)
            topology.add_atom(types[i], virtual_site, t_residue)

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

    Notes
    -----
    This function requires the NetworkX python package.
    """
    nx = import_('networkx')
    chains = []
    bond_list = np.asarray(bond_list)
    molecules = nx.Graph()
    molecules.add_nodes_from(set(bond_list.flatten()))
    molecules.add_edges_from(bond_list)
    return [sorted(x) for x in list(nx.connected_components(molecules))]

def _in_chain(chains, atom_index):
    """Check if an item is in a list of lists"""
    for chain in chains:
        if atom_index in chain:
            return chain
    return None
