##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: TJ Lane
# Contributors: Robert McGibbon, Jason Swails
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
"""Load an md.Topology from AMBER PRMTOP files
"""

# Written by: TJ Lane <tjlane@stanford.edu> 2/25/14
# This code was mostly stolen/stripped down from OpenMM code, specifically
# the files amber_file_parser.py and amberprmtopfile.py


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import re

from mdtraj.core import topology
from mdtraj.formats import pdb
from mdtraj.core import element as elem

FORMAT_RE_PATTERN = re.compile("([0-9]+)([a-zA-Z]+)([0-9]+)\.?([0-9]*)")

__all__ = ['load_prmtop']

##############################################################################
# Functions
##############################################################################


def _get_pointer_value(pointer_label, raw_data):
   
    POINTER_LABELS = """
    NATOM, NTYPES, NBONH, MBONA, NTHETH, MTHETA,
    NPHIH, MPHIA, NHPARM, NPARM, NEXT, NRES,
    NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA,
    NATYP, NPHB, IFPERT, NBPER, NGPER, NDPER,
    MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP
    """
    
    POINTER_LABEL_LIST = POINTER_LABELS.replace(',', '').split()
    
    index = POINTER_LABEL_LIST.index(pointer_label)
    return float(raw_data['POINTERS'][index])


def load_prmtop(filename, **kwargs):
    """Load an AMBER prmtop topology file from disk.

    Parameters
    ----------
    filename : path-like
        Path to the prmtop file on disk.

    Returns
    -------
    top : md.Topology
        The resulting topology, as an md.Topology object.

    Notes
    -----
    Deprecated fields in the prmtop file are not loaded. This includes the
    BOX dimensions, which should be stored in trajectory files instead of the
    prmtop for systems with periodic boundary conditions. Because '.binpos'
    files do not store box dimensions, this means that unitcell information
    will be lost if you use .binpos + .prmtop files with MDTraj.
    
    Examples
    --------
    >>> topology = md.load_prmtop('mysystem.prmtop')
    >>> # or
    >>> trajectory = md.load('trajectory.mdcrd', top='system.prmtop')
    """
    top = topology.Topology()

    prmtop_version = None
    flags      = []
    raw_format = {}
    raw_data   = {}
    ignoring = False

    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '%':
                if line.startswith('%VERSION'):
                    tag, prmtop_version = line.rstrip().split(None, 1)

                elif line.startswith('%FLAG'):
                    tag, flag = line.rstrip().split(None, 1)
                    flags.append(flag)
                    raw_data[flag] = []
                    ignoring = flag in ('TITLE', 'CTITLE')

                elif line.startswith('%FORMAT'):
                    format = line.rstrip()
                    index0=format.index('(')
                    index1=format.index(')')
                    format = format[index0+1:index1]
                    m = FORMAT_RE_PATTERN.search(format)
                    if m is None:
                        ignoring = True
                    else:
                        raw_format[flags[-1]] = (format, m.group(1), m.group(2), m.group(3), m.group(4))

                elif line.startswith('%COMMENT'):
                    continue

            elif not ignoring:
                flag=flags[-1]
                format, numItems, itemType, itemLength, itemPrecision = raw_format[flag]
                iLength=int(itemLength)
                line = line.rstrip()
                for index in range(0, len(line), iLength):
                    item = line[index:index+iLength]
                    if item:
                        raw_data[flag].append(item.strip())

    # Add atoms to the topology

    pdb.PDBTrajectoryFile._loadNameReplacementTables()
    previous_residue = None
    c = top.add_chain()

    n_atoms = int(_get_pointer_value('NATOM', raw_data))

    # built a dictionary telling us which atom belongs to which residue
    residue_pointer_dict = {}
    res_pointers = raw_data['RESIDUE_POINTER']        
    first_atom = [int(p)-1 for p in res_pointers] # minus 1 necessary
    first_atom.append(n_atoms)
    res = 0
    for i in range(n_atoms):
        while first_atom[res+1] <= i:
            res += 1
        residue_pointer_dict[i] = res

    # add each residue/atom to the topology object
    for index in range(n_atoms):

        res_number = residue_pointer_dict[index]
        if res_number != previous_residue:

            previous_residue = res_number

            # check
            res_name = raw_data['RESIDUE_LABEL'][residue_pointer_dict[index]].strip()
            if res_name in pdb.PDBTrajectoryFile._residueNameReplacements:
                res_name = pdb.PDBTrajectoryFile._residueNameReplacements[res_name]
            r = top.add_residue(res_name, c)

            if res_name in pdb.PDBTrajectoryFile._atomNameReplacements:
                atom_replacements = pdb.PDBTrajectoryFile._atomNameReplacements[res_name]
            else:
                atom_replacements = {}

        atom_name = raw_data['ATOM_NAME'][index].strip()
        if atom_name in atom_replacements:
            atom_name = atom_replacements[atom_name]

        # Get the element from the prmtop file if available
        if 'ATOMIC_NUMBER' in raw_data:
            try:
                element = elem.Element.getByAtomicNumber(int(raw_data['ATOMIC_NUMBER'][index]))
            except KeyError:
                element = elem.virtual
        else:
            # Try to guess the element from the atom name.

            upper = atom_name.upper()
            if upper.startswith('CL'):
                element = elem.chlorine
            elif upper.startswith('NA'):
                element = elem.sodium
            elif upper.startswith('MG'):
                element = elem.magnesium
            elif upper.startswith('ZN'):
                element = elem.zinc
            else:
                try:
                    element = elem.get_by_symbol(atom_name[0])
                except KeyError:
                    element = elem.virtual

        top.add_atom(atom_name, element, r)

    # Add bonds to the topology
    bond_pointers = raw_data["BONDS_INC_HYDROGEN"] + raw_data["BONDS_WITHOUT_HYDROGEN"]
    atoms = list(top.atoms)

    bond_list = []
    for ii in range(0,len(bond_pointers),3):
        if int(bond_pointers[ii])<0 or int(bond_pointers[ii+1])<0:
            raise Exception("Found negative bonded atom pointers %s"
                             % ((bond_pointers[ii],
                                 bond_pointers[ii+1]),))

        else:
            bond_list.append((int(bond_pointers[ii])//3, int(bond_pointers[ii+1])//3))

    for bond in bond_list:
        top.add_bond(atoms[bond[0]], atoms[bond[1]])

    return top
