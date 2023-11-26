##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors: Robert McGibbon, John D. Chodera
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
"""Load an md.Topology from tripos mol2 files.
"""

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np
import itertools
import re

from mdtraj.utils import import_
from mdtraj.utils.six.moves import cStringIO as StringIO
from mdtraj.formats.registry import FormatRegistry
from mdtraj.core import element as elem

__all__ = ['load_mol2', "mol2_to_dataframes"]

@FormatRegistry.register_loader('.mol2')
def load_mol2(filename):
    """Load a TRIPOS mol2 file from disk.

    Parameters
    ----------
    filename : path-like
        Path to the prmtop file on disk.

    Returns
    -------
    traj : md.Trajectory
        The resulting topology, as an md.Topology object.

    Notes
    -----
    This function should work on GAFF and sybyl style MOL2 files, but has
    been primarily tested on GAFF mol2 files.
    This function does NOT accept multi-structure MOL2 files!!!
    The elements are guessed using GAFF atom types or via the atype string.

    Examples
    --------
    >>> traj = md.load_mol2('mysystem.mol2')
    """
    from mdtraj.core.trajectory import Trajectory
    from mdtraj.core.topology import Topology, Single, Double, Triple, Aromatic, Amide

    atoms, bonds = mol2_to_dataframes(filename)

    atoms_mdtraj = atoms[["name", "resName"]].copy()
    atoms_mdtraj["serial"] = atoms.index

    #Figure out 1 letter element names

    # IF this is a GAFF mol2, this line should work without issues
    atoms_mdtraj["element"] = atoms.atype.map(gaff_elements)
    # If this is a sybyl mol2, there should be NAN (null) values
    if atoms_mdtraj.element.isnull().any():
        # If this is a sybyl mol2, I think this works generally.
        # Argument x is being passed as a list with only one element.
        def to_element(x):
            if isinstance(x, (list, tuple)):
                assert len(x) == 1
                x = x[0]

            if '.' in x:  # orbital-hybridizations in SYBL
                return x.split('.')[0]
            try:
                # check if we can convert the whole str to an Element,
                # if not, we only pass the first letter.
                from mdtraj.core.element import Element
                Element.getBySymbol(x)
            except KeyError:
                return x[0]
            return x
        atoms_mdtraj["element"] = atoms.atype.apply(to_element)

    # Check if elements inferred from atoms.atype are valid
    # If not, try to infer elements from atoms.name
    try:
        atoms_mdtraj['element'].apply(elem.get_by_symbol)
    except KeyError:
        try:
            atoms_mdtraj["element"] = atoms.name.apply(to_element)
            atoms_mdtraj['element'].apply(elem.get_by_symbol)
        except KeyError:
            raise KeyError('Invalid element passed to atoms DataFrame')

    atoms_mdtraj['resSeq'] = atoms['code']
    atoms_mdtraj["chainID"] = np.ones(len(atoms), 'int')

    bond_type_map = {
        '1': Single,
        '2': Double,
        '3': Triple,
        'am': Amide,
        'ar': Aromatic
    }
    if bonds is not None:
        bonds_mdtraj = bonds[["id0", "id1"]].values
        offset = bonds_mdtraj.min()  # Should this just be 1???
        bonds_mdtraj -= offset
        # Create the bond augment information
        n_bonds = bonds_mdtraj.shape[0]
        bond_augment = np.zeros([n_bonds, 2], dtype=float)
        # Add bond type information
        bond_augment[:, 0] = [float(bond_type_map[str(bond_value)]) for bond_value in bonds["bond_type"].values]
        # Add Bond "order" information, this is not known from Mol2 files
        bond_augment[:, 1] = [0.0 for _ in range(n_bonds)]
        # Augment array, dtype is cast to minimal representation of float
        bonds_mdtraj = np.append(bonds_mdtraj, bond_augment, axis=-1)
    else:
        bonds_mdtraj = None

    top = Topology.from_dataframe(atoms_mdtraj, bonds_mdtraj)

    xyzlist = np.array([atoms[["x", "y", "z"]].values])
    xyzlist /= 10.0  # Convert from angstrom to nanometer

    traj = Trajectory(xyzlist, top)

    return traj




def mol2_to_dataframes(filename):
    """Convert a GAFF (or sybyl) mol2 file to a pair of pandas dataframes.

    Parameters
    ----------
    filename : path-like
        Name of mol2 filename

    Returns
    -------
    atoms_frame : pd.DataFrame
        DataFrame containing atom information
    bonds_frame : pd.DataFrame
        DataFrame containing bond information

    Notes
    -----
    These dataframes may contain force field information as well as the
    information necessary for constructing the coordinates and molecular
    topology.  This function has been tested for GAFF and sybyl-style
    mol2 files but has been primarily tested on GAFF mol2 files.
    This function does NOT accept multi-structure MOL2 files!!!

    See Also
    --------
    If you just need the coordinates and bonds, use load_mol2(filename)
    to get a Trajectory object.
    """
    pd = import_('pandas')
    with open(filename) as f:
        data = dict((key, list(grp)) for key, grp in itertools.groupby(f, _parse_mol2_sections))

    # Mol2 can have "status bits" at the end of the bond lines. We don't care
    # about these, but they interfere with using pd_read_table because it looks
    # like one line has too many columns. So we just regex out the offending
    # text.
    status_bit_regex = r"BACKBONE|DICT|INTERRES|\|"
    data["@<TRIPOS>BOND\n"] = [re.sub(status_bit_regex, lambda _: "", s)
                               for s in data["@<TRIPOS>BOND\n"]]

    if len(data["@<TRIPOS>BOND\n"]) > 1:
        csv = StringIO()
        csv.writelines(data["@<TRIPOS>BOND\n"][1:])
        csv.seek(0)
        bonds_frame = pd.read_table(csv, names=["bond_id", "id0", "id1", "bond_type"],
            index_col=0, header=None, sep="\s+", engine='python')
    else:
        bonds_frame = None

    csv = StringIO()
    csv.writelines(data["@<TRIPOS>ATOM\n"][1:])
    csv.seek(0)
    atoms_frame = pd.read_csv(csv, sep="\s+", engine='python',  header=None)
    ncols = atoms_frame.shape[1]
    names=["serial", "name", "x", "y", "z", "atype", "code", "resName", "charge", "status"]
    atoms_frame.columns = names[:ncols]
    
    return atoms_frame, bonds_frame


def _parse_mol2_sections(x):
    """Helper function for parsing a section in a MOL2 file."""
    if x.startswith('@<TRIPOS>'):
        _parse_mol2_sections.key = x
    return _parse_mol2_sections.key




gaff_elements = {
    'br': 'Br',
    'c': 'C',
    'c1': 'C',
    'c2': 'C',
    'c3': 'C',
    'ca': 'C',
    'cc': 'C',
    'cd': 'C',
    'ce': 'C',
    'cf': 'C',
    'cg': 'C',
    'ch': 'C',
    'cl': 'Cl',
    'cp': 'C',
    'cq': 'C',
    'cu': 'C',
    'cv': 'C',
    'cx': 'C',
    'cy': 'C',
    'cz': 'C',
    'f': 'F',
    'h1': 'H',
    'h2': 'H',
    'h3': 'H',
    'h4': 'H',
    'h5': 'H',
    'ha': 'H',
    'hc': 'H',
    'hn': 'H',
    'ho': 'H',
    'hp': 'H',
    'hs': 'H',
    'hw': 'H',
    'hx': 'H',
    'i': 'I',
    'n': 'N',
    'n1': 'N',
    'n2': 'N',
    'n3': 'N',
    'n4': 'N',
    'na': 'N',
    'nb': 'N',
    'nc': 'N',
    'nd': 'N',
    'ne': 'N',
    'nf': 'N',
    'nh': 'N',
    'no': 'N',
    'o': 'O',
    'oh': 'O',
    'os': 'O',
    'ow': 'O',
    'p2': 'P',
    'p3': 'P',
    'p4': 'P',
    'p5': 'P',
    'pb': 'P',
    'px': 'P',
    'py': 'P',
    's': 'S',
    's2': 'S',
    's4': 'S',
    's6': 'S',
    'sh': 'S',
    'ss': 'S',
    'sx': 'S',
    'sy': 'S'}
