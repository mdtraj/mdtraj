##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Jason Swails
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
#
# Portions of this module originate from the ParmEd program, copyright (c) 2014
# Jason Swails, which is also distributed under the GNU Lesser General Public
# License
#
# Other portions of this code originate from the OpenMM molecular simulation
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
"""Load an md.Topology from CHARMM/XPLOR PSF files"""

# Written by Jason Swails <jason.swails@gmail.com> 9/8/2014
# This code was mostly stolen and stripped down from ParmEd

from mdtraj.core import element as elem
from mdtraj.core import topology
from mdtraj.formats import pdb
from mdtraj.utils.unit import unit_definitions as u

__all__ = ["load_psf"]


class PSFError(Exception):
    """Raised upon problems parsing a PSF file"""

    pass


class _PSFEOF(Exception):
    """Raised when EOF is hit on the parsed PSF"""

    pass


def _convert(string, type, message):
    """Converts a string to the desired data type, making sure to raise PSFError
    with the given message in the case of a failure.

    Parameters
    ----------
    string : str
        String to convert
    type : simple data type
        Either int, float, or str
    message : str
        Message to assign to the PSFError if the conversion fails

    Returns
    -------
    The converted string to the desired datatype
    """
    try:
        return type(string)
    except ValueError:
        raise PSFError(f"Could not convert {message} [{string}]")


def _parse_psf_section(psf):
    """Parses a section of the PSF file, returning the data as a list of integer
    pointers

    Parameters
    ----------
    psf : open file object
        The open PSF file whose pointer is at the beginning of the section to be
        parsed

    Returns
    -------
    (title, pointers, data)

    title : str
        The label of the PSF section we are parsing
    pointers : (int/tuple of ints)
        If one pointer is set, pointers is simply the integer that is the value
        of that pointer. Otherwise it is a tuple with every pointer value
        defined in the first line
    data : list of integers
        A list of all data in the parsed section as integers

    Raises
    ------
    PSFError upon any parsing errors and _PSFEOF on end-of-file
    """
    line = psf.readline()
    while not line.strip():
        if not line:
            raise _PSFEOF("Unexpected EOF in PSF file")
        else:
            line = psf.readline()
    if "!" in line:
        words = line[: line.index("!")].split()
        title = line[line.index("!") + 1 :].strip().upper()
        # Strip out the description
        if ":" in title:
            title = title[: title.index(":")]
    else:
        raise PSFError("Could not determine section title")
    if len(words) == 1:
        pointers = _convert(words[0], int, "pointer")
    else:
        pointers = tuple([_convert(w, int, "pointer") for w in words])
    line = psf.readline().strip()
    if not line and title.startswith("NNB"):
        # This will correctly handle the NNB section (which has a spurious blank
        # line)
        line = psf.readline().strip()
    data = []
    if title == "NATOM" or title == "NTITLE":
        # Store these two sections as strings (ATOM section we will parse
        # later). The rest of the sections are integer pointers
        while line:
            data.append(line)
            line = psf.readline().strip()
    else:
        while line:
            words = line.split()
            data.extend([_convert(w, int, "PSF data") for w in words])
            line = psf.readline().strip()
    return title, pointers, data


def load_psf(fname, **kwargs):
    """Load a CHARMM or XPLOR PSF file from disk

    Parameters
    ----------
    fname : path-like
        Path to the PSF file on disk

    Returns
    -------
    top : md.Topology
        The resulting topology as an md.Topology object

    Notes
    -----
    Only the bond and atom sections are read in, and all atoms are added to the
    same chain in the topology

    Raises
    ------
    PSFError if any parsing errors occur

    Examples
    --------
    >>> topology = md.load_psf('mysystem.psf')
    >>> # or
    >>> trajectory = md.load('trajectory.dcd', top='system.psf')
    """

    top = topology.Topology()

    with open(fname) as f:
        line = f.readline()
        if not line.startswith("PSF"):
            raise PSFError("Unrecognized PSF file.")

        # Store all of the sections and store them in a dict
        f.readline()
        psfsections = dict()
        while True:
            try:
                sec, ptr, data = _parse_psf_section(f)
            except _PSFEOF:
                break
            psfsections[sec] = (ptr, data)
            # We only have to parse up to the NBOND section
            if sec == "NBOND":
                break

    prev_residue = (None, None, None)

    pdb.PDBTrajectoryFile._loadNameReplacementTables()

    natom = _convert(psfsections["NATOM"][0], int, "natom")
    last_chain = None
    for i in range(natom):
        words = psfsections["NATOM"][1][i].split()
        atid = _convert(words[0], int, "atom index")
        if atid != i + 1:
            raise PSFError("Nonsequential atom indices detected!")
        segid = words[1]
        resid = _convert(words[2], int, "residue number")
        rname = words[3]
        name = words[4]
        #       attype = words[5]
        #       charge = _convert(words[6], float, 'partial atomic charge')
        mass = _convert(words[7], float, "atomic mass")
        if last_chain != segid:
            c = top.add_chain()
            last_chain = segid
        curr_residue = (resid, rname, segid)
        if prev_residue != curr_residue:
            prev_residue = curr_residue
            try:
                rname = pdb.PDBTrajectoryFile._residueNameReplacements[rname]
            except KeyError:
                pass
            r = top.add_residue(rname, c, resid, segid)

        try:
            name = pdb.PDBTrajectoryFile._atomNameReplacements[rname][name]
        except KeyError:
            pass

        # Try to guess the element from the atom name for some of the common
        # ions using the names that CHARMM assigns to ions. If it's not one of
        # these 'weird' ion names, look up the element by mass. If the mass is
        # 0, assume a lone pair
        upper = name.upper()
        if upper.startswith("CLA"):
            element = elem.chlorine
        elif upper.startswith("SOD"):
            element = elem.sodium
        elif upper.startswith("POT"):
            element = elem.potassium
        elif upper == "CAL":
            element = elem.calcium
        elif mass == 0:
            element = elem.virtual
        else:
            element = elem.Element.getByMass(mass * u.dalton)
        top.add_atom(name, element, r)

    # Add bonds to the topology
    atoms = list(top.atoms)
    bond_data = psfsections["NBOND"][1]
    nbond = _convert(psfsections["NBOND"][0], int, "number of bonds")
    if len(bond_data) != nbond * 2:
        raise PSFError("Got %d indexes for %d bonds" % (len(bond_data), nbond))

    for i in range(nbond):
        i2 = i * 2
        top.add_bond(atoms[bond_data[i2] - 1], atoms[bond_data[i2 + 1] - 1])

    return top
