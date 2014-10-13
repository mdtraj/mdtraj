##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Teng Lin
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

#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------
#  Open-Source PyMOL Copyright Notice
#  ==================================
# 
#  The Open-Source PyMOL source code is copyrighted, but you can freely
#  use and copy it as long as you don't change or remove any of the
#  Copyright notices. The Open-Source PyMOL product is made available
#  under the following open-source license terms:
# 
#  ----------------------------------------------------------------------
#  Open-Source PyMOL is Copyright (C) Schrodinger, LLC.
# 
#  All Rights Reserved
# 
#  Permission to use, copy, modify, distribute, and distribute modified
#  versions of this software and its built-in documentation for any
#  purpose and without fee is hereby granted, provided that the above
#  copyright notice appears in all copies and that both the copyright
#  notice and this permission notice appear in supporting documentation,
#  and that the name of Schrodinger, LLC not be used in advertising or
#  publicity pertaining to distribution of the software without specific,
#  written prior permission.
# 
#  SCHRODINGER, LLC DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN
#  NO EVENT SHALL SCHRODINGER, LLC BE LIABLE FOR ANY SPECIAL, INDIRECT OR
#  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
#  OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
#  OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE
#  USE OR PERFORMANCE OF THIS SOFTWARE.
#  ----------------------------------------------------------------------
# 
#  PyMOL Trademark Notice
#  ======================
# 
#  PyMOL(TM) is a trademark of Schrodinger, LLC. Derivative
#  software which contains PyMOL source code must be plainly
#  distinguished from any and all PyMOL products distributed by Schrodinger,
#  LLC in all publicity, advertising, and documentation.
# 
#  The slogans, "Includes PyMOL(TM).", "Based on PyMOL(TM) technology.",
#  "Contains PyMOL(TM) source code.", and "Built using PyMOL(TM).", may
#  be used in advertising, publicity, and documentation of derivative
#  software provided that the notice, "PyMOL is a trademark of Schrodinger,
#  LLC.", is included in a footnote or at the end of the
#  document.
# 
#  All other endorsements employing the PyMOL trademark require specific,
#  written prior permission.
# 

"""
Load an md.Topology from Schrodinger mae files.
"""

import re
import string
import copy
import numpy as np
from pandas import DataFrame
from mdtraj.utils import import_
from mdtraj.formats.registry import _FormatRegistry
from mdtraj.core.trajectory import Trajectory
from mdtraj.core.topology import Topology

__all__ = ['load_mae', "mae_to_dataframes"]


strip_re = re.compile(r'\#.*\#')
token_re = re.compile(r'"([^"]*)"|([^ ]+)')
array_re = re.compile(r'(.*)\[([0-9]+)\]')

coerce = {
    's' : str,
    'i' : int,
    'r' : float,
    'b': bool}

def _convert(func, tok):
    try:
        return apply(func,(tok,))
    except:
        return None

_atomnum2element = {
    1:"H",
    2:"He",
    3:"Li",
    4:"Be",
    5:"B",
    6:"C",
    7:"N",
    8:"O",
    9:"F",
    10:"Ne",
    11:"Na",
    12:"Mg",
    13:"Al",
    14:"Si",
    15:"P",
    16:"S",
    17:"Cl",
    18:"Ar",
    19:"K",
    20:"Ca",
    21:"Sc",
    22:"Ti",
    23:"V",
    24:"Cr",
    25:"Mn",
    26:"Fe",
    27:"Co",
    28:"Ni",
    29:"Cu",
    30:"Zn",
    31:"Ga",
    32:"Ge",
    33:"As",
    34:"Se",
    35:"Br",
    36:"Kr",
    37:"Rb",
    38:"Sr",
    39:"Y",
    40:"Zr",
    41:"Nb",
    42:"Mo",
    43:"Tc",
    44:"Ru",
    45:"Rh",
    46:"Pd",
    47:"Ag",
    48:"Cd",
    49:"In",
    50:"Sn",
    51:"Sb",
    52:"Te",
    53:"I",
    54:"Xe",
    55:"Cs",
    56:"Ba",
    57:"La",
    58:"Ce",
    59:"Pr",
    60:"Nd",
    61:"Pm",
    62:"Sm",
    63:"Eu",
    64:"Gd",
    65:"Tb",
    66:"Dy",
    67:"Ho",
    68:"Er",
    69:"Tm",
    70:"Yb",
    71:"Lu",
    72:"Hf",
    73:"Ta",
    74:"W",
    75:"Re",
    76:"Os",
    77:"Ir",
    78:"Pt",
    79:"Au",
    80:"Hg",
    81:"Tl",
    82:"Pb",
    83:"Bi",
    84:"Po",
    85:"At",
    86:"Rn",
    87:"Fr",
    88:"Ra",
    89:"Ac",
    90:"Th",
    91:"Pa",
    92:"U",
    93:"Np",
    94:"Pu",
    95:"Am",
    96:"Cm",
    97:"Bk",
    98:"Cf",
    99:"Es",
    100:"Fm",
    101:"Md",
    102:"No",
    103:"Lr",
    104:"Rf",
    105:"Db",
    106:"Sg",
    107:"Bh",
    108:"Hs",
    109:"Mt",
    110:"Ds",
    111:"Rg",
    112:"Uub",
    113:"Uut",
    114:"Uuq",
    115:"Uup",
    116:"Uuh",
}

class _MAEParser:

    def __init__(self,lst=None):
        self.i = 0 # index in list
        self.t = [] # token list
        self.d = [] # hiearchy of data read
        self.lst = lst
        self.lst_len = len(lst)

    def nxt_lin(self):
        if self.lst:
            if self.i<self.lst_len:
                self.i = self.i + 1
                return self.lst[self.i-1]
        return None

    def nxt_tok(self):
        while 1:
            if len(self.t):
                return self.t.pop(0)
            else:
                l = self.nxt_lin()
                if not l:
                    return None
                l = string.strip(strip_re.sub('',l))
                self.t = token_re.findall(l)
                self.t = map(lambda x:string.join(x,''),self.t)
        return None

    def push_tok(self,tok):
        self.t.insert(0,tok)

    def parse_top(self):
        dct = {}
        stk = [] # keyword stack
        mode = 0 # 0 = definition, 1 = data
        while 1:
            tok = self.nxt_tok()
            if tok==None:
                break
            if tok==':::':
                mode = 1
            if not len(stk):
                mode = 0
            if mode:
                lab = stk.pop(0)
                #dct[lab] = apply(coerce[lab[0]],(tok,))
                if tok==':::':
                    mode = 0
                else:
                    dct[lab] = _convert(coerce[lab[0]],tok)
            else:
                stk.append(tok)
            if tok=='}':
                break
        return dct

    def parse_array(self,n_rec): # creates a list of homogenous lists
                                          # each containing data for one field
        dct = {}
        data = [] # actual array data
        stk = [] # keyword stack
        coer = [] # coersion functions
        n_fld = 0
        mode = 0 # 0 = definition, 1 = data
        cc = 0
        while 1:
            tok = self.nxt_tok()
            if tok==None:
                break
            if tok=='}':
                break
            if tok==':::':
                if not mode:
                    mode = 1
                    n_fld = len(stk)
                    c = 0
                    for a in stk: # create row index for the array
                        dct[a] = c
                        data.append([]) # add row for each field
                        coer.append(coerce[a[0]])
                        c = c + 1
            elif not mode:
                stk.append(tok)
            else: # here we actually read the array
                self.push_tok(tok)
                for cc in range(n_rec):
                    tok = self.nxt_tok() # chuck index
                    if tok==':::': # truncated/incomplete array
                        break
                    c = 0
                    for c in range(n_fld):
                        tok = self.nxt_tok()
                        if tok==None:
                            break
                        #data[c].append(apply(coer[c],(tok,)))
                        data[c].append(_convert(coer[c],tok))
                        if tok=='}':
                            break
        return (n_rec,dct,data) # return a tuple

    def parse_m_ct(self):
        dct = {}
        stk = [] # keyword stack
        mode = 0 # 0 = definition, 1 = data
        while 1:
            tok = self.nxt_tok()
            if tok==None:
                break
            if tok==':::':
                mode = 1
            elif mode:
                if not len(stk):
                    mode = 0
                    self.push_tok(tok) # go around
                else:
                    dct[stk.pop(0)] = tok
            else:
                arm = array_re.findall(tok)
                if len(arm):
                    arm = arm[0]
                    n_rec=int(arm[1])
                    if arm[0] in ['m_atom','m_bond']:
                        self.nxt_tok() # skip '{'
                        dct[arm[0]]=self.parse_array(n_rec)
                else:
                    stk.append(tok)
            if tok=='}':
                break
        return dct

    def parse(self):
        while 1:
            tok = self.nxt_tok()
            if tok==None:
                break
            if tok=='{':
                self.d.append(['top',self.parse_top()])
            elif tok in ['f_m_ct','p_m_ct']:
                self.nxt_tok() # skip '{'
                self.d.append([tok,self.parse_m_ct()])
        return self.d

def _handle_atom_number(col):
    return [_atomnum2element[n] for n in col]



_atom_prop_map = {
    'x': ['r_m_x_coord'],
    'y': ['r_m_y_coord'],
    'z': ['r_m_z_coord'],
    'name' : ['s_m_atom_name', 's_m_pdb_atom_name'],
    'resName': ['s_m_pdb_residue_name'],
    'resSeq':['i_m_residue_number'],
    'chainID': ['s_m_chain_name'],
    'element': ['i_m_atomic_number']
}

_atom_handler = {
    'i_m_atomic_number': _handle_atom_number
}  

def _read_mae_atom( m_atom):
    """
    """
  
    at_ent = m_atom[1]
    at_dat = m_atom[2]
    d = {}
    for prop, mae_prop in _atom_prop_map.items():
        for p in mae_prop:
            if at_ent.has_key(p):
                val = at_dat[at_ent[p]]
                if p in _atom_handler:
                    func = _atom_handler[p]
                    val = func(val)
                d[prop] = val
                break

    data_frame = DataFrame(d)
    data_frame["serial"] = data_frame.index
    return data_frame



def _read_mae_bond(m_bond):
    bd_ent = m_bond[1]
    bd_dat = m_bond[2]
    d = {}
    for p in ['i_m_from', 'i_m_to', 'i_m_order']:
        d[p] = bd_dat[bd_ent[p]]

    bond_frame = DataFrame(d)
    bond_frame = bond_frame.query('i_m_from < i_m_to & i_m_order > 0')
    bonds = bond_frame[["i_m_from", "i_m_to"]].values
    bonds -= 1
    return bonds


def _traj_append(traj, xyz_list):

    xyz = np.array(xyz_list)
    xyz_shape = xyz.shape
    traj_shape = traj.xyz.shape

    if xyz_shape[1] != traj_shape[1] or xyz_shape[2] != traj_shape[2]:
        return

    traj.xyz.reshape((traj_shape[0]+xyz_shape[0], traj_shape[1], traj_shape[2]))
    traj.xyz[-1-xyz_shape[0]:-1]=xyz


def parse_mae(file_object):
    """
    
    """

    mp = _MAEParser(file_object.readlines())
    mp_rec = mp.parse()

    full_traj = None
    result = []
    xyz_list = []
    for mp_ent in mp_rec:

        if mp_ent[0] == 'f_m_ct':
            f_m_ct = mp_ent[1]

            if not f_m_ct.has_key('m_atom') or not f_m_ct.has_key('m_bond'):
                continue

            if xyz_list and full_traj:
                _traj_append(full_traj, xyz_list)
                xyz_list = []

            m_atom = f_m_ct['m_atom']
            atom_frame = _read_mae_atom(m_atom)

            m_bond = f_m_ct['m_bond']
            bonds = _read_mae_bond(m_bond)
            
            top = Topology.from_dataframe(atom_frame, bonds)
        
            xyz = np.array([atom_frame[["x", "y", "z"]].values])
            xyz /= 10.0  # Convert from angstrom to nanometer
        
            traj = Trajectory(xyz, top)            
            full_traj = traj
            result.append(traj)

        elif mp_ent[0]=='p_m_ct' and full_traj != None:

            f_m_ct = mp_ent[1]

            if f_m_ct.has_key('m_atom'):
                m_atom = f_m_ct['m_atom']
                atom_frame = _read_mae_atom(m_atom)
                xyz = np.array([atom_frame[["x", "y", "z"]].values])
                xyz_list.append(xyz)

            else:
                pass

    if xyz_list and full_traj:
        _traj_append(full_traj, xyz_list)

    return result


def parse_cms(self, file_object):
    """
    """
    mp = _MAEParser(file_object.readlines())
    mp_rec = mp.parse()
    
    pass
    

@_FormatRegistry.register_loader('.mae')
def load_mae(filename):
    """Load a Schrodinger maestro file from disk.

    Parameters
    ----------
    filename : str
        Path to the maestro file on disk.

    Returns
    -------
    traj : md.Trajectory
        The resulting topology, as an md.Topology object.

    Notes
    -----
    This function should work on GAFF and sybyl style mae files, but has
    been primarily tested on GAFF mae files.
    This function does NOT accept multi-structure mae files!!!
    The elements are guessed using GAFF atom types or via the atype string.

    Examples
    --------
    >>> traj = md.load_mae('mysystem.mae')
    """

    f = open(filename, 'r')
    traj = parse_mae(f)

    return traj


@_FormatRegistry.register_loader('.cms')
def load_cms(filename):
    """Load a Schrodinger maestro file from disk.

    Parameters
    ----------
    filename : str
        Path to the maestro file on disk.

    Returns
    -------
    traj : md.Trajectory
        The resulting topology, as an md.Topology object.

    Notes
    -----
    This function should work on GAFF and sybyl style mae files, but has
    been primarily tested on GAFF mae files.
    This function does NOT accept multi-structure mae files!!!
    The elements are guessed using GAFF atom types or via the atype string.

    Examples
    --------
    >>> traj = md.load_mae('mysystem.mae')
    """
    pass


if __name__ == '__main__':


    from mdtraj.testing import get_fn
    fn_dtr = get_fn('1h1q_ligand.mae')
    load_mae(fn_dtr)
