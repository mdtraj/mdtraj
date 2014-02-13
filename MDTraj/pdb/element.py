##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christopher M. Bruns
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
##############################################################################


"""
element.py: Used for managing elements.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
Contributors: Robert T. McGibbon

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
from __future__ import print_function, division
import numpy as np

class Element(tuple):
    """An Element represents a chemical element.

    The mdtraj.pdb.element module contains objects for all the standard chemical elements,
    such as element.hydrogen or element.carbon.  You can also call the static method
    Element.getBySymbol() to look up the Element with a particular chemical symbol."""
    __slots__ = []
    _elements_by_symbol = {}
    _elements_by_atomic_number = {}

    def __new__(cls, number, name, symbol, mass):
        """Create a new element

        Parameters
        ----------
        number : int
            The atomic number of the element
        name : str
            The name of the element
        symbol : str
            The chemical symbol of the element
        mass : float
            The atomic mass of the element
        """

        newobj = tuple.__new__(cls, (number, name, symbol, mass))

        # Index this element in a global table
        s = symbol.strip().upper()
        assert s not in Element._elements_by_symbol
        Element._elements_by_symbol[s] = newobj
        if number in Element._elements_by_atomic_number:
            other_element = Element._elements_by_atomic_number[number]
            if mass < other_element.mass:
                # If two "elements" share the same atomic number, they're
                # probably hydrogen and deuterium, and we want to choose
                # the lighter one to put in the table by atomic_number,
                # since it's the "canonical" element.
                Element._elements_by_atomic_number[number] = newobj
        else:
            Element._elements_by_atomic_number[number] = newobj

        return newobj


    @staticmethod
    def getBySymbol(symbol):
        """Get the Element with a particular chemical symbol."""
        s = symbol.strip().upper()
        return Element._elements_by_symbol[s]

    @property
    def number(self):
        return tuple.__getitem__(self, 0)

    @property
    def name(self):
        return tuple.__getitem__(self, 1)

    @property
    def symbol(self):
        return tuple.__getitem__(self, 2)

    @property
    def mass(self):
        return tuple.__getitem__(self, 3)

    def __getitem__(self, item):
        raise TypeError
    
    def __str__(self):
        return self.name

    def atomic_number(self):
        return self.number

# This is for backward compatibility.
def get_by_symbol(symbol):
    s = symbol.strip().upper()
    return Element._elements_by_symbol[s]


hydrogen =       Element(  1, "hydrogen", "H", 1.007947)
deuterium =      Element(  1, "deuterium", "D", 2.01355321270)
helium =         Element(  2, "helium", "He", 4.003)
lithium =        Element(  3, "lithium", "Li", 6.9412)
beryllium =      Element(  4, "beryllium", "Be", 9.0121823)
boron =          Element(  5, "boron", "B", 10.8117)
carbon =         Element(  6, "carbon", "C", 12.01078)
nitrogen =       Element(  7, "nitrogen", "N", 14.00672)
oxygen =         Element(  8, "oxygen", "O", 15.99943)
fluorine =       Element(  9, "fluorine", "F", 18.99840325)
neon =           Element( 10, "neon", "Ne", 20.17976)
sodium =         Element( 11, "sodium", "Na", 22.989769282)
magnesium =      Element( 12, "magnesium", "Mg", 24.30506)
aluminum =       Element( 13, "aluminum", "Al", 26.98153868)
silicon =        Element( 14, "silicon", "Si", 28.08553)
phosphorus =     Element( 15, "phosphorus", "P", 30.9737622)
sulfur =         Element( 16, "sulfur", "S", 32.0655)
chlorine =       Element( 17, "chlorine", "Cl", 35.4532)
argon =          Element( 18, "argon", "Ar", 39.9481)
potassium =      Element( 19, "potassium", "K", 39.09831)
calcium =        Element( 20, "calcium", "Ca", 40.0784)
scandium =       Element( 21, "scandium", "Sc", 44.9559126)
titanium =       Element( 22, "titanium", "Ti", 47.8671)
vanadium =       Element( 23, "vanadium", "V", 50.94151)
chromium =       Element( 24, "chromium", "Cr", 51.99616)
manganese =      Element( 25, "manganese", "Mn", 54.9380455)
iron =           Element( 26, "iron", "Fe", 55.8452)
cobalt =         Element( 27, "cobalt", "Co", 58.9331955)
nickel =         Element( 28, "nickel", "Ni", 58.69342)
copper =         Element( 29, "copper", "Cu", 63.5463)
zinc =           Element( 30, "zinc", "Zn", 65.4094)
gallium =        Element( 31, "gallium", "Ga", 69.7231)
germanium =      Element( 32, "germanium", "Ge", 72.641)
arsenic =        Element( 33, "arsenic", "As", 74.921602)
selenium =       Element( 34, "selenium", "Se", 78.963)
bromine =        Element( 35, "bromine", "Br", 79.9041)
krypton =        Element( 36, "krypton", "Kr", 83.7982)
rubidium =       Element( 37, "rubidium", "Rb", 85.46783)
strontium =      Element( 38, "strontium", "Sr", 87.621)
yttrium =        Element( 39, "yttrium", "Y", 88.905852)
zirconium =      Element( 40, "zirconium", "Zr", 91.2242)
niobium =        Element( 41, "niobium", "Nb", 92.906382)
molybdenum =     Element( 42, "molybdenum", "Mo", 95.942)
technetium =     Element( 43, "technetium", "Tc", 98)
ruthenium =      Element( 44, "ruthenium", "Ru", 101.072)
rhodium =        Element( 45, "rhodium", "Rh", 102.905502)
palladium =      Element( 46, "palladium", "Pd", 106.421)
silver =         Element( 47, "silver", "Ag", 107.86822)
cadmium =        Element( 48, "cadmium", "Cd", 112.4118)
indium =         Element( 49, "indium", "In", 114.8183)
tin =            Element( 50, "tin", "Sn", 118.7107)
antimony =       Element( 51, "antimony", "Sb", 121.7601)
tellurium =      Element( 52, "tellurium", "Te", 127.603)
iodine =         Element( 53, "iodine", "I", 126.904473)
xenon =          Element( 54, "xenon", "Xe", 131.2936)
cesium =         Element( 55, "cesium", "Cs", 132.90545192)
barium =         Element( 56, "barium", "Ba", 137.3277)
lanthanum =      Element( 57, "lanthanum", "La", 138.905477)
cerium =         Element( 58, "cerium", "Ce", 140.1161)
praseodymium =   Element( 59, "praseodymium", "Pr", 140.907652)
neodymium =      Element( 60, "neodymium", "Nd", 144.2423)
promethium =     Element( 61, "promethium", "Pm", 145)
samarium =       Element( 62, "samarium", "Sm", 150.362)
europium =       Element( 63, "europium", "Eu", 151.9641)
gadolinium =     Element( 64, "gadolinium", "Gd", 157.253)
terbium =        Element( 65, "terbium", "Tb", 158.925352)
dysprosium =     Element( 66, "dysprosium", "Dy", 162.5001)
holmium =        Element( 67, "holmium", "Ho", 164.930322)
erbium =         Element( 68, "erbium", "Er", 167.2593)
thulium =        Element( 69, "thulium", "Tm", 168.934212)
ytterbium =      Element( 70, "ytterbium", "Yb", 173.043)
lutetium =       Element( 71, "lutetium", "Lu", 174.9671)
hafnium =        Element( 72, "hafnium", "Hf", 178.492)
tantalum =       Element( 73, "tantalum", "Ta", 180.947882)
tungsten =       Element( 74, "tungsten", "W", 183.841)
rhenium =        Element( 75, "rhenium", "Re", 186.2071)
osmium =         Element( 76, "osmium", "Os", 190.233)
iridium =        Element( 77, "iridium", "Ir", 192.2173)
platinum =       Element( 78, "platinum", "Pt", 195.0849)
gold =           Element( 79, "gold", "Au", 196.9665694)
mercury =        Element( 80, "mercury", "Hg", 200.592)
thallium =       Element( 81, "thallium", "Tl", 204.38332)
lead =           Element( 82, "lead", "Pb", 207.21)
bismuth =        Element( 83, "bismuth", "Bi", 208.980401)
polonium =       Element( 84, "polonium", "Po", 209)
astatine =       Element( 85, "astatine", "At", 210)
radon =          Element( 86, "radon", "Rn", 222.018)
francium =       Element( 87, "francium", "Fr", 223)
radium =         Element( 88, "radium", "Ra", 226)
actinium =       Element( 89, "actinium", "Ac", 227)
thorium =        Element( 90, "thorium", "Th", 232.038062)
protactinium =   Element( 91, "protactinium", "Pa", 231.035882)
uranium =        Element( 92, "uranium", "U", 238.028913)
neptunium =      Element( 93, "neptunium", "Np", 237)
plutonium =      Element( 94, "plutonium", "Pu", 244)
americium =      Element( 95, "americium", "Am", 243)
curium =         Element( 96, "curium", "Cm", 247)
berkelium =      Element( 97, "berkelium", "Bk", 247)
californium =    Element( 98, "californium", "Cf", 251)
einsteinium =    Element( 99, "einsteinium", "Es", 252)
fermium =        Element(100, "fermium", "Fm", 257)
mendelevium =    Element(101, "mendelevium", "Md", 258)
nobelium =       Element(102, "nobelium", "No", 259)
lawrencium =     Element(103, "lawrencium",     "Lr", 262)
rutherfordium =  Element(104, "rutherfordium",  "Rf", 261)
dubnium =        Element(105, "dubnium",        "Db", 262)
seaborgium =     Element(106, "seaborgium",     "Sg", 266)
bohrium =        Element(107, "bohrium",        "Bh", 264)
hassium =        Element(108, "hassium",        "Hs", 269)
meitnerium =     Element(109, "meitnerium",     "Mt", 268)
darmstadtium =   Element(110, "darmstadtium",   "Ds", 281)
roentgenium =    Element(111, "roentgenium",    "Rg", 272)
ununbium =       Element(112, "ununbium",       "Uub", 285)
ununtrium =      Element(113, "ununtrium",      "Uut", 284)
ununquadium =    Element(114, "ununquadium",    "Uuq", 289)
ununpentium =    Element(115, "ununpentium",    "Uup", 288)
ununhexium =     Element(116, "ununhexium",     "Uuh", 292)

# Aliases to recognize common alternative spellings. Both the '==' and 'is'
# relational operators will work with any chosen name
sulphur = sulfur
aluminium = aluminum
