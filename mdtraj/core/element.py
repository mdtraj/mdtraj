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
from mdtraj.utils.unit.quantity import is_quantity
from mdtraj.utils.unit.unit_definitions import daltons

class Element(tuple):
    """An Element represents a chemical element.

    The mdtraj.pdb.element module contains objects for all the standard chemical elements,
    such as element.hydrogen or element.carbon.  You can also call the static method
    Element.getBySymbol() to look up the Element with a particular chemical symbol."""
    __slots__ = []
    _elements_by_symbol = {}
    _elements_by_atomic_number = {}

    def __new__(cls, number, name, symbol, mass, radius):
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
        radius : float
            The van der Waals radius of the element, in nm.
        """

        newobj = tuple.__new__(cls, (number, name, symbol, mass, radius))

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

    def __reduce__(self):
        # __reduce__ is part of the pickle protocol. we need to make sure that
        # elements still act as singletons after a pickle load/save cycle --
        # so load() has to *not create* a new object.
        # see http://docs.python.org/3.3/library/pickle.html#object.__reduce__

        # relevant test:
        # >>> cPickle.loads(cPickle.dumps(md.load(get_fn('bpti.pdb')).topology))
        return str(self.name)

    @staticmethod
    def getBySymbol(symbol):
        """Get the Element with a particular chemical symbol

        Parameters
        ----------
        symbol : str

        Returns
        -------
        element : Element
        """
        s = symbol.strip().upper()
        return Element._elements_by_symbol[s]

    @staticmethod
    def getByAtomicNumber(number):
        """ Get the element with a particular atomic number

        Parameters
        ----------
        number : int

        Returns
        -------
        element : Element
        """
        return Element._elements_by_atomic_number[number]

    @staticmethod
    def getByMass(mass):
        """Get the element whose mass is CLOSEST to the requested mass. This
        method should not be used for repartitioned masses.

        Parameters
        ----------
        mass : float

        Returns
        -------
        element : Element
        """
        # Convert any masses to daltons
        if is_quantity(mass):
            mass = mass.value_in_unit(daltons)
        diff = mass
        best_guess = None

        for key in Element._elements_by_atomic_number:
            element = Element._elements_by_atomic_number[key]
            massdiff = abs(element.mass - mass)
            if massdiff < diff:
                best_guess = element
                diff = massdiff

        return best_guess

    @property
    def number(self):
        """Atomic number."""
        return tuple.__getitem__(self, 0)

    @property
    def name(self):
        """Element name"""
        return tuple.__getitem__(self, 1)

    @property
    def symbol(self):
        """Element symbol"""
        return tuple.__getitem__(self, 2)

    @property
    def mass(self):
        """Element mass"""
        return tuple.__getitem__(self, 3)

    @property
    def radius(self):
        """Element atomic radius

        van der Waals radii are taken from A. Bondi, J. Phys. Chem., 68, 441 -
        452, 1964, except the value for H, which is taken from R.S.
        Rowland & R. Taylor, J.Phys.Chem., 100, 7384 - 7391, 1996. Radii
        that are not available in either of these publications have RvdW =
        2.00 A. The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on
        the CHARMM27 Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by
        default.
        """
        return tuple.__getitem__(self, 4)

    def __getitem__(self, item):
        raise TypeError

    def __str__(self):
        return self.name

    @property
    def atomic_number(self):
        """Atomic number"""
        return tuple.__getitem__(self, 0)

    # Make it so only virtual sites evaluate to boolean False (since it's really
    # *not* an element)
    def __bool__(self):
        return bool(self.mass)

    def __nonzero__(self):
        return bool(self.mass)



# This is for backward compatibility.
def get_by_symbol(symbol):
    s = symbol.strip().upper()
    return Element._elements_by_symbol[s]



# van der Waals radii are taken from A. Bondi,
# J. Phys. Chem., 68, 441 - 452, 1964,
# except the value for H, which is taken from R.S. Rowland & R. Taylor,
# J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
# either of these publications have RvdW = 2.00 A
# The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
# Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.

virtual =        Element(  0,"virtual_site","VS", 0.0, 0.0)
hydrogen =       Element(  1,"hydrogen","H", 1.007947, 0.12)
deuterium =      Element(  1,"deuterium","D", 2.0135532127, 0.12)
helium =         Element(  2,"helium","He", 4.003, 0.14)
lithium =        Element(  3,"lithium","Li", 6.9412, 0.182)
beryllium =      Element(  4,"beryllium","Be", 9.0121823, 0.2)
boron =          Element(  5,"boron","B", 10.8117, 0.2)
carbon =         Element(  6,"carbon","C", 12.01078, 0.17)
nitrogen =       Element(  7,"nitrogen","N", 14.00672, 0.155)
oxygen =         Element(  8,"oxygen","O", 15.99943, 0.152)
fluorine =       Element(  9,"fluorine","F", 18.99840325, 0.147)
neon =           Element( 10,"neon","Ne", 20.17976, 0.154)
sodium =         Element( 11,"sodium","Na", 22.989769282, 0.136)
magnesium =      Element( 12,"magnesium","Mg", 24.30506, 0.118)
aluminum =       Element( 13,"aluminum","Al", 26.98153868, 0.2)
silicon =        Element( 14,"silicon","Si", 28.08553, 0.21)
phosphorus =     Element( 15,"phosphorus","P", 30.9737622, 0.18)
sulfur =         Element( 16,"sulfur","S", 32.0655, 0.18)
chlorine =       Element( 17,"chlorine","Cl", 35.4532, 0.227)
argon =          Element( 18,"argon","Ar", 39.9481, 0.188)
potassium =      Element( 19,"potassium","K", 39.09831, 0.176)
calcium =        Element( 20,"calcium","Ca", 40.0784, 0.137)
scandium =       Element( 21,"scandium","Sc", 44.9559126, 0.2)
titanium =       Element( 22,"titanium","Ti", 47.8671, 0.2)
vanadium =       Element( 23,"vanadium","V", 50.94151, 0.2)
chromium =       Element( 24,"chromium","Cr", 51.99616, 0.2)
manganese =      Element( 25,"manganese","Mn", 54.9380455, 0.2)
iron =           Element( 26,"iron","Fe", 55.8452, 0.2)
cobalt =         Element( 27,"cobalt","Co", 58.9331955, 0.2)
nickel =         Element( 28,"nickel","Ni", 58.69342, 0.163)
copper =         Element( 29,"copper","Cu", 63.5463, 0.14)
zinc =           Element( 30,"zinc","Zn", 65.4094, 0.139)
gallium =        Element( 31,"gallium","Ga", 69.7231, 0.107)
germanium =      Element( 32,"germanium","Ge", 72.641, 0.2)
arsenic =        Element( 33,"arsenic","As", 74.921602, 0.185)
selenium =       Element( 34,"selenium","Se", 78.963, 0.19)
bromine =        Element( 35,"bromine","Br", 79.9041, 0.185)
krypton =        Element( 36,"krypton","Kr", 83.7982, 0.202)
rubidium =       Element( 37,"rubidium","Rb", 85.46783, 0.2)
strontium =      Element( 38,"strontium","Sr", 87.621, 0.2)
yttrium =        Element( 39,"yttrium","Y", 88.905852, 0.2)
zirconium =      Element( 40,"zirconium","Zr", 91.2242, 0.2)
niobium =        Element( 41,"niobium","Nb", 92.906382, 0.2)
molybdenum =     Element( 42,"molybdenum","Mo", 95.942, 0.2)
technetium =     Element( 43,"technetium","Tc", 98, 0.2)
ruthenium =      Element( 44,"ruthenium","Ru", 101.072, 0.2)
rhodium =        Element( 45,"rhodium","Rh", 102.905502, 0.2)
palladium =      Element( 46,"palladium","Pd", 106.421, 0.163)
silver =         Element( 47,"silver","Ag", 107.86822, 0.172)
cadmium =        Element( 48,"cadmium","Cd", 112.4118, 0.158)
indium =         Element( 49,"indium","In", 114.8183, 0.193)
tin =            Element( 50,"tin","Sn", 118.7107, 0.217)
antimony =       Element( 51,"antimony","Sb", 121.7601, 0.2)
tellurium =      Element( 52,"tellurium","Te", 127.603, 0.206)
iodine =         Element( 53,"iodine","I", 126.904473, 0.198)
xenon =          Element( 54,"xenon","Xe", 131.2936, 0.216)
cesium =         Element( 55,"cesium","Cs", 132.90545192, 0.21)
barium =         Element( 56,"barium","Ba", 137.3277, 0.2)
lanthanum =      Element( 57,"lanthanum","La", 138.905477, 0.2)
cerium =         Element( 58,"cerium","Ce", 140.1161, 0.2)
praseodymium =   Element( 59,"praseodymium","Pr", 140.907652, 0.2)
neodymium =      Element( 60,"neodymium","Nd", 144.2423, 0.2)
promethium =     Element( 61,"promethium","Pm", 145, 0.2)
samarium =       Element( 62,"samarium","Sm", 150.362, 0.2)
europium =       Element( 63,"europium","Eu", 151.9641, 0.2)
gadolinium =     Element( 64,"gadolinium","Gd", 157.253, 0.2)
terbium =        Element( 65,"terbium","Tb", 158.925352, 0.2)
dysprosium =     Element( 66,"dysprosium","Dy", 162.5001, 0.2)
holmium =        Element( 67,"holmium","Ho", 164.930322, 0.2)
erbium =         Element( 68,"erbium","Er", 167.2593, 0.2)
thulium =        Element( 69,"thulium","Tm", 168.934212, 0.2)
ytterbium =      Element( 70,"ytterbium","Yb", 173.043, 0.2)
lutetium =       Element( 71,"lutetium","Lu", 174.9671, 0.2)
hafnium =        Element( 72,"hafnium","Hf", 178.492, 0.2)
tantalum =       Element( 73,"tantalum","Ta", 180.947882, 0.2)
tungsten =       Element( 74,"tungsten","W", 183.841, 0.2)
rhenium =        Element( 75,"rhenium","Re", 186.2071, 0.2)
osmium =         Element( 76,"osmium","Os", 190.233, 0.2)
iridium =        Element( 77,"iridium","Ir", 192.2173, 0.2)
platinum =       Element( 78,"platinum","Pt", 195.0849, 0.172)
gold =           Element( 79,"gold","Au", 196.9665694, 0.166)
mercury =        Element( 80,"mercury","Hg", 200.592, 0.155)
thallium =       Element( 81,"thallium","Tl", 204.38332, 0.196)
lead =           Element( 82,"lead","Pb", 207.21, 0.202)
bismuth =        Element( 83,"bismuth","Bi", 208.980401, 0.2)
polonium =       Element( 84,"polonium","Po", 209, 0.2)
astatine =       Element( 85,"astatine","At", 210, 0.2)
radon =          Element( 86,"radon","Rn", 222.018, 0.2)
francium =       Element( 87,"francium","Fr", 223, 0.2)
radium =         Element( 88,"radium","Ra", 226, 0.2)
actinium =       Element( 89,"actinium","Ac", 227, 0.2)
thorium =        Element( 90,"thorium","Th", 232.038062, 0.2)
protactinium =   Element( 91,"protactinium","Pa", 231.035882, 0.2)
uranium =        Element( 92,"uranium","U", 238.028913, 0.186)
neptunium =      Element( 93,"neptunium","Np", 237, 0.2)
plutonium =      Element( 94,"plutonium","Pu", 244, 0.2)
americium =      Element( 95,"americium","Am", 243, 0.2)
curium =         Element( 96,"curium","Cm", 247, 0.2)
berkelium =      Element( 97,"berkelium","Bk", 247, 0.2)
californium =    Element( 98,"californium","Cf", 251, 0.2)
einsteinium =    Element( 99,"einsteinium","Es", 252, 0.2)
fermium =        Element(100,"fermium","Fm", 257, 0.2)
mendelevium =    Element(101,"mendelevium","Md", 258, 0.2)
nobelium =       Element(102,"nobelium","No", 259, 0.2)
lawrencium =     Element(103,"lawrencium","Lr", 262, 0.2)
rutherfordium =  Element(104,"rutherfordium","Rf", 261, 0.2)
dubnium =        Element(105,"dubnium","Db", 262, 0.2)
seaborgium =     Element(106,"seaborgium","Sg", 266, 0.2)
bohrium =        Element(107,"bohrium","Bh", 264, 0.2)
hassium =        Element(108,"hassium","Hs", 269, 0.2)
meitnerium =     Element(109,"meitnerium","Mt", 268, 0.2)
darmstadtium =   Element(110,"darmstadtium","Ds", 281, 0.2)
roentgenium =    Element(111,"roentgenium","Rg", 272, 0.2)
ununbium =       Element(112,"ununbium","Uub", 285, 0.2)
ununtrium =      Element(113,"ununtrium","Uut", 284, 0.2)
ununquadium =    Element(114,"ununquadium","Uuq", 289, 0.2)
ununpentium =    Element(115,"ununpentium","Uup", 288, 0.2)
ununhexium =     Element(116,"ununhexium","Uuh", 292, 0.2)

# Aliases to recognize common alternative spellings. Both the '==' and 'is'
# relational operators will work with any chosen name
sulphur = sulfur
aluminium = aluminum
virtual_site = virtual
