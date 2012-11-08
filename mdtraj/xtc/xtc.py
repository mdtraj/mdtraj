#!/usr/bin/env python
# Copyright 2011 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
A simple interface to the Gromacs XDRLIB library.
"""

# June 2010: Initial version (Lutz)

import numpy as np
from numpy.ctypeslib import ndpointer
from ctypes import *
from ctypes.util import find_library
import os.path
import sys
import imp

# define handle to xdr library as global variable (but it should only be used within this module)
_xdrlib = None

# return codes of xdr library (defined in xdrfile.h, descriptions in xdrfile.c)
_EXDROK = 0             # OK
_EXDRHEADER = 1         # Header
_EXDRSTRING = 2         # String
_EXDRDOUBLE = 3         # Double
_EXDRINT = 4            # Integer
_EXDRFLOAT = 5          # Float
_EXDRUINT = 6           # Unsigned integer
_EXDR3DX = 7            # Compressed 3d coordinate
_EXDRCLOSE = 8          # Closing file
_EXDRMAGIC = 9          # Magic number
_EXDRNOMEM = 10         # Not enough memory
_EXDRENDOFFILE = 11     # End of file
_EXDRFILENOTFOUND = 12  # File not found


def loadXDRLibrary(LoadDirectFromPackage=True):
    global _xdrlib

    xdr_library_path = find_library("xdrfile")
    if xdr_library_path and not LoadDirectFromPackage:
        _xdrlib = CDLL(xdr_library_path)
    elif LoadDirectFromPackage:
        PackagePath=imp.find_module("mdtraj")[1]
        _xdrlib=CDLL(os.path.join(PackagePath, "_xtc.so"))
    else:
        # Lutz: find_library does not look in LD_LIBRARY_PATH on linux machines (see http://bugs.python.org/issue2936), even though cdll.LoadLibrary does
        # as a workaround, we try to load the shared object file manually
        if sys.platform.startswith("linux"):
            try:
                _xdrlib = CDLL("_xtc.so")
            except:
                raise RuntimeError("Unable to find xdrfile library. Make sure that it is compiled as a shared library, and that its location is known to the dynamic linker (e.g., set LD_LIBRARY_PATH on Linux).")
        else:
            raise RuntimeError("Unable to find xdrfile library. Make sure that it is compiled as a shared library, and that its location is known to the dynamic linker (e.g., set DYLD_LIBRARY_PATH on Mac OS X, LD_LIBRARY_PATH on Linux).")

    # declare interface for functions in the xdr library to insure some amount of type safety
    # for use of numpy arrays in ctypes, see http://thread.gmane.org/gmane.comp.python.numeric.general/7418/focus=7418
    _xdrlib.read_xtc_natoms.argtypes = [c_char_p, POINTER(c_int)]
    _xdrlib.read_trr_natoms.argtypes = [c_char_p, POINTER(c_int)]
    _xdrlib.xdrfile_open.argtypes = [c_char_p, c_char_p]
    _xdrlib.xdrfile_open.restype = c_void_p
    _xdrlib.xdrfile_close.argtypes = [c_void_p]
    _xdrlib.read_xtc.argtypes = [c_void_p, c_int, POINTER(c_int), POINTER(c_float), ndpointer(dtype="single",shape=(3,3),flags="C_CONTIGUOUS"), ndpointer(dtype="single",ndim=2,flags="C_CONTIGUOUS"), POINTER(c_float)]
    _xdrlib.read_trr.argtypes = [c_void_p, c_int, POINTER(c_int), POINTER(c_float),POINTER(c_float), ndpointer(dtype="single",shape=(3,3),flags="C_CONTIGUOUS"), ndpointer(dtype="single",ndim=2,flags="C_CONTIGUOUS"),ndpointer(dtype="single",ndim=2,flags="C_CONTIGUOUS"),ndpointer(dtype="single",ndim=2,flags="C_CONTIGUOUS")]
    _xdrlib.write_xtc.argtypes = [c_void_p, c_int, c_int, c_float, ndpointer(dtype="single",shape=(3,3),flags="C_CONTIGUOUS"), ndpointer(dtype="single",ndim=2,flags="C_CONTIGUOUS"), c_float]

    
def number_of_atoms(filename):
    """Returns the number of atoms in specified xtc file."""
    n = c_int()
    if _xdrlib.read_xtc_natoms(filename, byref(n)) != _EXDROK:
        raise IOError("Unable to determine number of atoms in file " + str(filename) + ".")
    return int(n.value)

def number_of_atoms_trr(filename):
    """Returns the number of atoms in specified xtc file."""
    n = c_int()
    if _xdrlib.read_trr_natoms(filename, byref(n)) != _EXDROK:
        raise IOError("Unable to determine number of atoms in file " + str(filename) + ".")
    return int(n.value)

class XTCWriter:
    def __init__(self, filename, overwrite = False):
        filename = str(filename)
        if os.path.exists(filename) and not overwrite:
            raise IOError("File \"" + filename + "\" already exists.")
        
        self.xdr = _xdrlib.xdrfile_open(filename, "w")
        if not self.xdr:
            raise IOError("Unable to open xtc file " + str(filename) + ".")

    def __del__(self):
        if hasattr(self, "xdr") and self.xdr:
            if _xdrlib.xdrfile_close(self.xdr) != 0:
                raise IOError("Unable to close xtc file.")
            else:
                self.xdr = None

    def write(self, coords, step, time, box, precision):
        if self.xdr:
            if _xdrlib.write_xtc(self.xdr, len(coords), step, time, box, coords, precision) != _EXDROK:
                raise IOError("An error occured while writing to xtc file.")
        else:
            raise IOError("Trying to write to xtc file that has not been opened.")

class Configuration:
    """Structure containing a single frame from a xtc file.    
    Elements:
    coords:    array of atom coordinates
    step:      timestep
    time:      time
    precision: precision
    box:       simulation box"""
    def __init__(self,step,time,precision,box,allcoords,atomindices=None):
        self.step=step
        self.time=time
        self.precision=precision
        self.box=box
        if atomindices==None:
            self.coords=allcoords
        else:
            self.coords=allcoords[atomindices]

class XTCReader:
    """Object that allows iteration over the configurations in a xtc trajectory file.

    Note that the returned configuration *might* contain references to an instance-owned memory block
    that contains the atom coordinates and box geometries. You have to copy them if you need
    this data beyond the next iteration over this class.
    
    Examples:

    for c in XTCReader("filename.xtc"):
        print c.step

    for c in XTCReader(["part1.xtc","part2.xtc"], firstframe=100, lastframe = 200, stepframe = 2, atomindices = [0,1,2]):
        print c.step
    """
        


    def _open(self, filename):
        """Opens the xtc file with the specified name."""        
        self.xdr = _xdrlib.xdrfile_open(filename, "r")
        if not self.xdr:
            raise IOError("Unable to open xtc file " + str(filename) + ".")

    def _close(self):
        """Closes currently open xtc file."""
        if hasattr(self, "xdr") and self.xdr:
            if _xdrlib.xdrfile_close(self.xdr) != 0:
                raise IOError("Unable to close xtc file.")
            else:
                self.xdr = None

    def __init__(self, filenames, firstframe = 0, lastframe = None, stepframe = 1, atomindices = None, skipcont = True):
        self._firstframe = firstframe
        self._lastframe = lastframe
        self._stepframe = stepframe
        self._nextframe = self._firstframe
        self._atomindices = atomindices
        self._skipcont = skipcont
        
        if isinstance(filenames, str):
            self._filenames = [filenames]
        else:
            self._filenames = filenames[:]#Want to copy the list so as not to destroy it.

        #print(self._filenames)
        self.natoms = number_of_atoms(self._filenames[0])
        self._allcoords = np.empty([self.natoms,3],dtype='single',order='C')
        if self._atomindices != None:
            self._coords = np.empty([len(self._atomindices),3],dtype='single',order='C')
        self._box = np.empty([3,3],dtype='single',order='C')
        self._step = c_int()
        self._time = c_float()
        self._precision = c_float()
        self._frame = -1

        self._open(self._filenames.pop(0))

    def __del__(self):
        self._close()
        
    def __iter__(self):
        return self

    def next(self):
        while self._frame < self._nextframe:
            if self._lastframe != None and self._frame >= self._lastframe:
                raise StopIteration
            result = _xdrlib.read_xtc(self.xdr, self.natoms, byref(self._step), byref(self._time), self._box, self._allcoords, byref(self._precision))
            if result == _EXDROK:
                self._frame += 1
            elif result == _EXDRENDOFFILE:
                self._close()
                if not self._filenames:
                    raise StopIteration
                else:
                    self._open(self._filenames.pop(0))                   
                    if self._skipcont:
                        self._nextframe += 1
            else:
                raise IOError("An error occured while reading xtc file.")

        # now generate configuration structure that will be returned to caller
        config=Configuration(self._step.value,self._time.value,self._precision.value,self._box,self._allcoords,self._atomindices)
        self._nextframe += self._stepframe
        return config

def readxtc(*arg, **karg):
    """Returns a list of coordinates corresponding to the atom positions in an xtc file.

    Arguments are the same as those of the XTCReader class.
    """
    
    erg = []
    for c in XTCReader(*arg, **karg):
        erg.append(c.coords.copy())
    return erg

class TRRConfiguration:
    """Structure containing a single frame from a xtc file.    
    Elements:
    coords:    array of atom coordinates
    velocities:    array of atom velocities
    forces:    array of atom forces
    step:      timestep
    time:      time
    precision: precision
    box:       simulation box"""
    def __init__(self,step,time,box,allcoords,velocities,forces,atomindices=None):
        self.step=step
        self.time=time
        self.box=box
        if atomindices==None:
            self.coords=allcoords
            self.forces=forces
            self.velocities=velocities
        else:
            self.coords=allcoords[atomindices]
            self.forces=forces[atomindices]
            self.velocities=velocities[atomindices]

class TRRReader:
    """Object that allows iteration over the configurations in a xtc trajectory file.

    Note that the returned configuration *might* contain references to an instance-owned memory block
    that contains the atom coordinates and box geometries. You have to copy them if you need
    this data beyond the next iteration over this class.
    
    Examples:

    for c in TRRReader("filename.trr"):
        print c.step

    for c in TRRReader(["part1.trr","part2.trr"], firstframe=100, lastframe = 200, stepframe = 2, atomindices = [0,1,2]):
        print c.step
    """
        


    def _open(self, filename):
        """Opens the xtc file with the specified name."""        
        self.xdr = _xdrlib.xdrfile_open(filename, "r")
        if not self.xdr:
            raise IOError("Unable to open trr file " + str(filename) + ".")

    def _close(self):
        """Closes currently open xtc file."""
        if hasattr(self, "xdr") and self.xdr:
            if _xdrlib.xdrfile_close(self.xdr) != 0:
                raise IOError("Unable to close trr file.")
            else:
                self.xdr = None

    def __init__(self, filenames, firstframe = 0, lastframe = None, stepframe = 1, atomindices = None, skipcont = True):
        self._firstframe = firstframe
        self._lastframe = lastframe
        self._stepframe = stepframe
        self._nextframe = self._firstframe
        self._atomindices = atomindices
        self._skipcont = skipcont
        
        if isinstance(filenames, str):
            self._filenames = [filenames]
        else:
            self._filenames = filenames[:]#Want to copy the list so as not to destroy it.

        print(self._filenames)
        self.natoms = number_of_atoms_trr(self._filenames[0])
        self._allcoords = np.empty([self.natoms,3],dtype='single',order='C')
        if self._atomindices != None:
            self._coords = np.empty([len(self._atomindices),3],dtype='single',order='C')

        self._allvelocities = np.empty([self.natoms,3],dtype='single',order='C')
        if self._atomindices != None:
            self._velocities = np.empty([len(self._atomindices),3],dtype='single',order='C')

        self._allforces = np.empty([self.natoms,3],dtype='single',order='C')
        if self._atomindices != None:
            self._forces = np.empty([len(self._atomindices),3],dtype='single',order='C')
        
        self._box = np.empty([3,3],dtype='single',order='C')
        self._step = c_int()
        self._time = c_float()
        self._Lambda = c_float()
        self._frame = -1


        self._open(self._filenames.pop(0))

    def __del__(self):
        self._close()
        
    def __iter__(self):
        return self

    def next(self):
        while self._frame < self._nextframe:
            if self._lastframe != None and self._frame >= self._lastframe:
                raise StopIteration
            result = _xdrlib.read_trr(self.xdr, self.natoms, byref(self._step), byref(self._time),byref(self._Lambda), self._box, self._allcoords,self._allvelocities,self._allforces)
            print(result)
            if result == _EXDROK:
                self._frame += 1
            elif result == _EXDRENDOFFILE or result== 4:#HACK by KAB--no idea why 4 is getting returned
                self._close()
                if not self._filenames:
                    raise StopIteration
                else:
                    self._open(self._filenames.pop(0))                   
                    if self._skipcont:
                        self._nextframe += 1
            else:
                raise IOError("An error occured while reading trr file.")

        # now generate configuration structure that will be returned to caller
        config=TRRConfiguration(self._step.value,self._time.value,self._box,self._allcoords,self._allvelocities,self._allforces,self._atomindices)
        self._nextframe += self._stepframe
        return config

try:
    loadXDRLibrary()
except OSError:
    # HACK for building the documentation on READTHEDOCS
    # this is a workaround because DLL open is hard to mock
    # but basically, if we're on the readthedocs.org build server
    # we need to be able to import packages to introspect their docstrings,
    # but we don't actually have any of the C libraries installed
    if os.environ.get('READTHEDOCS', None) == 'True':
        pass
    else:
        raise

if __name__ == "__main__":
    if _xdrlib:
        print "Successfully loaded xdrfile library."
    else:
        print "Unable to load xdrfile library."

    
