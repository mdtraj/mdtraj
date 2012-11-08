#!/usr/bin/env python
# This file is part of MSMBuilder.
#
# Toni Giorgino at Upf dot Edu
# Copyright 2011 Universitat Pompeu Fabra
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
A simple interface to VMD's DCD molfile plugin.
"""

# July 2011: Initial version (Toni)

import os, sys
import numpy as np
from ctypes import Structure, POINTER, c_float, c_double, CDLL, c_char_p, c_int
from ctypes import c_void_p, byref
from ctypes.util import find_library
import imp
import warnings

# define handle to dcd library as global variable (but it should only be used within this module)
_dcdlib = None


class MolfileTimestep(Structure):
    """Wrapper for the timestep C structure used in molfile plugins."""
    _fields_ = [
        ("coords", POINTER(c_float)),
        ("velocities", POINTER(c_float)),
        ("A", c_float),
        ("B", c_float),
        ("C", c_float),
        ("alpha", c_float),
        ("beta", c_float),
        ("gamma", c_float),
        ("physical_time", c_double)
        ]


def loadDCDLibrary(LoadDirectFromPackage=True):
    global _dcdlib

    dcd_library_path = find_library("_dcd.so")
    if dcd_library_path and not LoadDirectFromPackage:
        _dcdlib = CDLL(dcd_library_path)
    elif LoadDirectFromPackage == True:
        PackagePath = imp.find_module("mdtraj")[1]
        _dcdlib = CDLL(os.path.join(PackagePath, '_dcd.so'))
    else:
        # Lutz: find_library does not look in LD_LIBRARY_PATH on linux
        # machines (see http://bugs.python.org/issue2936), even though
        # cdll.LoadLibrary does as a workaround, we try to load the
        # shared object file manually
        if sys.platform.startswith("linux"):
            try:
                _dcdlib = CDLL("_dcd.so")
            except:
                raise RuntimeError("Unable to find dcdplugin_s library. Make sure that it is compiled as a shared library, and that its location is known to the dynamic linker (e.g., set LD_LIBRARY_PATH on Linux).")
        else:
            raise RuntimeError("Unable to find libdcdfile.so library. Make sure that it is compiled as a shared library, and that its location is known to the dynamic linker (e.g., set DYLD_LIBRARY_PATH on Mac OS X, LD_LIBRARY_PATH on Linux).")

    # declare interface for functions in the dcd library to insure
    # some amount of type safety for use of numpy arrays in ctypes,
    # (DCD writing functions could also be wrapped)

    _dcdlib.open_dcd_read.argtypes = [c_char_p, c_char_p, POINTER(c_int)]
    _dcdlib.open_dcd_read.restype = c_void_p
    
    # 0 is OK,  -1 is EOF, else error
    _dcdlib.read_next_timestep.argtypes = [c_void_p, c_int, POINTER(MolfileTimestep)]
    _dcdlib.read_next_timestep.restype = c_int

    _dcdlib.close_file_read.argtypes = [c_void_p]
    _dcdlib.close_file_read.restype = None


class DCDReader:
    """Object that allows iteration over the configurations in a dcd trajectory file.

    Note that the returned configuration *might* contain references to an instance-owned memory block
    that contains the atom coordinates and box geometries. You have to copy them if you need
    this data beyond the next iteration over this class.
    
    Examples:

    for c in DCDReader("filename.dcd"):
        print c.step

    for c in DCDReader(["part1.dcd","part2.dcd"], firstframe=100, lastframe = 200, stepframe = 2, atomindices = [0,1,2]):
        print c.step
    """
        
    def _open(self, filename):
        """Opens the dcd file with the specified name."""
        self._filename = filename
        self.dcd = _dcdlib.open_dcd_read(filename, "dcd", byref(self.natoms))
        if not self.dcd:
            raise IOError("Unable to open dcd file " + str(filename) + ".")
        # print "Opened file %s with %d atoms" % (self._filename,self.natoms.value)
        if self._atomindices != None:
            print "Keeping %d atoms." % len(self._atomindices)

    def _close(self):
        """Closes currently open dcd file."""
        if hasattr(self, "dcd") and self.dcd:
            _dcdlib.close_file_read(self.dcd)
            self.dcd = None

    def __init__(self, filenames, firstframe=0, lastframe=None,
                 stepframe=1, atomindices=None, skipcont=True):
        self._firstframe = firstframe
        self._lastframe = lastframe
        self._stepframe = stepframe
        self._nextframe = self._firstframe
        self._atomindices = atomindices
        self._skipcont = skipcont
        
        if isinstance(filenames, str):
            self._filenames = [filenames]
        else:
            self._filenames = filenames[:]  # Want to copy the list so as not to destroy it.

        print(self._filenames)

        # Init the plugin
        if _dcdlib.vmdplugin_init() != 0:
            raise IOError("Unable to init DCD plugin")

        self._ts = MolfileTimestep()
        self.natoms = c_int(-1)
        self._frame = -1         # why was -1?

        # self.natoms = number_of_atoms(self._filenames[0])
#         self._allcoords = np.empty([self.natoms,3],dtype='single',order='C')
#         if self._atomindices != None:
#             self._coords = np.empty([len(self._atomindices),3],dtype='single',order='C')
#         self._box = np.empty([3,3],dtype='single',order='C')
#         self._step = c_int()
#         self._time = c_float()
#         self._precision = c_float()

        self._open(self._filenames.pop(0))

    def __del__(self):
        self._close()
        
    def __iter__(self):
        return self

    def next(self):
        Coords = c_float * (self.natoms.value * 3)
        xyzvec = Coords()
        self._ts.coords = xyzvec

        while self._frame < self._nextframe:
            if self._lastframe != None and self._frame >= self._lastframe:
                raise StopIteration
#            result = _dcdlib.read_dcd(self.dcd, self.natoms, byref(self._step), byref(self._time), self._box, self._allcoords, byref(self._precision))
            result = _dcdlib.read_next_timestep(self.dcd, self.natoms, byref(self._ts))
            if result == 0:
                self._frame += 1
            elif result == -1:  # EOF
                print "Finished with file %s, %d atoms, at frame %d" % \
                    (self._filename, self.natoms.value, self._frame + 1)
                self._close()
                if not self._filenames:
                    raise StopIteration
                else:
                    self._open(self._filenames.pop(0))
                    if self._skipcont:
                        self._nextframe += 1
            else:
                raise IOError("An error occured while reading dcd file.")

        # now generate the structure that will be returned to caller
           # config=Configuration(self._step.value,self._time.value,
           #       self._precision.value,self._box,self._allcoords,self._atomindices)

        self._nextframe += self._stepframe

        # create a "coords" numpy array for returning
        with warnings.catch_warnings():
            # supress the PEP3118 warning
            # http://stackoverflow.com/questions/4964101/pep-3118-warning-when-using-ctypes-array-as-numpy-array
            warnings.simplefilter("ignore")
            coords = np.asfarray(np.array(xyzvec).reshape(self.natoms.value, 3))
            
        if self._atomindices != None:
            coords = coords[self._atomindices, ]
        coords = coords * 0.1       # \AA -> nm
#        print coords
        return coords

try:
    loadDCDLibrary()
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
    if _dcdlib:
        print "Successfully loaded dcdfile library."
    else:
        print "Unable to load dcdfile library."
