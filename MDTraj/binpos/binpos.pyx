# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

###############################################################################
# Imports
###############################################################################

import cython
import warnings
cimport cython
import os
import numpy as np
cimport numpy as np
np.import_array()
from mdtraj.utils.arrays import ensure_type
from libc.stdlib cimport malloc, free
from binposlib cimport molfile_timestep_t
from binposlib cimport open_binpos_read, close_file_read, read_next_timestep
from binposlib cimport open_binpos_write, close_file_write, write_timestep

###############################################################################
# Globals
###############################################################################

# codes that indicate status on return from library
cdef int _BINPOS_SUCESS = 0  # regular exit code
cdef int _BINPOS_EOF = -1  # end of file (or error)

###############################################################################
# Classes
###############################################################################

cdef class BINPOSTrajectoryFile:
    cdef int n_atoms
    cdef int approx_n_frames
    cdef int min_chunk_size
    cdef int chunk_size_multiplier
    cdef int is_open
    cdef int frame_counter

    cdef void* fh
    cdef molfile_timestep_t* timestep

    def __cinit__(self, char* filename, char* mode=b'r', force_overwrite=True):
        """Open an AMBER BINPOS file for reading/writing.

        Parameters
        ----------
        filename : str
            The filename to open. A path to a file on disk.
        mode : {'r', 'w'}
            The mode in which to open the file, either 'r' for read or 'w' for write.
        force_overwrite : bool
            If opened in write mode, and a file by the name of `filename` already exists on disk, should we overwrite it?

        Other Parameters
        ----------------
        min_chunk_size : int, default=100
            In read mode, we need to allocate a buffer in which to store the data without knowing how many frames are
            in the file. This parameter is the minimum size of the buffer to allocate.
        chunk_size_multiplier, int, default=1.5
            In read mode, we need to allocate a buffer in which to store the data without knowing how many frames are in
            the file. We can *guess* this information based on the size of the file on disk, but it's not perfect. This
            parameter inflates the guess by a multiplicative factor.
        """
        self.is_open = False
        self.frame_counter = 0

        if mode == b'r':
            self.n_atoms = 0
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self.fh = open_binpos_read(filename, "binpos", &self.n_atoms)
            if self.n_atoms <= 0:
                raise IOError('Malformed BINPOS file. Number of atoms <= 0. '
                              'Are you sure this is a valid AMBER BINPOS file?')

            # binpos stores the 8 bytes per float, 3*n_atoms floats per frame, directly
            # with no compression
            self.approx_n_frames = os.stat(filename).st_size / (self.n_atoms * 8 * 3)

            self.min_chunk_size = kwargs.pop('min_chunk_size', 100)
            self.chunk_size_multiplier = kwargs.pop('chunk_size_multiplier', 1.5)

        elif mode == b'w':
            pass

        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.timestep is NULL:
            raise MemoryError('There was an error allocating working space')

        if self.fh is NULL:
            raise IOError('There was an error opening the binpos file: %s' % filename)

        for key in kwargs.keys():
            warnings.warn('kwarg "%s" was not recognized or processed' % key)

        self.is_open = True
        self.mode = mode

    def close(self):
        if self.is_open:
            close_file_read(self.fh)
            self.is_open = False

    def __dealloc__(self):
        free(self.timestep)
        self.close()

    def read(self, n_frames=None):
        """Read data from a BINPOS file

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
            The cartesian coordinates
        """
        if not self.mode == b'r':
            raise ValueError('read() is only available when file is opened in mode="r"')

    def _read(self, n_frames):
        pass

    def write(self, xyz):
        pass




cdef class BINPOSReader:
    cdef int n_atoms
    cdef void* fh
    cdef molfile_timestep_t* timestep
    cdef int status
    cdef int chunk

    def __cinit__(self, char* filename, int chunk=1):
        self.chunk = chunk
        #open the file
        self.fh = open_binpos_read(filename, "binpos", &self.n_atoms)
        if self.fh is NULL:
            raise IOError('There was an error opening the binpos file: %s' % filename)

        # alloc the molfile_timestep, which is the struct that the library is
        # going to deposit its data into each timestep
        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.fh is  NULL:
            raise MemoryError


    def __dealloc__(self):
        "Shut this whole thing down"

        # free whatever we malloced
        free(self.timestep)

        # close the file
        if self.fh is not NULL:
            close_file_read(self.fh)


    def __iter__(self):
        return self

    @cython.boundscheck(False)
    def __next__(self):
        cdef int i
        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((self.chunk, self.n_atoms, 3), dtype=np.float32)

        for i in range(self.chunk):
            self.timestep.coords = &xyz[i,0,0]
            status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            if status != _BINPOS_SUCESS:
                if i == 0:
                    raise StopIteration
                break

        if status != _BINPOS_SUCESS:
            return xyz[0:i]

        return xyz


cdef class BINPOSWriter:
    cdef char* filename
    cdef void* fh
    cdef np.ndarray xyz
    cdef int n_atoms
    cdef molfile_timestep_t* timestep

    def __cinit__(self, char* filename, np.ndarray[np.float32_t, ndim=3, mode="c"] xyz):
        """
        Set up the BINPOS writer

        Nothing will be written to disk until you call write()

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
            The xyz coordinates of each atom in each frame
        """
        self.filename = filename
        self.xyz = xyz
        self.n_atoms = self.xyz.shape[1]

        self.fh = open_binpos_write(filename, "binpos", self.n_atoms)
        if self.fh is NULL:
            raise IOError('There was an error opening the file: %s' % filename)

        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.timestep is NULL:
            raise MemoryError

    def __dealloc__(self):
        close_file_write(self.fh)
        if self.timestep is not NULL:
            free(self.timestep)

    @cython.boundscheck(False)
    def write(self):
        "Write all the data"
        cdef int i, status
        cdef n_frames = len(self.xyz)
        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = self.xyz

        for i in range(n_frames):
            self.timestep.coords = &xyz[i, 0, 0]
            status = write_timestep(self.fh, self.timestep)

            if status != _BINPOS_SUCESS:
                raise RuntimeError("BINPOS Error: %s" % status)
