#cython: c_string_type=str, c_string_encoding=ascii
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

import cython
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

# codes that indicate status on return from library
cdef int _BINPOS_SUCESS = 0  # regular exit code
cdef int _BINPOS_EOF = -1  # end of file (or error)

def read(filename, chunk=1):
    """Read the data from a AMBER binpos file

    Parameters
    ----------
    filename : str
        The filename of the binpos file to read from
    chunk : int
        Size of the chunks to read

    Returns
    -------
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        The xyz coordinates of each atom in each frame
    """
    xyz = np.vstack(tuple(BINPOSReader(filename, chunk)))

    return xyz


def write(filename, xyz, force_overwrite=True):
    """Write xyz coordinates to a AMBER binpos file

    Note that the box size entries in the BINPOS file will be left blank (zeros)

    Parameters
    ----------
    filename : str
        The path to the binpos file to write to
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        The xyz coordinates of each atom in each frame
    force_overwrite : bool, default=False
        Overwrite anything that exists at filename, if its already there
    """

    if force_overwrite and os.path.exists(filename):
        os.unlink(filename)

    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    #make sure all the arrays are the right shape
    xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False)

    writer = BINPOSWriter(filename, xyz)
    writer.write()


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
