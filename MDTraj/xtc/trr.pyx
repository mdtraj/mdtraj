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

import os
import cython
cimport cython
try:
    from itertools import izip
    zip = izip
except ImportError:  # python3
    pass
import numpy as np
cimport numpy as np
import numpy as np
np.import_array()
from mdtraj.utils.arrays import ensure_type
from trrlib cimport (XDRFILE, read_trr_natoms, xdrfile_open, \
    xdrfile_close, write_trr, read_trr)
from collections import namedtuple
TRRFile = namedtuple('TRRFile', ['xyz', 'time', 'step', 'box', 'lambd'])

# numpy variable types include the specific numpy of bytes of each, but the c
# variables in our interface file don't. this could get bad if we're on a wierd
# machine, so lets make sure first
if sizeof(int) != sizeof(np.int32_t):
    raise RuntimeError('Integers on your compiler are not 32 bits. This is not good.')
if sizeof(float) != sizeof(np.float32_t):
    raise RuntimeError('Floats on your compiler are not 32 bits. This is not good')


def read(filename, chunk=1000):
    """
    Read the data from a Gromacs TRR file

    Parameters
    ----------
    filename : str
        The filename of the xtc file to read from
    chunk : int
        Size of the chunks to read

    Notes
    -----
    The velocity and forces entries in the TRR file will not be read

    Returns
    -------
    xyz : np.ndarray, dtype=float32, shape=(n_frames, n_atoms, 3)
        The xyz coordinates
    time : np.ndarray, dtype=float32, shape=(n_frames)
    step :  np.ndarray, dtype=int32, shape=(n_frames)
    box : np.ndarray, dtype=float32, shape=(n_frames, 3, 3)
    lambd :  np.ndarray, dtype=float32, shape=(n_frames)
    """
    zipper = tuple(zip(*TRRReader(filename, chunk)))
    xyz = np.vstack(zipper[0])
    time = np.concatenate(zipper[1])
    step = np.concatenate(zipper[2])
    box = np.vstack(zipper[3])
    lambd = np.concatenate(zipper[4])

    return TRRFile(xyz, time, step, box, lambd)


def write(filename, xyz, time=None, step=None, box=None, lambd=None,
    force_overwrite=True):
    """
    Write a Gromacs TRR file

    Notes
    -----
    The velocities and the forces entries in the trr will be set to zeros

    Parameters
    ----------
    filename : str
    xyz : np.ndarray, dtype=float32, shape=(n_frames, n_atoms, 3)
        The xyz coordinates
    time : np.ndarray, dtype=float32, shape=(n_frames)
    step :  np.ndarray, dtype=int32, shape=(n_frames)
    box : np.ndarray, dtype=float32, shape=(n_frames, 3, 3)
    lambd :  np.ndarray, dtype=float32, shape=(n_frames)
    force_overwrite : bool
    """
    if force_overwrite and os.path.exists(filename):
        os.unlink(filename)

    # only overwrite if you really want to
    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    # make sure all the arrays are the right shape
    xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False)
    n_frames = len(xyz)

    step = ensure_type(step, dtype=np.int32, ndim=1, name='step', can_be_none=True,
        length=n_frames)
    if step is None:
        step = np.ones(n_frames, dtype=np.int32)

    time = ensure_type(time, dtype=np.float32, ndim=1, name='time', can_be_none=True,
        length=n_frames)
    if time is None:
        time = np.arange(n_frames, dtype=np.float32)

    box = ensure_type(box, dtype=np.float32, ndim=3, name='box', can_be_none=True,
        length=n_frames, shape=(n_frames, 3, 3))
    if box is None:
        # make each box[i] be the identity matrix
        box = np.zeros((n_frames, 3, 3), dtype=np.float32)
        box[:,0,0] = np.ones(n_frames, dtype=np.float32)
        box[:,1,1] = np.ones(n_frames, dtype=np.float32)
        box[:,2,2] = np.ones(n_frames, dtype=np.float32)

    lambd = ensure_type(lambd, dtype=np.float32, ndim=1, name='lambd', can_be_none=True,
        length=n_frames)
    if lambd is None:
        lambd = np.ones(n_frames, dtype=np.float32)

    writer = TRRWriter(filename)
    writer.write(xyz, time, step, box, lambd)


# code that indicates a sucessful return from the library
cdef int _EXDROK = 0             # OK
cdef int _EXDRHEADER = 1         # Header
cdef int _EXDRSTRING = 2         # String
cdef int _EXDRDOUBLE = 3         # Double
cdef int _EXDRINT = 4            # Integer
cdef int _EXDRFLOAT = 5          # Float
cdef int _EXDRUINT = 6           # Unsigned integer
cdef int _EXDR3DX = 7            # Compressed 3d coordinate
cdef int _EXDRCLOSE = 8          # Closing file
cdef int _EXDRMAGIC = 9          # Magic number
cdef int _EXDRNOMEM = 10         # Not enough memory
cdef int _EXDRENDOFFILE = 11     # End of file
cdef int _EXDRFILENOTFOUND = 12  # File not found

cdef class TRRReader:
    cdef XDRFILE *_xd
    cdef public int n_atoms
    cdef int chunk

    def __cinit__(self, filename, int chunk=1):
        # set self.n_atoms
        self.n_atoms = 0
        read_trr_natoms(filename, &self.n_atoms)

        # open file descriptor
        self._xd = xdrfile_open(filename, 'r')
        if self._xd is NULL:
            raise IOError("File not found: %s" % filename)

        self.chunk = chunk

    def __dealloc__(self):
        if self._xd is not NULL:
            xdrfile_close(self._xd)

    def __iter__(self):
        return self

    @cython.boundscheck(False)
    def __next__(self):
        cdef int status = _EXDROK
        cdef int i = 0

        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz  = np.empty((self.chunk, self.n_atoms, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time = np.empty(self.chunk, dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.int32_t, mode='c']   step = np.empty(self.chunk, dtype=np.int32)
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box  = np.empty((self.chunk, 3, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] lambd = np.empty(self.chunk, dtype=np.float32)

        while (i < self.chunk) and (status != _EXDRENDOFFILE):
            # the cython compiler seems to want use to explicitly cast from
            # int32_t* to int*. I know it's a little ugly, but I did assert
            # at the top of this file that they're the same size
            status = read_trr(self._xd, self.n_atoms, <int*> &step[i], &time[i],
                &lambd[i], &box[i,0,0], &xyz[i, 0, 0], NULL, NULL)

            if (status != _EXDRENDOFFILE) and status != _EXDROK:
                raise RuntimeError("TRR Read error: %s." % status)

            i += 1

        if status == _EXDRENDOFFILE:
            # if the file is over and we didn't read any data, raise
            # the stop itetation
            if i == 1:
                raise StopIteration
            # otherwise, return the data we have. since no data was read in
            # the last iteration, we need to chop that off
            xyz = xyz[0:i-1]
            box = box[0:i-1]
            time = time[0:i-1]
            box = box[0:i-1]
            lambd = lambd[0:i-1]

        return xyz, time, step, box, lambd


cdef class TRRWriter:
    cdef XDRFILE* fh

    def __cinit__(self, char* filename):
        self.fh = xdrfile_open(filename, 'w')
        if self.fh == NULL:
            raise IOError("Unable to open file %s" % filename)


    def __dealloc__(self):
        xdrfile_close(self.fh)


    @cython.boundscheck(False)
    def write(self,     np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz not None,
                        np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time not None,
                        np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step not None,
                        np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box not None,
                        np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] lambd not None):
        cdef int n_frames = len(xyz)
        cdef int n_atoms = xyz.shape[1]
        cdef int status, i

        # all same shape
        assert n_frames == len(box) == len(step) == len(time) == len(lambd)

        cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] force = np.zeros((n_atoms, 3), dtype=np.float32)
        cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] vel = np.zeros((n_atoms, 3), dtype=np.float32)

        for i in range(n_frames):
            status = write_trr(self.fh, n_atoms, step[i], time[i], lambd[i],
                &box[i, 0, 0], &xyz[i, 0, 0], &force[0,0], &vel[0,0])
            if status != _EXDROK:
                raise RuntimeError('TRR write error: %s' % status)

        return status
