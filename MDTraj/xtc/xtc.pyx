#cython: boundscheck=False
#cython: wraparound=False

import os
from itertools import izip
import numpy as np
cimport numpy as np
import numpy as np
np.import_array()
cimport xdrlib
from cython.operator cimport dereference as deref

# numpy variable types include the specific numpy of bytes of each, but the c
# variables in our interface file don't. this could get bad if we're on a wierd
# machine, so lets make sure first
if sizeof(int) != sizeof(np.int32_t):
    raise RuntimeError('Integers on your compiler are not 32 bits. This is not good.')
if sizeof(float) != sizeof(np.float32_t):
    raise RuntimeError('Floats on your compiler are not 32 bits. This is not good')


def read(filename, chunk=1):
    """
    Read the xyz coordinates from a Gromacs XTC file

    Parameters
    ----------
    filename : str
        The filename of the xtc file to read from
    chunk : int
        Size of the chunks to read

    Returns
    -------
    xyz : np.ndarray, dtype=float32, shape=(n_frames, n_atoms, 3)
        The xyz coordinates
    box : np.ndarray, dtype=float32, shape=(n_frames, 3, 3)
    time : np.ndarray, dtype=float32, shape=(n_frames)
    prec :  np.ndarray, dtype=float32, shape=(n_frames)
    step :  np.ndarray, dtype=int32, shape=(n_frames)
    """
    zipper = tuple(izip(*XTCReader(filename, chunk)))
    xyz = np.vstack(zipper[0])
    time = np.concatenate(zipper[1])
    step = np.concatenate(zipper[2])
    box = np.vstack(zipper[3])
    prec = np.concatenate(zipper[4])

    return xyz, time, step, box, prec


def write(filename, xyz, time=None, step=None, box=None, prec=None, force_overwrite=False):
    """
    Write a Gromacs XTC file
    
    Parameters
    ----------
    filename : str
    xyz : np.ndarray
    time : np.ndarray, optional
    step : np.ndarray, optional
    box : np.ndarray, optional
    prec : np.ndarray, optional
    force_overwrite : bool
    """
    if force_overwrite and os.path.exists(filename):
        os.unlink(filename)

    # only overwrite if you really want to
    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    def ensure_type(val, dtype, ndim, name, length=None, can_be_none=False, shape=None):
        "Ensure dtype and shape of an ndarray"
        if can_be_none and val is None:
            return None
        if not isinstance(val, np.ndarray):
            raise TypeError("%s must be numpy array. You supplied type %s" % (name, type(val)))
        val = np.ascontiguousarray(val, dtype=dtype)
        if not val.ndim == ndim:
            raise ValueError('%s must be ndim %s. You supplied %s' % (name, ndim, val.ndim))
        if length is not None and len(val) != length:
            raise ValueError('%s must be length %s. You supplied %s' % (name, length, len(val)))
        if shape is not None and val.shape != shape:
            raise ValueError('%s must be shape %s. You supplied %s' % (name, shape, val.shape))

        return val

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

    prec = ensure_type(prec, dtype=np.float32, ndim=1, name='prec', can_be_none=True,
        length=n_frames)
    if prec is None:
        prec = 1000.0 * np.ones(n_frames, dtype=np.float32)


    writer = XTCWriter(filename)
    writer.write(xyz, time, step, box, prec)


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

cdef class XTCReader:
    cdef xdrlib.XDRFILE *_xd
    cdef public int n_atoms
    cdef int chunk

    def __cinit__(self, filename, int chunk=1):
        # set self.n_atoms
        self.n_atoms = 0
        xdrlib.read_xtc_natoms(filename, &self.n_atoms)

        # open file descriptor
        self._xd = xdrlib.xdrfile_open(filename, 'r')
        if self._xd is NULL:
            raise IOError("File not found: %s" % filename)

        self.chunk = chunk

    def __dealloc__(self):
        if self._xd is not NULL:
            xdrlib.xdrfile_close(self._xd)

    def __iter__(self):
        return self

    def __next__(self):
        cdef int status = _EXDROK
        cdef int i = 0

        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz  = np.empty((self.chunk, self.n_atoms, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time = np.empty(self.chunk, dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.int32_t, mode='c']   step = np.empty(self.chunk, dtype=np.int32)
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box  = np.empty((self.chunk, 3, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec = np.empty(self.chunk, dtype=np.float32)

        while (i < self.chunk) and (status != _EXDRENDOFFILE):
            # the cython compiler seems to want use to explicitly cast from
            # int32_t* to int*. I know it's a little ugly, but I did assert
            # at the top of this file that they're the same size
            status = xdrlib.read_xtc(self._xd, self.n_atoms, <int*> &step[i], &time[i],
                &box[i,0,0], &xyz[i, 0, 0], &prec[i])

            if status != _EXDRENDOFFILE and status != _EXDROK:
                raise RuntimeError("XTC Read error: %s." % status)
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
            prec = prec[0:i-1]
            step = step[0:i-1]

        return xyz, time, step, box, prec


cdef class XTCWriter:
    cdef xdrlib.XDRFILE* fh


    def __cinit__(self, char* filename):
        self.fh = xdrlib.xdrfile_open(filename, 'w')
        if self.fh == NULL:
            raise IOError("Unable to open file %s" % filename)


    def __dealloc__(self):
        xdrlib.xdrfile_close(self.fh)


    def write(self,     np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz not None,
                        np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time not None,
                        np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step not None,
                        np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box not None,
                        np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec not None):
        cdef int n_frames = len(xyz)
        cdef int n_atoms = xyz.shape[1]
        cdef int status, i

        # all same shape
        assert n_frames == len(box) == len(step) == len(time) == len(prec)


        for i in range(n_frames):
            status = xdrlib.write_xtc(self.fh, n_atoms, step[i], time[i], &box[i, 0, 0], &xyz[i, 0, 0], prec[i])
            if status != _EXDROK:
                raise RuntimeError('XTC write error: %s' % status)

        return status
