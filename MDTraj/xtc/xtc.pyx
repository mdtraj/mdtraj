import cython
import os
from itertools import izip
import numpy as np
cimport numpy as np
import numpy as np
np.import_array()
cimport xdrlib

def read(filename, chunk=1):
    """Read the xyz coordinates from a Gromacs XTC file

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
    box = np.vstack(zipper[1])
    time = np.concatenate(zipper[2])
    prec = np.concatenate(zipper[3])
    step = np.concatenate(zipper[4])
    
    return xyz, time, step, box, prec


def write(filename, xyz, time=None, step=None, box=None, prec=None, force_overwrite=False):
    """Write a Gromacs XTC file
    
    """
    if force_overwrite and os.path.exists(filename):
        os.unlink(filename)
        
    # only overwrite if you really want to
    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    def ensure_type(val, dtype, ndim, length=None, can_be_none=False, shape=None):
        "Ensure dtype and shape of an ndarray"
        if can_be_none and val is None:
            return None
        if not isinstance(val, np.ndarray):
            raise TypeError("Must be numpy array")
        val = np.ascontiguousarray(val, dtype=dtype)
        if not val.ndim == ndim:
            raise ValueError('ndim is wrong')
        if length is not None and len(val) != length:
            raise ValueError('Length is not right. Got %s, should be %s' % (len(val), length))
        if shape is not None and val.shape != shape:
            raise ValueError('Wrong shape. Got %s, should be %s' % (val.shape, shape))
            
        return val
    
    # make sure all the arrays are the right shape
    xyz = ensure_type(xyz, dtype=np.float32, ndim=3, can_be_none=False)
    n_frames = len(xyz)
    
    step = ensure_type(step, dtype=np.int32, ndim=1, can_be_none=True,
        length=n_frames)
    if step is None:
        step = np.arange(n_frames, dtype=np.int32)
        
    time = ensure_type(time, dtype=np.float32, ndim=1, can_be_none=True,
        length=n_frames)
    if time is None:
        time = np.arange(n_frames, dtype=np.float32)
        
    box = ensure_type(box, dtype=np.float32, ndim=3, can_be_none=True,
        length=n_frames, shape=(n_frames, 3, 3))
    if box is None:
        box = np.zeros((n_frames, 3, 3), dtype=np.float32)
        
    prec = ensure_type(prec, dtype=np.float32, ndim=1, can_be_none=True,
        length=n_frames)
    if prec is None:
        prec = np.zeros(n_frames, dtype=np.float32)
                    
                
    writer = XTCWriter(filename, xyz, step, time, box, prec)
    writer.write()


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
        xdrlib.read_xtc_natoms(filename, cython.address(self.n_atoms))

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
        cdef int position = 0

        cdef np.ndarray xyz = np.zeros((self.chunk, self.n_atoms, 3), dtype=np.float32)
        cdef np.ndarray box = np.zeros((self.chunk, 3, 3), dtype=np.float32)
        cdef np.ndarray time = np.zeros(self.chunk, dtype=np.float32)
        cdef np.ndarray prec = np.zeros(self.chunk, dtype=np.float32)                      
        cdef np.ndarray step = np.zeros(self.chunk, dtype=np.int32)

        cdef float* xyz_p = <float*> xyz.data
        cdef float* box_p = <float*> box.data
        cdef float* time_p = <float*> time.data
        cdef float* prec_p = <float*> prec.data
        cdef int* step_p = <int*> step.data

        while (position < self.chunk) and (status != _EXDRENDOFFILE):
            status = xdrlib.read_xtc(self._xd, self.n_atoms, step_p,
                time_p, box_p, xyz_p, prec_p)

            if status == _EXDR3DX:
                raise RuntimeError("Decompression error in xrd lib")
            #print 'pos', position, 'status', status
            
            position += 1                
            xyz_p += self.n_atoms * 3
            box_p += 9
            time_p += 1
            step_p += 1
            prec_p += 1
            
        
        if status == _EXDRENDOFFILE:
            # if the file is over and we didn't read any data, raise
            # the stop itetation
            if position == 1:
                raise StopIteration
            # otherwise, return the data we have. since no data was read in
            # the last iteration, we need to chop that off
            xyz = xyz[0:position-1]
            box = box[0:position-1]
            time = time[0:position-1]
            prec = prec[0:position-1]
            step = step[0:position-1]
        
        
        return xyz, box, time, prec, step


cdef class XTCWriter:
    cdef xdrlib.XDRFILE* fh
    cdef np.ndarray xyz
    cdef np.ndarray step
    cdef np.ndarray time
    cdef np.ndarray box
    cdef np.ndarray prec
    
    def __cinit__(self, char* filename, np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz,
                    np.ndarray[ndim=1, dtype=np.int32_t] step,
                    np.ndarray[ndim=1, dtype=np.float32_t] time,
                    np.ndarray[ndim=3, dtype=np.float32_t] box,
                    np.ndarray[ndim=1, dtype=np.float32_t] prec):
        self.fh = xdrlib.xdrfile_open(filename, 'w')

        self.xyz = xyz
        self.step = step
        self.time = time
        self.box = box
        self.prec = prec


    def __dealloc(self):
        xdrlib.xdrfile_close(self.fh)


    def write(self):
        cdef int n_frames = len(self.xyz)
        cdef int status
        cdef int n_atoms = self.xyz.shape[1]
        cdef float* xyz_p = <float*> self.xyz.data
        cdef float* box_p = <float*> self.box.data

        cdef int i
        for i in range(n_frames):
            status = xdrlib.write_xtc(self.fh, n_atoms, self.step[i],
                self.time[i], box_p, xyz_p, self.prec[i])
            if status != _EXDROK:
                raise RuntimeError('XTC error: %s' % status)

            box_p += 9
            xyz_p += n_atoms * 3
        