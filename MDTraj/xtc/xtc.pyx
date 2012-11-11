import cython
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
    
    return xyz, box, time, prec, step


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
            
        self.chunk = chunk

    def __dealloc__(self):
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

        while (status == _EXDROK and position < self.chunk):
            status = xdrlib.read_xtc(self._xd, self.n_atoms, step_p,
                time_p, box_p, xyz_p, prec_p)

            position += 1                
            xyz_p += self.n_atoms * 3
            box_p += 9
            time_p += 1
            step_p += 1
            prec_p += 1
        
        if status != _EXDROK:
            if position == 1:
                raise StopIteration
            
            xyz = xyz[0:position-1]
            box = box[0:position-1]
            time = time[0:position-1]
            prec = prec[0:position-1]
            step = step[0:position-1]
        
        return xyz, box, time, prec, step