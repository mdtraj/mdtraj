import cython
import numpy as np
cimport numpy as np
cimport xdrlib

# code that indicates a sucessful return from the library
cdef int _EXDROK = 0

def load_xtc(filename, chunk_len=1000):
    
    cdef _XTCReader x = _XTCReader(filename, chunk_len)
    
    cdef int result = _EXDROK
    while result == _EXDROK: 
        result = x.get_next()
    
    # append the last chunk to the list
    for li, el in [(x.xyzl, x.xyz), (x.boxl, x.box), (x.timel, x.time),
                   (x.precl, x.prec), (x.stepl, x.step)]:
        # elements that never got filled in shouldn't be taken
        li.append(el[0:x.array_position-1])

    # concatenate the lists onto 5 numpy arrays to return
    return tuple(np.concatenate(li) for li in [x.xyzl, x.boxl, x.timel, x.precl, x.stepl])

    
cdef class _XTCReader:
    cdef xdrlib.XDRFILE *_xd  # the file handle
    cdef int n_atoms
    
    # there are five arrays in the XTC file itself. For each frame, the files lists
    # the coordinates, the time, the step, the box size and the precision
    
    # we don't know the number of frames in the data file a priori, so we use
    # a chunking strategy that combines the numpy arrays and python lists.
    # the idea is to keep a preallocated numpy array as a *buffer* that we add
    # each frames data to step by step. When the buffer fills up, we add the "full"
    # buffer to a python list and then allocate a new buffer. at the end, we concatenate
    # all of the numpy arrays in the list together

    # declare the buffers
    cdef np.ndarray xyz
    cdef np.ndarray time
    cdef np.ndarray step
    cdef np.ndarray box
    cdef np.ndarray prec
    
    # this is going to hold our current position in the buffer
    cdef int array_position
    
    # declare the lists
    cdef list xyzl
    cdef list boxl
    cdef list timel
    cdef list precl
    cdef list stepl
    
    # length of the buffer (major axis)
    cdef int chunk_len
    
    def __cinit__(self, filename, chunk_len):
        # set self.n_atoms
        xdrlib.read_xtc_natoms(filename, cython.address(self.n_atoms))
        # open file descriptor
        self._xd = xdrlib.xdrfile_open(filename, 'r')
        
        self.chunk_len = chunk_len
        
        # allocate working buffers
        self.xyz = np.zeros((self.chunk_len, self.n_atoms, 3), dtype=np.float32)
        self.box = np.zeros((self.chunk_len, 3, 3), dtype=np.float32)
        self.time = np.zeros(self.chunk_len, dtype=np.float32)
        self.prec = np.zeros(self.chunk_len, dtype=np.float32)
        self.step = np.zeros(self.chunk_len, dtype=np.int32)
        self.array_position = 0

        # allocate lists in which we put the filled buffers
        self.xyzl = []
        self.boxl = []
        self.timel = []
        self.precl = []
        self.stepl = []


    def __dealloc__(self):
        xdrlib.xdrfile_close(self._xd)


    cdef int get_next(self):
        # get pointers to the memory locations where the library will deposit the data
        cdef float* xyz_p = (<float*> self.xyz.data) + (self.array_position * self.n_atoms * 3)
        cdef float* box_p = (<float*> self.box.data) + (self.array_position * 9)
        cdef float* time_p = (<float*> self.time.data) + self.array_position
        cdef float* prec_p = (<float*> self.prec.data) + self.array_position
        cdef int* step_p = (<int*> self.step.data) + self.array_position

        cdef int result = xdrlib.read_xtc(self._xd, self.n_atoms, step_p, time_p, box_p, xyz_p, prec_p)

        self.array_position += 1
        
        # check if the buffers are full
        if (result == _EXDROK and self.array_position >= self.chunk_len):
            # push the now fill working buffers into their stacks
            self.xyzl.append(self.xyz)
            self.boxl.append(self.box)
            self.timel.append(self.time)
            self.precl.append(self.prec)
            self.stepl.append(self.step)

            # allocate new working buffers
            self.xyz = np.zeros((self.chunk_len, self.n_atoms, 3), dtype=np.float32)
            self.box = np.zeros((self.chunk_len, 3, 3), dtype=np.float32)
            self.time = np.zeros(self.chunk_len, dtype=np.float32)
            self.prec = np.zeros(self.chunk_len, dtype=np.float32)
            self.step = np.zeros(self.chunk_len, dtype=np.int32)
            self.array_position = 0

        return result