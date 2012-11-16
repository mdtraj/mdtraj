import cython
import os, warnings
import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from binposlib cimport molfile_timestep_t
from binposlib cimport open_binpos_read, close_file_read, read_next_timestep
from binposlib cimport open_binpos_write, close_file_write, write_timestep


# codes that indicate status on return from library
cdef int _BINPOS_SUCESS = 0  # regular exit code
cdef int _BINPOS_EOF = -1  # end of file (or error)

def read_xyz(filename, chunk=1):
    """Read the xyz coordinates from a NAMD/CHARMM DCD File

    Parameters
    ----------
    filename : str
        The filename of the dcd file to read from
    chunk : int
        Size of the chunks to read

    Returns
    -------
    xyz : np.ndarray, dtype=float32, shape=(n_frames, n_atoms, 3)
        The xyz coordinates
    """
    return np.concatenate(tuple(BINPOSReader(filename, chunk)))


#
# def write_xyz(filename, xyz, force_overwrite=False):
#     """Write xyz coordinates to a NAMD/CHARMM DCD File
#
#     Note that the box size entries in the DCD file will be left blank (zeros)
#
#     Parameters
#     ----------
#     filename : str
#         The filename of the dcd file to write to
#     xyz : np.ndarray, ndim=3, dtype=np.float32
#         The xyz coordinates
#     """
#
#     if not force_overwrite and os.path.exists(filename):
#         raise IOError('The file already exists: %s' % filename)
#
#     if not isinstance(xyz, np.ndarray):
#         raise TypeError("Must be numpy array")
#     if xyz.dtype != np.float32:
#         warnings.warn('Casting to float32')
#         xyz = np.array(xyz, dtype=np.float32)
#     if not xyz.flags.c_contiguous:
#         warnings.warn('Casting to contiguous')
#         xyz = np.ascontiguousarray(xyz, dtype=np.float32)
#
#     writer = BINPOSWriter(filename, xyz)
#     writer.write()

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


    def __next__(self):
        cdef int i
        cdef np.ndarray[dtype=np.float32_t, ndim=3, mode='c'] xyz = np.zeros((self.chunk, self.n_atoms, 3), dtype=np.float32)
        self.timestep.coords = <float*> xyz.data

        for i in range(self.chunk):
            self.timestep.coords = &xyz[i,0,0]
            status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            print i, status
            if status != _BINPOS_SUCESS:
                if i == 0:
                    raise StopIteration
                break

        if status != _BINPOS_SUCESS:
            return xyz[0:i]

        return xyz

#
# cdef class DCDWriter:
#     cdef char* filename
#     cdef dcdhandle* fh
#     cdef np.ndarray xyz
#     cdef int n_atoms
#     cdef molfile_timestep_t* timestep
#
#     def __cinit__(self, char* filename, np.ndarray[np.float32_t, ndim=3, mode="c"] xyz):
#         self.filename = filename
#         self.xyz = xyz
#         self.n_atoms = self.xyz.shape[1]
#
#         self.fh = open_dcd_write(filename, "dcd", self.n_atoms)
#         if self.fh is NULL:
#             raise IOError('There was an error opening the file: %s' % filename)
#
#         self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
#         if self.timestep is NULL:
#             raise MemoryError
#
#     def __dealloc__(self):
#         close_file_write(self.fh)
#         if self.timestep is not NULL:
#             free(self.timestep)
#
#     def write(self):
#         cdef int i, status
#         cdef n_frames = len(self.xyz)
#
#         self.timestep.A = 0
#         self.timestep.B = 0
#         self.timestep.C = 0
#         self.timestep.alpha = 0
#         self.timestep.beta = 0
#         self.timestep.gamma = 0
#
#         self.timestep.coords = <float*> self.xyz.data
#
#         for i in range(n_frames):
#             status = write_timestep(self.fh, self.timestep)
#             self.timestep.coords += self.n_atoms * 3
#             if status != _DCD_SUCCESS:
#                 raise RuntimeError("DCD Error: %s" % status)
#
#
#
#
