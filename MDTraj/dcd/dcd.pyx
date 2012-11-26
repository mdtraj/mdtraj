import cython
import os, warnings
import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from dcdlib cimport molfile_timestep_t, dcdhandle
from dcdlib cimport open_dcd_read, close_file_read, read_next_timestep
from dcdlib cimport open_dcd_write, close_file_write, write_timestep


# codes that indicate status on return from library
cdef int _DCD_SUCCESS    = 0   # No problems
cdef int _DCD_EOF        = -1  # Normal EOF
cdef int _DCD_DNE        = -2  # DCD file does not exist
cdef int _DCD_OPENFAILED = -3  # Open of DCD file failed
cdef int _DCD_BADREAD    = -4  # read call on DCD file failed
cdef int _DCD_BADEOF     = -5  # premature EOF found in DCD file
cdef int _DCD_BADFORMAT  = -6  # format of DCD file is wrong
cdef int _DCD_FILEEXISTS = -7  # output file already exists
cdef int _DCD_BADMALLOC  = -8  # malloc failed
cdef int _DCD_BADWRITE   = -9  # write call on DCD file failed

def read_xyz(filename):
    """Read the xyz coordinates from a NAMD/CHARMM DCD file

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
    return DCDReader(filename).read()

def write_xyz(filename, xyz, force_overwrite=False):
    """Write xyz coordinates to a NAMD/CHARMM DCD file
    
    Note that the box size entries in the DCD file will be left blank (zeros)
    
    Parameters
    ----------
    filename : str
        The filename of the dcd file to write to
    xyz : np.ndarray, ndim=3, dtype=np.float32
        The xyz coordinates
    force_overwrite : bool
        Overwrite anything that exists at filename, if its already there
    """
    
    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    if not isinstance(xyz, np.ndarray):
        raise TypeError("Must be numpy array")
    if xyz.dtype != np.float32:
        warnings.warn('Casting to float32')
        xyz = np.array(xyz, dtype=np.float32)
    if not xyz.flags.c_contiguous:
        warnings.warn('Casting to contiguous')
        xyz = np.ascontiguousarray(xyz, dtype=np.float32)

    writer = DCDWriter(filename, xyz)
    writer.write()

cdef class DCDReader:
    cdef int n_atoms, n_frames
    cdef dcdhandle* fh
    cdef molfile_timestep_t* timestep
    cdef int status
    
    def __cinit__(self, char* filename):
        #open the file
        self.fh = open_dcd_read(filename, "dcd", &self.n_atoms, &self.n_frames)
        if self.fh is NULL:
            raise IOError('There was an error opening the dcd file: %s' % filename)
        
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
        

    def read(self):
        """Extract the XYZ coordinates from the file
        
        Note: if we wanted to, we could get the box dimensions out as well. But
        they're kindof boring.
        """

        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((self.n_frames, self.n_atoms, 3), dtype=np.float32)
        self.timestep.coords = <float*> xyz.data
        cdef int position = 0
        cdef int status = _DCD_SUCCESS
        
        while (status == _DCD_SUCCESS and position < self.n_frames):
            status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            self.timestep.coords += 3*self.n_atoms
            position += 1
        
        if status != _DCD_SUCCESS:
            raise IOError("Error: %s", status)
        
        return xyz


cdef class DCDWriter:
    cdef char* filename
    cdef dcdhandle* fh
    cdef np.ndarray xyz
    cdef int n_atoms
    cdef molfile_timestep_t* timestep
    
    def __cinit__(self, char* filename, np.ndarray[np.float32_t, ndim=3, mode="c"] xyz):
        self.filename = filename
        self.xyz = xyz
        self.n_atoms = self.xyz.shape[1]
        
        self.fh = open_dcd_write(filename, "dcd", self.n_atoms)
        if self.fh is NULL:
            raise IOError('There was an error opening the file: %s' % filename)
            
        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.timestep is NULL:
            raise MemoryError
        
    def __dealloc__(self):
        close_file_write(self.fh)
        if self.timestep is not NULL:
            free(self.timestep)
    
    def write(self):
        cdef int i, status
        cdef n_frames = len(self.xyz)
        
        self.timestep.A = 0
        self.timestep.B = 0
        self.timestep.C = 0
        self.timestep.alpha = 0
        self.timestep.beta = 0 
        self.timestep.gamma = 0

        self.timestep.coords = <float*> self.xyz.data
        
        for i in range(n_frames):
            status = write_timestep(self.fh, self.timestep)
            self.timestep.coords += self.n_atoms * 3
            if status != _DCD_SUCCESS:
                raise RuntimeError("DCD Error: %s" % status)
        
        
        
        
