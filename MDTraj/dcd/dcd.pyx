import cython
import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from dcdlib cimport molfile_timestep_t, dcdhandle
from dcdlib cimport open_dcd_read, close_file_read, read_next_timestep

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
    """Read the xyz coordinates from a NAMD/CHARMM DCD File

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
    """
    return DCDReader(filename).read()


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
