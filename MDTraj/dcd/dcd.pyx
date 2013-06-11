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
from dcdlib cimport molfile_timestep_t, dcdhandle
from dcdlib cimport open_dcd_read, close_file_read, read_next_timestep
from dcdlib cimport open_dcd_write, close_file_write, write_timestep
from collections import namedtuple
DCDFile = namedtuple('DCDFile', ['xyz', 'box_lengths', 'box_angles'])

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

cdef ERROR_MESSAGES = {
    -1: 'Normal EOF',
    -2: 'DCD file does not exist',
    -3: 'Open of DCD file failed',
    -4: 'Read call on DCD file failed',
    -5: 'Premature EOF found in DCD file',
    -6: 'Format of DCD file is wrong',
    -7: 'Output file already exists',
    -8: 'Malloc failed',
    -9: 'Write call on DCD file failed',
}


def read(filename):
    """Read the data from a NAMD/CHARMM DCD file

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
    xyz, box_lengths, box_angles = DCDReader(filename).read()

    return DCDFile(xyz, box_lengths, box_angles)


def write(filename, xyz, box_lengths=None, box_angles=None, force_overwrite=True):
    """Write data to a NAMD/CHARMM DCD file

    Note that the box size entries in the DCD file will be left blank (zeros)

    Parameters
    ----------
    filename : str
        The filename of the dcd file to write to
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
        The xyz coordinates of each atom in each frame
    box_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
        The length of the box in each frame. `box_lengths[i,0]` is the length
        of the A axis (in frame i), and `box_lengths[i,1]` and
        `box_lengths[i,2]` are the B and C axis respectively. If unsupplied,
        all of the box lengths will be set to 1
    box_angles : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
        Organized analogously to box_lengths. Gives the alpha, beta and
        gamma angles respectively in entries `box_angles[i,0]`,
        `box_angles[i,1]`, `box_angles[i,2]`. If unsupplied, all of the box
        angles will be set to 90.0.
    force_overwrite : bool, default=False
        Overwrite anything that exists at filename, if its already there
    """

    if not force_overwrite and os.path.exists(filename):
        raise IOError('The file already exists: %s' % filename)

    # make sure all the arrays are the right shape
    xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False)
    n_frames = len(xyz)

    box_lengths = ensure_type(box_lengths, dtype=np.float32, ndim=2, name='box_lengths',
        can_be_none=True, shape=(n_frames, 3))
    if box_lengths is None:
        box_lengths = np.ones((n_frames, 3), dtype=np.float32)

    box_angles = ensure_type(box_angles, dtype=np.float32, ndim=2, name='box_angles',
        can_be_none=True, shape=(n_frames, 3))
    if box_angles is None:
        box_angles = 90.0 * np.ones((n_frames, 3), dtype=np.float32)

    writer = DCDWriter(filename, xyz, box_lengths, box_angles)
    writer.write()


cdef class DCDReader:
    # these variables hold the number of atoms and number of frames in the
    # file, as read off the header of the DCD file
    cdef int n_atoms, n_frames
    # the number of frames we're currently read off the file.
    cdef int frame_counter
    # this is the filehandle
    cdef dcdhandle* fh
    # C struct into which the library deposits all the data from the
    # most recent timestep
    cdef molfile_timestep_t* timestep

    def __cinit__(self, char* filename):
        # open the file handle and read n_atoms and n_frames
        self.fh = open_dcd_read(filename, "dcd", &self.n_atoms, &self.n_frames)
        if self.fh is NULL:
            raise IOError('There was an error opening the dcd file: %s' % filename)
        assert self.n_atoms > 0, 'DCD Corruption: n_atoms was not positive'
        assert self.n_frames >= 0, 'DCD corruption: n_frames < 0'
        # alloc the molfile_timestep, which is the struct that the library is
        # going to deposit its data into each timestep
        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.fh is  NULL:
            raise MemoryError('Malloc failed.')

        # we've currently read zero frames...
        self.frame_counter = 0


    def __dealloc__(self):
        "Shut this whole thing down"

        # free whatever we malloced
        free(self.timestep)

        # close the file
        if self.fh is not NULL:
            close_file_read(self.fh)

    @cython.boundscheck(False)
    def read(self, n_frames=None):
        """Extract the coordinates from a DCD file

        Parameters
        ----------
        n_frames : int, optional
            If not None, then read only the next `n_frames` frames.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
            The xyz coordinates of each atom in each frame
        box_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32
            The length of the box in each frame. `box_lengths[i,0]` is the length
            of the A axis (in frame i), and `box_lengths[i,1]` and
            `box_lengths[i,2]` are the B and C axis respectively.
        box_angles : np.ndarray, shape=(n_frames, 3), dtype=float32
            Organized analogously to box_lengths. Gives the alpha, beta and
            gamma angles respectively in entries `box_angles[i,0]`,
            `box_angles[i,1]`, `box_angles[i,2]`.
        """
        
        if n_frames is None:
            # if the user specifies n_frames=None, they want to read to the
            # end of the file
            n_frames = self.n_frames - self.frame_counter

        # malloc space to put the data that we're going to read off the disk
        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((n_frames, self.n_atoms, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] box_lengths = np.zeros((n_frames, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] box_angles = np.zeros((n_frames, 3), dtype=np.float32)

        cdef int i = 0
        cdef int status = _DCD_SUCCESS

        for i in range(n_frames):
            self.timestep.coords = &xyz[i,0,0]
            status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            self.frame_counter += 1
            box_lengths[i, 0] = self.timestep.A
            box_lengths[i, 1] = self.timestep.B
            box_lengths[i, 2] = self.timestep.C
            box_angles[i, 0] = self.timestep.alpha
            box_angles[i, 1] = self.timestep.beta
            box_angles[i, 2] = self.timestep.gamma

            if status != _DCD_SUCCESS:
                # if the frame was not successfully read, then we're done
                break

        if status == _DCD_SUCCESS:
            # if we're done either because of we read all of the n_frames
            # requested succcessfully, return
            return xyz, box_lengths, box_angles

        if status == _DCD_EOF:
            # if we're doing because we reached a normal EOF (perhaps the)
            # user asked to read more frames than were in the file, we need
            # to truncate the return arrays -- we don't want to return them
            # a big stack of zeros.
            xyz = xyz[0:i]
            box_lengths = box_lengths[0:i]
            box_angles = box_angles[0:i]
            return xyz, box_lengths, box_angles
    
        # If we got some other status, thats a "real" error.
        raise IOError("Error: %s", ERROR_MESSAGES(status))

        


cdef class DCDWriter:
    cdef char* filename
    cdef dcdhandle* fh
    cdef np.ndarray xyz
    cdef np.ndarray box_lengths
    cdef np.ndarray box_angles
    cdef int n_atoms
    cdef molfile_timestep_t* timestep

    def __cinit__(self, char* filename,
        np.ndarray[np.float32_t, ndim=3, mode="c"] xyz,
        np.ndarray[np.float32_t, ndim=2, mode="c"] box_lengths,
        np.ndarray[np.float32_t, ndim=2, mode="c"] box_angles):
        """
        Set up the DCD writer

        Nothing will be written to disk until you call write()

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
            The xyz coordinates of each atom in each frame
        box_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32
            The length of the box in each frame. `box_lengths[i,0]` is the length
            of the A axis (in frame i), and `box_lengths[i,1]` and
            `box_lengths[i,2]` are the B and C axis respectively.
        box_angles : np.ndarray, shape=(n_frames, 3), dtype=float32
            Organized analogously to box_lengths. Gives the alpha, beta and
            gamma angles respectively in entries `box_angles[i,0]`,
            `box_angles[i,1]`, `box_angles[i,2]`.
        """
        self.filename = filename
        self.xyz = xyz
        self.box_lengths = box_lengths
        self.box_angles = box_angles
        self.n_atoms = self.xyz.shape[1]

        assert self.box_lengths.shape[0] == len(self.xyz)
        assert self.box_lengths.shape[1] == 3
        assert self.box_angles.shape[0] == len(self.xyz)
        assert self.box_angles.shape[1] == 3

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

    @cython.boundscheck(False)
    def write(self):
        "Write all the data"
        cdef int i, status
        cdef n_frames = len(self.xyz)

        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = self.xyz

        for i in range(n_frames):
            self.timestep.coords = &xyz[i, 0, 0]
            self.timestep.A = self.box_lengths[i, 0]
            self.timestep.B = self.box_lengths[i, 1]
            self.timestep.C = self.box_lengths[i, 2]
            self.timestep.alpha = self.box_angles[i, 0]
            self.timestep.beta  = self.box_angles[i, 1]
            self.timestep.gamma = self.box_angles[i, 2]

            status = write_timestep(self.fh, self.timestep)

            if status != _DCD_SUCCESS:
                raise IOError("DCD Error: %s" % ERROR_MESSAGES(status))
