# cython: c_string_type=str, c_string_encoding=ascii

from libc.stdlib cimport malloc, free
from libc.string cimport strcpy, strlen
from mdtraj.formats.registry import FormatRegistry
import numpy as np
cimport numpy as np

# TODO: Refactor molfile
cdef extern from "../dcd/include/molfile_plugin.h":
    ctypedef struct molfile_timestep_t:
        float *coords  # coordinates of all atoms, arranged xyzxyzxyz
        float *velocities  # space for velocities of all atoms; same layout
        float A, B, C, alpha, beta, gamma

    ctypedef struct molfile_atom_t:
        char name[16]
        char type[16]
        char resname[8]
        int resid
        char segid[8]
        char chain[2]

        char altloc[2]
        char insertion[2]
        float occupancy
        float bfactor
        float mass
        float charge
        float radius
        int atomicnumber

cdef int MOLFILE_NOOPTIONS=0x0000
cdef int MOLFILE_INSERTION=     0x0001
cdef int MOLFILE_OCCUPANCY=   0x0002
cdef int MOLFILE_BFACTOR=     0x0004
cdef int MOLFILE_MASS=        0x0008
cdef int MOLFILE_CHARGE=      0x0010
cdef int MOLFILE_RADIUS=      0x0020
cdef int MOLFILE_ALTLOC=      0x0040
cdef int MOLFILE_ATOMICNUMBER=0x0080
cdef int MOLFILE_BONDSSPECIAL=0x0100

cdef extern from 'include/maeffplugin.h':
    void * open_file_read(const char *fname, const char *ftype, int *vmdatoms)
    int read_structure(void *, int *optflags, molfile_atom_t *atoms)
    void close_file_read(void *v)
    int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts)

cdef class MAETrajectoryFile:

    cdef int n_atoms
    cdef void* fh
    cdef char* mode
    cdef char* filename
    cdef int is_open, _needs_write_initialization
    cdef readonly char* distance_unit

    cdef molfile_atom_t * structure

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True):
        """Open a DCD Trajectory File
        """
        self.distance_unit = 'angstroms'
        self.is_open = False
        self.mode = mode

        # Note: We copy this string because we do so for good reason in the
        # dcd plugin. Might be a problem here too. See the note there for details
        self.filename = <char*>malloc(strlen(filename)+1)
        if self.filename is NULL:
            raise MemoryError()
        strcpy(self.filename, filename)


        if str(mode) == 'r':
            self.fh = open_file_read(filename, "dcd", &self.n_atoms)
            if self.fh is NULL:
                raise IOError("Could not open file: %s" % filename)
            assert self.n_atoms > 0, 'MAE Corruption: n_atoms was not positive'
            # we're at the beginning of the file now
            self.is_open = True
        elif str(mode) == 'w':
            raise NotImplementedError("Can't write mae yet")
        else:
            raise ValueError("most must be one of ['r', 'w']")

    def __dealloc__(self):
        # free whatever we malloced
        free(self.filename)
        self.close()

    def close(self):
        """Close the DCD file handle
        """
        if self.is_open and self.fh is not NULL:
            if str(self.mode) == 'r':
                close_file_read(self.fh)
                pass
            else:
                raise NotImplementedError()
            self.is_open = False

        self._needs_write_initialization = False

    def read_topology(self):
        from mdtraj.core.topology import Topology
        from mdtraj.core.element import Element
        cdef molfile_atom_t * atoms
        atoms = <molfile_atom_t*>malloc(sizeof(molfile_atom_t)*self.n_atoms)
        if atoms is NULL:
            raise MemoryError()

        cdef int optflags;
        cdef int ret = read_structure(self.fh, &optflags, atoms)
        if ret:
            raise RuntimeError(ret)
        assert optflags & MOLFILE_ATOMICNUMBER > 0, 'Need atomic number'

        chains = {}
        residues = {}
        top = Topology()

        cdef int i
        for i in range(self.n_atoms):
            chain = atoms[i].chain
            resid = atoms[i].resid

            if chain not in chains:
                chains[chain] = top.add_chain()

            if (chain, resid) not in residues:
                residues[chain, resid] = top.add_residue(atoms[i].resname, chains[chain], resid, atoms[i].segid)

            top.add_atom(atoms[i].name, Element.getByAtomicNumber(atoms[i].atomicnumber), residues[chain,resid])
        free(atoms)
        return top

    def read_xyz(self):
        cdef molfile_timestep_t ts
        cdef float[:,::1] coords = np.empty((self.n_atoms, 3), dtype=np.float32)
        cdef float[:,::1] vels = np.empty((self.n_atoms, 3), dtype=np.float32)
        ts.coords = &coords[0,0]
        ts.velocities = &vels[0,0]
        cdef int ret = read_next_timestep(self.fh, self.n_atoms, &ts)
        if ret:
            raise RuntimeError(ret)
        return coords, vels, ts.A, ts.B, ts.C, ts.alpha, ts.beta, ts.gamma

    def read(self):
        from mdtraj.core.trajectory import Trajectory
        top = self.read_topology()
        coords, vels, A, B, C, alpha, beta, gamma = self.read_xyz()
        ang_to_nm = 10
        return Trajectory(np.array(coords)[np.newaxis, :, :] / ang_to_nm,
                          top,
                          time=np.array([0]),
                          unitcell_lengths=np.array([A, B, C]),
                          unitcell_angles=np.array([alpha, beta, gamma]) * np.pi / 180)

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()


FormatRegistry.register_fileobject('.mae')(MAETrajectoryFile)

@FormatRegistry.register_loader('.mae')
def load_mae(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_mae(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a MAE file from disk.

    Parameters
    ----------
    filename : str
        String filename of MAE file.
    top : unused
        MAE containes topology information. Unused. Must be None
    stride : unused
        MAE can only hold one frame. Must be None or 1
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it
        requires an extra copy, but will save memory.
    frame : unused
        MAE can only hold one frame. Must be None or 0

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.MAETrajectoryFile :  Low level interface to MAE files
    """
    if top is not None:
        raise ValueError("MAE contains topology information")
    if not (stride is None or stride == 1):
        raise ValueError("Don't stride! MAE only contains one frame")
    if not (frame is None or frame == 0):
        raise ValueError("MAE only contains one frame")

    with MAETrajectoryFile(filename, 'r') as f:
        traj =  f.read()
        if atom_indices is not None:
            traj.atom_slice(atom_indices, inplace=True)
        return traj
