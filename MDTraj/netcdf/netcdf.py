import os
import numpy as np
from mdtraj import version
from mdtraj.utils.arrays import ensure_type

try:
    import netCDF4 as netcdf
except ImportError:
    print '#'*73
    print 'ERROR'
    print '#'*73
    print 'MDTraj\'s AMBER NetCDF bindings require the python-netcdf'
    print 'library. You can install the library using the python package'
    print 'managers "pip" or "easy_install" with'
    print ''
    print 'pip install netcdf4'
    print 'or'
    print 'easy_install netcdf4'
    print ''
    print 'You can also download the library directly and find more documentation at'
    print 'https://pypi.python.org/pypi/netCDF4'
    print '#'*73
    raise


class NetCDFFile(object):
    def __init__(self, filename, mode='r', force_overwrite=False):
        if mode not in ['r', 'w', 'a', 'ws', 'as']:
            raise ValueError(("mode must be one of ['r', 'w', 'a', 'ws', 'as']"
                " 'r' indicates read, 'w' indicates write, and 'a' indicates"
                " append. 'a' and 'w' can be appended with 's', which turns "
                " off buffering"))

        # AMBER uses the NetCDF3 format, with 64 bit encodings
        self._handle = netcdf.Dataset(filename, mode=mode, format='NETCDF3_64BIT',
                                     clobber=force_overwrite)
        
        # the current frame that we're at in the file
        self._frame_index = 0
        
        # _needs_initialization indicates whether we need to set the
        # global properties of the file. This is required before the first
        # write operation on a new file
        self._needs_initialization = False

        if mode in ['w', 'ws']:
            self._needs_initialization = True

    def __del__(self):
        if hasattr(self, '_handle'):
            self._handle.close()
        
    def flush(self):
        "Write all buffered data in the to the disk file."
        self._handle.sync()

    def _initialize_headers(self, n_atoms):
        """Initialize the NetCDF file according to the AMBER NetCDF Convention,
        Version 1.0, revision B.

        The convention is defined here: http://ambermd.org/netcdf/nctraj.html
        """

        # Set attributes.
        setattr(self._handle, 'title', '')
        setattr(self._handle, 'application', 'AMBER')
        setattr(self._handle, 'program', 'MDTraj')
        setattr(self._handle, 'programVersion', version.short_version)
        setattr(self._handle, 'Conventions', 'AMBER')
        setattr(self._handle, 'ConventionVersion', '1.0')


        # set the dimensions
        # unlimited number of frames in trajectory
        self._handle.createDimension('frame', 0)
        # number of spatial coordinates
        self._handle.createDimension('spatial', 3)
        # number of atoms
        self._handle.createDimension('atom', n_atoms)
        # three spatial coordinates for the length of the unit cell
        self._handle.createDimension('cell_spatial', 3)
        # three spatial coordinates for the angles that define the shape
        # of the unit cell
        self._handle.createDimension('cell_angular', 3)
        # length of the longest string used for a label
        self._handle.createDimension('label', 5)
        
        # Define variables to store unit cell data
        cell_spatial = self._handle.createVariable('cell_spatial', 'c', ('cell_spatial',))
        cell_angular = self._handle.createVariable('cell_angular', 'c', ('cell_spatial', 'label'))
        cell_lengths = self._handle.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
        setattr(cell_lengths, 'units', 'angstrom')
        cell_angles = self._handle.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
        setattr(cell_angles, 'units', 'degree')

        self._handle.variables['cell_spatial'][0] = 'x'
        self._handle.variables['cell_spatial'][1] = 'y'
        self._handle.variables['cell_spatial'][2] = 'z'

        self._handle.variables['cell_angular'][0] = 'alpha'
        self._handle.variables['cell_angular'][1] = 'beta '
        self._handle.variables['cell_angular'][2] = 'gamma'

        # Define coordinates and snapshot times.
        frame_times = self._handle.createVariable('time', 'f', ('frame',))
        setattr(frame_times, 'units', 'picosecond')
        frame_coordinates = self._handle.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
        setattr(frame_coordinates, 'units', 'angstrom')


    def read(self, n_frames=None):
        pass


    def write(self, coordinates, time=None, cell_lengths=None, cell_angles=None):

        # typecheck all of the input arguments rigorously
        coordinates = ensure_type(coordinates, np.float32, 3, 'coordinates', length=None,
            can_be_none=False, shape=(None, None, 3), warn_on_cast=True, add_newaxis_on_deficient_ndim=True)
        n_frames, n_atoms = coordinates.shape[0], coordinates.shape[1]

        time = ensure_type(time, np.float32, 1, 'time', length=n_frames,
            can_be_none=True, warn_on_cast=True, add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, np.float32, 2, 'cell_lengths', length=n_frames,
            can_be_none=True, shape=(n_frames, 3), warn_on_cast=True, add_newaxis_on_deficient_ndim=True)
        cell_angles = ensure_type(cell_angles, np.float32, 2, 'cell_angles', length=n_frames,
            can_be_none=True, shape=(n_frames, 3), warn_on_cast=True, add_newaxis_on_deficient_ndim=True)

        # are we dealing with a periodic system?
        if cell_lengths is None and cell_angles is None:
            is_periodic = False
        elif cell_lengths is not None and cell_angles is not None:
            is_periodic = True
        else:
            provided, neglected = 'cell_lengths', 'cell_angles'
            if cell_lengths is None:
                provided, neglected = neglected, provided
            raise ValueError('You provided the variable "%s", but neglected to '
                             'provide "%s". They either BOTH must be provided, or '
                             'neither. Having one without the other is meaningless' % (
                                provided, neglected))

        if self._needs_initialization:
            self._initialize_headers(n_atoms)

        # this slice object says where we're going to put the data in the
        # arrays
        frameslice = slice(self._frame_index, self._frame_index + n_frames)
        
        # update the frame index pointers
        self._frame_index += n_frames
        
        # deposit the data
        self._handle.variables['coordinates'][frameslice, : :] = coordinates
        if time is not None:
            self._handle.variables['time'][frameslice] = time
        if cell_lengths is not None:
            self._handle.variables['cell_lengths'][frameslice, :] = cell_lengths
        if cell_angles is not None:
            self._handle.variables['cell_angles'][frameslice, :] = cell_angles

        return
        

if __name__ == '__main__':
    f = NetCDFFile('file.nc', 'w', force_overwrite=True)
    f.write(coordinates=np.random.randn(3,3), cell_angles=np.random.randn(3),
        cell_lengths=np.random.randn(3), time=1)
