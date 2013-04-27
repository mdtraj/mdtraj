# Copyright 2013 mdtraj developers
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
"""
This module implements the MDTraj HDF5 format described in
https://github.com/rmcgibbo/mdtraj/issues/36
"""

##############################################################################
# Imports
##############################################################################

import sys
import operator
try:
    import simplejson as json
except ImportError:
    import json
import numpy as np


from mdtraj import version
import mdtraj.pdb.element as elem
from mdtraj.topology import Topology
from mdtraj.utils.unit import in_units_of
from mdtraj.utils.arrays import ensure_type

try:
    import simtk.unit as units
    HAVE_UNIT = True
except ImportError:
    HAVE_UNIT = False
    warnings.warn('The package simtk.unit was not imported, which means that'
                  'no unit processing will be done. Please install OpenMM to '
                  'get access to automatic unit conversion and validation. It '
                  'is highly recommended.')


def import_tables():
    """Delayed import of the tables (HDF5) module
    """

    try:
        if 'tables' not in sys.modules:
            import tables
        return sys.modules['tables']
    except ImportError:
        print '#'*73
        print 'ERROR'
        print '#'*73
        print 'MDTraj\'s HDF5 bindings require the tables library. You can'
        print 'install the library using the python package managers "pip" or '
        print '"easy_install" with'
        print ''
        print 'pip install netcdf4'
        print 'or'
        print 'easy_install netcdf4'
        print ''
        print 'You can also download the library directly and find more '
        print 'documentation at http://www.pytables.org'
        print '#'*73
        raise


##############################################################################
# Utilities
##############################################################################


def ensure_mode(*m):
    """This is a little decorator that is used inside HDF5Trajectory
    to validate that the file is open in the correct mode before doing
    a a method

    Parameters
    ----------
    m : str or list
        One or more of ['w', 'r', 'a'], giving the allowable modes
        for the method

    Example
    -------
    class HDF5Trajectory:
        @ensure_mode('w')
        def method_that_is_only_allowed_to_be_called_in_write_mode(self):
            print 'i must be in write mode!'
    """

    def inner(f):
        def wrapper(*args, **kwargs):
            # args[0] is self on the method
            if args[0]._mode in m:
                return f(*args)
            raise ValueError('This operation is only available when a file '
                             'is open in mode="%s".' % args[0]._mode)
        return wrapper
    return inner


##############################################################################
# Classes
##############################################################################


class HDF5Trajectory(object):
    def __init__(self, filename, mode='r', force_overwrite=False, compression='zlib'):
        self._open = False  # is the file handle currently open?
        self._mode = mode  # the mode in which the file was opened?

        if not mode in ['r', 'w', 'a']:
            raise ValueError("mode must be one of ['r', 'w', 'a']")

        tables = import_tables()

        if compression == 'zlib':
            compression = tables.Filters(complib='zlib', shuffle=True, complevel=1)
        elif compression is None:
            compression = None
        else:
            raise ValueError('compression must be either "zlib" or None')

        self._handle = tables.openFile(filename, mode=mode, filters=compression)
        self._open = True

        if mode == 'w':
            # what frame are we currently reading or writing at?
            self._frame_index = 0
            # do we need to write the header information?
            self._needs_initialization = True
        elif mode == 'a':
            self._frame_index = len(self._handle.root.coordinates)
            self._needs_initialization = False
        elif mode == 'r':
            self._frame_index = 0
            self._needs_initialization = False

    @ensure_mode('r')
    @property
    def root(self):
        """Direct access to the root group of the underlying Tables HDF5 file.

        This can be used for random or specific access to the underlying arrays
        on disk
        """
        return self._handle.root

    #####################################################
    # title global attribute (optional)
    #####################################################

    @property
    def title(self):
        if hasattr(self._handle.root._v_attrs, 'title'):
            return self._handle.root._v_attrs.title
        return None

    @title.setter
    @ensure_mode('w', 'a')
    def title(self, value):
        self._handle.root._v_attrs.title = str(value)

    #####################################################
    # application global attribute (optional)
    #####################################################

    @property
    def application(self):
        if hasattr(self._handle.root._v_attrs, 'application'):
            return self._handle.root._v_attrs.application
        return None

    @application.setter
    @ensure_mode('w', 'a')
    def application(self, value):
        self._handle.root._v_attrs.application = str(value)

    #####################################################
    # topology global attribute (optional, recommended)
    #####################################################

    @property
    def topology(self):
        """Get the topology out from the file

        Returns
        -------
        topology : mdtraj.Topology
            A topology object
        """
        if not hasattr(self._handle.root._v_attrs, 'topology'):
            return None
        topology_dict = json.loads(self._handle.root._v_attrs.topology)
        topology = Topology()

        for chain_dict in sorted(topology_dict['chains'], key=operator.attrgetter('index')):
            chain = topology.addChain()
            for residue_dict in sorted(chain_dict['residues'], key=operator.attrgetter('index')):
                residue = topology.addResidue(residue_dict['name'], chain)
                for atom_dict in sorted(residue_dict['residue_dict'], key=operator.attrgetter('index')):
                    try:
                        element =  elem.get_by_symbol(atom_dict['element'])
                    except KeyError:
                        raise ValueError('The symbol %s isn\'t a valid element' % atom_dict['element'])
                    topology.addAtom(atom_dict['name'], element, residue)

        atoms = list(topology.atoms())
        for index1, index2 in topology_dict['atoms']:
            topology.addBond(atoms[index1], atoms[index2])

        return topology

    @topology.setter
    @ensure_mode('w', 'a')
    def topology(self, topology_object):
        """Set the topology in the file

        Parameters
        -------
        topology_object : mdtraj.Topology
            A topology object
        """

        try:
            topology_dict = {
                'chains': [],
                'bonds': []
            }
            for chain in topology_object.chains():
                chain_dict = {
                    'residues': [],
                    'index': int(chain.index)
                }
                for residue in chain.residues():
                    residue_dict = {
                        'index': int(residue.index),
                        'name': str(residue.name),
                        'atoms': []
                    }
                    for atom in residue.atoms():
                        residue_dict['atoms'].append({
                            'index': int(atom.index),
                            'name': str(atom.name),
                            'element': str(atom.element.symbol)
                        })
                    chain_dict['residues'].append(residue_dict)
                topology_dict['chains'].append(chain_dict)

            for atom1, atom2 in topology_object.bonds():
                topology_dict['bonds'].append([
                    int(atom1.index),
                    int(atom2.index)
                ])

        except AttributeError as e:
            raise AttributeError('topology_object fails to implement the'
                'chains() -> residue() -> atoms() and bond() protocol. '
                'Specifically, we encountered the following %s' % e)

        # actually set the attribute
        self._handle.root._v_attrs.topology = json.dumps(topology_dict)

    #####################################################
    # randomState global attribute (optional)
    #####################################################

    @property
    def randomState(self):
        if hasattr(self._handle.root._v_attrs, 'randomState'):
            return self._handle.root._v_attrs.randomState
        return None

    @randomState.setter
    @ensure_mode('w', 'a')
    def randomState(self, value):
        self._handle.root._v_attrs.randomState = str(value)

    #####################################################
    # forcefield global attribute (optional)
    #####################################################

    @property
    def forcefield(self):
        if hasattr(self._handle.root._v_attrs, 'forcefield'):
            return self._handle.root._v_attrs.forcefield
        return None

    @forcefield.setter
    @ensure_mode('w', 'a')
    def forcefield(self, value):
        self._handle.root._v_attrs.forcefield = str(value)

    #####################################################
    # reference global attribute (optional)
    #####################################################

    @property
    def reference(self):
        if hasattr(self._handle.root._v_attrs, 'reference'):
            return self._handle.root._v_attrs.reference
        return None

    @reference.setter
    @ensure_mode('w', 'a')
    def reference(self, value):
        self._handle.root._v_attrs.reference = str(value)

    #####################################################
    # constraints array (optional)
    #####################################################

    @property
    def constraints(self):
        if hasattr(self._handle.root, 'constraints'):
            return self._handle.root.constraints[:]
        return None

    @constraints.setter
    @ensure_mode('w', 'a')
    def constraints(self, value):
        dtype = np.dtype([
                ('atom1', np.int32),
                ('atom2', np.int32),
                ('distance', np.float32)
        ])
        if not value.dtype == dtype:
            raise ValueError('Constraints must be an array with dtype=%s. '
                             'currently, I don\'t do any casting' % dtype)

        if not hasattr(self._handle.root, 'constraints'):
            self._handle.createTable(where='/', name='constraints',
                                     description=dtype)

        self._handle.root.constraints[:] = value

    #####################################################
    # read/write methods for file-like behavior
    #####################################################

    @ensure_mode('r')
    def read(self, n_frames=None, stride=None, atom_indices=None):
        raise NotImplementedError

    @ensure_mode('w', 'a')
    def write(self, coordinates, time=None, cell_lengths=None, cell_angles=None,
                    velocities=None, kineticEnergy=None, potentialEnergy=None,
                    temperature=None, lambdaValues=None):
        """Write one or more frames of data to a HDF5 trajectory file

        This method saves data that is associated with a single simulation
        frame. To save global attributes, use the methods

        Parameters
        ----------
        coordinates : np.ndarray, shape=(n_frames, n_atoms, 3)
        """
        tables = import_tables()  # our delayed import

        # these must be either both present or both absent. since
        # we're going to throw an error if one is present w/o the other,
        # lets do it now.
        if cell_lengths is None and cell_angles is not None:
            raise ValueError('cell_lengths were given, but no cell_angles')
        if cell_lengths is not None and cell_angles is None:
            raise ValueError('cell_angles were given, but no cell_lengths')

        # if the input arrays are simtk.unit.Quantities, convert them
        # into md units. Note that this acts as a no-op if the user doesn't
        # have simtk.unit installed (e.g. they didn't install OpenMM)
        coordinates = in_units_of(coordinates, 'nanometers')
        time = in_units_of(time, 'picoseconds')
        cell_lengths = in_units_of(cell_lengths, 'nanometers')
        cell_angles = in_units_of(cell_angles, 'degrees')
        velocoties = in_units_of(velocoties, 'nanometers/picosecond')
        kineticEnergy = in_units_of(kineticEnergy, 'kilojoules_per_mole')
        potentialEnergy = in_units_of(potentialEnergy, 'kilojoules_per_mole')
        temperature = in_units_of(temperature, 'kelvin')
        lambdaValues = in_units_of(lambdaValues, 'dimensionless')

        # do typechecking and shapechecking on the arrays
        # this ensure_type method has a lot of options, but basically it lets
        # us validate most aspects of the array. Also, we can upconvert
        # on defficent ndim, which means that if the user sends in a single
        # frame of data (i.e. coordinates is shape=(n_atoms, 3)), we can
        # realize that. obviously the default mode is that they want to
        # write multiple frames at a time, so the coordinate shape is
        # (n_frames, n_atoms, 3)
        coordinates = ensure_type(coordinates, dtype=np.float32, ndim=3,
            name='coordinates', shape=(None, None, 3), can_be_none=False,
            warn_on_cast=True, add_newaxis_on_deficient_ndim=True)
        n_frames, n_atoms, = coordinates.shape[0:2]
        time = ensure_type(time, dtype=np.float32, ndim=1,
            name='time', shape=(n_frames,), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, dtype=np.float32, ndim=2,
            name='cell_lengths', shape=(n_frames, 3), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        cell_angles = ensure_type(cell_angles, dtype=np.float32, ndim=1,
            name='cell_angles', shape=(n_frames, 3), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        velocoties = ensure_type(velocoties, dtype=np.float32, ndim=3,
            name='velocoties', shape=(n_frames, n_atoms, 3), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        kineticEnergy = ensure_type(kineticEnergy, dtype=np.float32, ndim=1,
            name='kineticEnergy', shape=(n_frames,), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        potentialEnergy = ensure_type(potentialEnergy, dtype=np.float32, ndim=1,
            name='potentialEnergy', shape=(n_frames,), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        temperature = ensure_type(temperature, dtype=np.float32, ndim=1,
            name='temperature', shape=(n_frames,), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)
        lambdaValues = ensure_type(lambdaValues, dtype=np.float32, ndim=1,
            name='lambdaValues', shape=(n_frames,), can_be_none=True,
            warn_on_cast=False, add_newaxis_on_deficient_ndim=True)

        # if this is our first call to write(), we need to create the headers
        # and the arrays in the underlying HDF5 file
        if self._needs_initialization:
            self.__initialize_headers(
                n_atoms=n_atoms,
                set_coordinates=True,
                set_time=(time is not None),
                set_cell=(cell_lengths is not None or cell_angles is not None),
                set_velocities=(velocities is not None),
                set_kineticEnergy=(kineticEnergy is not None),
                set_potentialEnergy=(potentialEnergy is not None),
                set_temperature=(temperature is not None),
                set_lambdaValues=(lambdaValues is not None))
            self._needs_initialization = False

            # we need to check that that the entries that the user is trying
            # to save are actually fields in OUR file

        try:
            # try to get the nodes for all of the fields that we have
            # which are not None
            for name in ['coordinates', 'time', 'cell_angles', 'cell_lengths',
                         'velocities', 'kineticEnergy', 'potentialEnergy', 'temperature']:
                contents = locals([name])
                if contents is not None:
                    self._handle.getNode(where='/', name=name).append(contents)

            # lambda is different, since the name in the file is lambda
            # but the name in this python function is lambdaValues
            if lambdaValues is not None:
                name = 'lambda'
                self._handle.getNode(where='/', name=name).append(lambdaValues)

        except tables.NoSuchNodeError:
            raise ValueError("The file that you're trying to save to doesn't "
                "contain the field %s. You can always save a new trajectory "
                "and have it contain this information, but I don't allow 'ragged' "
                "arrays. If one frame is going to have %s information, then I expect "
                "all of them to. So I can't save it for just these frames. Sorry "
                "about that :)" % (name, name))

        self._frame_index += n_frames
        self.flush()

    def _initialize_headers(self, n_atoms, set_coordinates, set_time, set_cell,
                            set_velocoties, set_kineticEnergy, set_potentialEnergy,
                            set_temperature, set_lambdaValues):
        tables = import_tables()
        self._n_atoms = n_atoms

        self._handle.root._v_attrs.conventions = 'Pande'
        self._handle.root._v_attrs.conventionVersion = '1.0'
        self._handle.root._v_attrs.program = 'MDTraj'
        self._handle.root._v_atttrs.programVersion = version.short_version
        self._handle.root._v_attrs.title = 'title'

        # if the client has not the title attribute themselves, we'll
        # set it to MDTraj as a default option.
        if not hasattr(self._handle.root._v_attrs, 'application'):
            self._handle.root._v_attrs.application = 'MDTraj'

        # create arrays that store frame level informat
        if set_coordinates:
            self._handle.createEArray(where='/', name='coordinates',
                atom=tables.Float32Atom(), shape=(0, self._n_atoms, 3))
            self._handle.root.coordinates.attrs['units'] = 'nanometers'

        if set_time:
            self._handle.createEArray(where='/', name='time',
                atom=tables.Float32Atom(), shape=(0))
            self._handle.root.time.attrs['units'] = 'picoseconds'

        if set_cell:
            self._handle.createEArray(where='/', name='cell_lengths',
                atom=tables.Float32Atom(), shape=(0, 3))
            self._handle.createEArray(where='/', name='cell_angles',
                atom=tables.Float32Atom(), shape=(0, 3))
            self._handle.root.cell_lengths.attrs['units'] = 'nanometers'
            self._handle.root.cell_angles.attrs['units'] = 'nanometers'

        if set_velocoties:
            self._handle.createEArray(where='/', name='velocities',
                atom=tables.Float32Atom(), shape=(0, self._n_atoms, 3))
            self._handle.root.velocities.attrs['units'] = 'nanometers/picosecond'

        if set_kineticEnergy:
            self._handle.createEArray(where='/', name='kineticEnergy',
                atom=tables.Float32Atom(), shape=(0,))
            self._handle.root.kineticEnergy.attrs['units'] = 'kJ/mol'

        if set_potentialEnergy:
            self._handle.createEArray(where='/', name='potentialEnergy',
                atom=tables.Float32Atom(), shape=(0,))
            self._handle.root.potentialEnergy.attrs['units'] = 'kJ/mol'

        if set_temperature:
            self._handle.createEArray(where='/', name='temperature',
                atom=tables.Float32Atom(), shape=(0,))
            self._handle.root.temperature.attrs['units'] = 'Kelvin'

        if set_lambdaValues:
            self._handle.createEArray(where='/', name='lambda',
                atom=tables.Float32Atom(), shape=(0,))
            self._handle.getNode('/', name='lambda').attrs['units'] = 'dimensionless'

    def close(self):
        "Close the file handle"
        if self._open:
            self._handle.close()
            self._open = False

    def flush(self):
        "Flush all information to disk"
        if self._open:
            self._handle.flush()

    def __del__(self):
        self.close()

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()
