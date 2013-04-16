"""OpenMM Reporter for saving the positions of a molecular dynamics simulation
in the HDF5 format.
"""

##############################################################################
# Imports
##############################################################################
# stdlib
import os
import warnings

# 3rd party
import tables
import numpy as np

try:
    # openmm
    import simtk.unit as units
    OPENMM_IMPORTED = True
except ImportError:
    # if someone tries to import all of mdtraj but doesn't
    # OpenMM installed, we don't want that to choke. It should
    # only choke if they actually try to USE the HDF5Reporter
    OPENMM_IMPORTED = False

# local
import mdtraj.topology
import mdtraj.trajectory
import mdtraj.io


##############################################################################
# Classes
##############################################################################


class HDF5Reporter(object):
    """HDF5Reporter stores a molecular dynamics trajectory in the HDF5 format.
    The atomic positions, periodic box vectors, and time index are saved.

    Example
    -------
    >>> simulation = Simulation(topology, system, integrator)
    >>> h5_reporter = HDF5Reporter('traj.lh5', 100)
    >>> simulation.reporters.append(h5_reporter)
    >>> simulation.step(10000)

    >>> traj = mdtraj.load('traj.lh5')
    """

    def __init__(self, file, reportInterval, lossy=True, **kwargs):
        """Create a LH5Reporter.

        Parameters
        ----------
        file : string
            The name of the file to save. If using lossy encoding, the '.lh5'
            extension is recommended. If not using lossy encoding, the '.hdf5'
            or '.h5' extensions are recommended.
        reportInterval : int
            The interval (in time steps) at which to write frames
        lossy : bool
            Use a fixed-width compression to store each coordinate with
            16 bits on disk. This is the same compression algorithm used by
            the gromacs XTC file format.

        Other Parameters
        -----------------
        append : bool, default=False
            If you'd like to append to an existing trajectory, you must
            explicitly activate this flag. To use this safely, you need to
            know what you're doing.
        filter : tables.Filter, or None
            A compression filter used to minimize storage. You can supply
            your own filter if you would like. See the pytables documentation
            (http://pytables.github.io/usersguide/libref/helper_classes.html#the-filters-class)
            for details on how to construct a filter. By default, we use the
            BLOSC filter (level 9) with shuffling enabled. If you'd like to do
            analysis of your trajectory using other HDF5-enabled programs that
            do not link against BLOSC, you might consider using a different
            compression library, or using None.
        n_expected_frames : int, default=10000
            An estimate about the number of frames that will be added to
            the trajectory. If you plan to create either a much smaller or a
            much bigger trajectory, try providing a guess; this will optimize
            the HDF5 B-Tree creation and management process time and the
            amount of memory used.
        """
        append = kwargs.pop('append', False)
        if os.path.exists(file) and not append:
            raise OSError('The file specified, "%s", already exists. If you\d '
                          'like to append to it, you need to explicitly use '
                          'the append_to_file keyword argument' % file)
        self._handle = tables.openFile(file, 'a')
        if (file.endswith('.h5') or file.endswith('.hdf5')) and lossy:
            warnings.warn('''Extension warning: The .h5 extension is generally
            used for trajectories *without* the lossy encoding. Consider using
            .lh5 instead, if using the lossy encoding. ''')
        if (file.endswith('.lh5') or file.endswith('.lhdf5')) and (not lossy):
            warnings.warn('''Extension warning: The .lh5 extension is generally
            used for trajectories *with* the lossy encoding. Consider using
            the .h5 extension instead, if you prefer not to use the lossy
            encoding''')

        self._reportInterval = reportInterval
        self._lossy = True
        self._filter = kwargs.pop('filter', mdtraj.io.COMPRESSION)
        self._n_expected_frames = kwargs.pop('n_expected_frames', 10000)
        self._is_intialized = False
        self._n_particles = None

        if len(kwargs) > 0:
            raise KeyError('Unrecognized keyword arg %s' % kwargs.keys()[0])

        if not OPENMM_IMPORTED:
            raise ImportError('OpenMM not found.')

    def _initialize(self, simulation):
        """Deferred initialization of the reporter, which happens before
        processing the first report.

        At the time that the first report is processed, we now have access
        to the simulation object, which we don't have at the point when the
        reporter is instantiated

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
            The Simulation to generate a report for
        """
        ######################################################################
        # "/xyz" node, to store the coordinates
        ######################################################################
        if self._lossy:
            xyz_atom = tables.Int16Atom()
        else:
            xyz_atom = tables.Float32Atom()

        self._n_particles = simulation.system.getNumParticles()
        try:
            xyz_node = self._handle.getNode(where='/', name='xyz')
            if xyz_node.shape[1] != self._n_particles:
                raise ValueError('The number of atoms in the trajectory '
                    'doesn\'t match the number of atoms in your system. '
                    'Perhaps you should just save your trajectory to a '
                    'different filename?')
            if xyz_node.atom != xyz_atom:
                raise ValueError('Atom mismatch. Perhaps you should just save '
                    'your trajectory to a different filename?')
        except tables.NoSuchNodeError:
            # if it doesn't contain the xyz node, make it
            self._handle.createEArray(where='/', name='xyz', atom=xyz_atom,
                shape=[0, self._n_particles, 3], filters=self._filter,
                expectedrows=self._n_expected_frames)

        ######################################################################
        # "/topology" node, to store the topology.
        # this one is stored only once -- its not extended during each report.
        ######################################################################
        topology = mdtraj.topology.to_bytearray(simulation.topology)
        try:
            topology_node = self._handle.getNode(where='/', name='topology')
            # if the file does contain a topology, make sure its the same
            if not np.all(topology == topology_node[:]):
                raise ValueError('The topology already in the file does not '
                    'your system\'s topology. You should save your trajectory '
                    'under a unique filename')
        except tables.NoSuchNodeError:
            # if the file doesn't contain a topology, save ours
            topology_node = self._handle.createCArray(where='/', name='topology',
                atom=tables.Atom.from_dtype(topology.dtype), shape=topology.shape,
                filters=self._filter)
            topology_node[:] = topology

        ######################################################################
        # "/time" node, to store the time
        ######################################################################
        try:
            self._handle.getNode(where='/', name='time')
        except tables.NoSuchNodeError:
            self._handle.createEArray(where='/', name='time',
                atom=tables.Float32Atom(), shape=[0], filters=self._filter,
                expectedrows=self._n_expected_frames)

        ######################################################################
        # "/box" node, to store the box vectors
        ######################################################################
        try:
            self._handle.getNode(where='/', name='box')
        except tables.NoSuchNodeError:
            self._handle.createEArray(where='/', name='box',
                atom=tables.Float32Atom(), shape=[0, 3, 3],
                filters=self._filter, expectedrows=self._n_expected_frames)

        self._handle.flush()

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
            The Simulation to generate a report for

        Returns
        -------
        report_description : tuple
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
            The Simulation to generate a report for
        state : simtk.openmm.State
            The current state of the simulation
        """
        if not self._is_intialized:
            self._initialize(simulation)
            self._is_intialized = True

        # add the positions to the file
        xyz = state.getPositions(asNumpy=True).value_in_unit(units.nanometers)
        if self._lossy:
            xyz = mdtraj.trajectory._convert_to_lossy_integers(xyz)
        self._handle.root.xyz.append(xyz.reshape(1, self._n_particles, 3))

        # add the box vectors to the file
        box = state.getPeriodicBoxVectors().value_in_unit(units.nanometers)
        self._handle.root.box.append(np.array(box).reshape(1, 3, 3))

        # add the time to the file
        time = state.getTime().value_in_unit(units.picosecond)
        self._handle.root.time.append([time])

        # flush the file to disk. it might not be necessary to do this every
        # report, but this is the most proactive solution. We don't want to
        # accumulate a lot of data in memory only to find out, at the very
        # end of the run, that there wasn't enough space on disk to hold the
        # data.
        self._handle.flush()

    def __del__(self):
        if hasattr(self, '_handle'):
            self._handle.close()
