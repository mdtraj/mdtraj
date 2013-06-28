Loading Trajectories
====================

MDTraj supports loading (and saving) molecular dynamics trajectories in a
variety of formats. Each of these functions loads a file from disk, returning
a ``Trajectory`` object. Low level details like unit conversion into MDTraj 
standard units (nanometers, degrees, picoseconds) are handled for you.

Only the PDB and HDF5 formats contains the topology information in the trajectory file on disk. Thus to load a trajectory in another format, it's obligatory that you also supply information from which a topology can be extracted, such as a corresponding PDB file.

For low-level access to the files on disk, without the amenities like unit conversion, see the ``XXXTrajectoryFile`` classes, e.g. :ref:`HDF5TrajectoryFile`, :ref:`PDBTrajectoryFile`,  :ref:`XTCTrajectoryFile`, :ref:`DCDTrajectoryFile`, :ref:`TRRTrajectoryFile`, :ref:`NetCDFTrajectoryFile`, :ref:`BINPOSTrajectoryFile`, 


.. currentmodule:: mdtraj

.. autofunction:: load

.. autofunction:: load_pdb

.. autofunction:: load_hdf5

.. autofunction:: load_xtc

.. autofunction:: load_trr

.. autofunction:: load_dcd

.. autofunction:: load_netcdf

.. autofunction:: load_binpos