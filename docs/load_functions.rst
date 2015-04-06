.. _loading:
Trajectories
------------

Loading molecular dynamics trajectories is easy with MDTraj. ::

    >>> import mdtraj as md
    >>> t = md.load('trajectory.pdb')

The function :func:`mdtraj.load` will automatically detect the appropriate format
based on the filename extension. The supported formats are all listed below in
the format-specific loaders. The return value of :func:`mdtraj.load` is an
:class:`mdtraj.Trajectory` object.

While :func:`mdtraj.load` will read an entire trajectory into memory, other
functions like :func:`mdtraj.load` and :func:`mdtraj.open` are available for
working with trajectories in chunks, without loading them entirely into memory
all at once.

The trajectory object
*********************

.. currentmodule:: mdtraj
.. autosummary::
    :toctree: api/generated/

    Trajectory
    Topology

Cross-format loaders
********************

.. currentmodule:: mdtraj
.. autosummary::
    :toctree: api/generated/

    load
    iterload
    load_frame
    open

Format-specific loaders
***********************

.. autosummary::
    :toctree: api/generated/

    load_binpos
    load_lh5
    load_pdb
    load_xml
    load_arc
    load_dcd
    load_hdf5
    load_netcdf
    load_trr
    load_xtc
    load_prmtop
    load_xyz
    load_lammpstrj
    load_hoomdxml
