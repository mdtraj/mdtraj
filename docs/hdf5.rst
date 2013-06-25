HDF5TrajectoryFile
==================

.. currentmodule:: mdtraj

.. autoclass:: HDF5TrajectoryFile
   :members: flush, close
   
   .. automethod:: read(n_frames=None, stride=None, atom_indices=None)

   .. automethod:: write(coordinates, time=None, cell_lengths=None, cell_angles=None, velocities=None, kineticEnergy=None, potentialEnergy=None, temperature=None, alchemicalLambda=None)