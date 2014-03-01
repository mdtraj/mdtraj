OpenMM Reporters
================

MDTraj provides a few flexible reporters for use with the `OpenMM <https://simtk.org/home/openmm>`_ python app. OpenMM is a toolkit for molecular simulation using high performance GPU code. OpenMM itself ships with a `DCD reporter <https://simtk.org/api_docs/openmm/api5_0/python/classsimtk_1_1openmm_1_1app_1_1dcdreporter_1_1DCDReporter.html>`_, but it lacks the ability to, for instance, report on only a subset of the atoms, which might be desired to print only the protein coordinates and discard water during a simulation.

MDTraj currently provides three reporters, ``HDF5Reporter``, ``NetCDFReporter`` and ``DCDReporter``. Of these, ```HDF5Reporter`` is the most flexible, because the :ref:`HDF5 Format <HDF5FormatSpec>` is the most full-featured trajectory file format available.


Example Usage
-------------

.. code-block:: python

   from simtk.openmm.app import *
   from simtk.openmm import *
   from simtk.unit import *
   from mdtraj.reporters import NetCDFReporter         # <-- new import from mdtraj

   pdb = PDBFile('input.pdb')
   forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                    nonbondedCutoff=1*nanometer, constraints=HBonds)
   integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
   simulation = Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)
   simulation.minimizeEnergy()
   simulation.reporters.append(NetCDFReporter('output.nc'), 1000)   # <-- AMBER compatible
   simulation.step(10000)

.. currentmodule:: mdtraj.reporters
.. autosummary::
    :toctree: generated/
    
    HDF5Reporter
    NetCDFReporter
    DCDReporter
