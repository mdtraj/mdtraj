.. currentmodule:: mdtraj.reporters

OpenMM Reporters : :class:`md.reporters`
========================================

MDTraj provides a few flexible reporters for use with the `OpenMM <https://simtk.org/home/openmm>`_ python app. OpenMM is a toolkit for molecular simulation using high performance GPU code. OpenMM itself ships with a `DCD reporter <https://simtk.org/api_docs/openmm/api5_0/python/classsimtk_1_1openmm_1_1app_1_1dcdreporter_1_1DCDReporter.html>`_, but it lacks the ability to, for instance, report on only a subset of the atoms, which might be desired to print only the protein coordinates and discard water during a simulation.

MDTraj currently provides three reporters, ``HDF5Reporter``, ``NetCDFReporter`` and ``DCDReporter``. Of these, ```HDF5Reporter`` is the most flexible, because the :ref:`HDF5 Format <HDF5FormatSpec>` is the most full-featured trajectory file format available.


Example Usage
-------------

:: 
 
 >>> # todo: add more examples of how to use these classes


Classes
-------


.. _HDF5Reporter:
.. autoclass:: HDF5Reporter

.. _NetCDFReporter:
.. autoclass:: NetCDFReporter

.. _DCDReporter:
.. autoclass:: DCDReporter

