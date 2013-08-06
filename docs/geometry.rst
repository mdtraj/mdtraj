Geometry : :class:`mdtraj.geometry`
===================================

:class:`mdtraj.geometry` contains functions to extract geometry observables
from molecular dynamics trajectories, including distances, angles and
dihedrals. Also included is code for global structural properties, like
the radius of gyration, and for 3d alignment and rotation.

The core geometry code for calculating distances, angles and dihedrals
is written in both SIMD optimized C and standard numpy. The implementations
are tested against each other (as well as external references) and the
fastest platform available on your machine will be used.


Distances
---------

.. autofunction:: mdtraj.geometry.compute_distances

Angles
------

.. autofunction:: mdtraj.geometry.compute_angles


Dihedrals
---------

.. automodule:: mdtraj.geometry.dihedral

Radius of gyration
------------------

.. automodule:: mdtraj.geometry.rg

Alignment
---------
Note: this code is not an interface to the IRMSD fast SSE2/SSE3
:ref:`RMSD library <RMSD>`. This code provides more flexibility, at
the cost of some speed (IRMSD can only compute distances).

.. automodule:: mdtraj.geometry.alignment