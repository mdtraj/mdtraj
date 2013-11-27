Geometry : :class:`mdtraj.geometry`
===================================

:class:`mdtraj.geometry` contains functions to extract geometry observables
from molecular dynamics trajectories, including distances, angles,
dihedrals and solvent accessible surface areas. Also included is
code for global structural properties, like the radius of gyration,
3d alignment and rotation. Finally, code for more complex structural analysis,
including the identification of hydrogen bonds based on the distance and angle
constraints or Kabsch-Sander hydrogen bond energy (used by DSSP) is included.

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

Hydrogen Bonds
--------------
.. autofunction:: mdtraj.geometry.kabsch_sander
.. autofunction:: mdtraj.geometry.baker_hubbard


Solvent Accessible Surface Area
-------------------------------
.. autofunction:: mdtraj.geometry.shrake_rupley
