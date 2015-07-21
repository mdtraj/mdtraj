.. _analysis:

.. currentmodule:: mdtraj

Analysis Functions
==================

Trajectory analysis is the heart of MDTraj. These functions can be used to run
a variety of analyses on :class:`mdtraj.Trajectory` objects.
It's usually as simple as ::

    >>> import mdtraj as md
    >>> t = md.load('trajectory.pdb')
    >>> print(md.compute_phi(t))

Root-mean-square deviation (RMSD)
---------------------------------
.. autosummary::
    :toctree: api/generated/

    rmsd
    lprmsd
    Trajectory.superpose


Hydrogen Bonding
----------------
.. autosummary::
    :toctree: api/generated/

    baker_hubbard
    kabsch_sander
    wernet_nilsson


Secondary Structure
-------------------
.. autosummary::
    :toctree: api/generated/

    compute_dssp


Surface Area and Radius of Gyration
-----------------------------------
.. autosummary::
    :toctree: api/generated/

    shrake_rupley
    compute_rg


Distances
---------
.. autosummary::
    :toctree: api/generated/

    compute_distances
    compute_displacements
    compute_neighbors
    compute_contacts
    compute_drid
    compute_center_of_mass
    geometry.squareform
    compute_rdf


Bond Angles and Dihedrals
-------------------------
.. autosummary::
    :toctree: api/generated/

    compute_angles
    compute_dihedrals
    compute_phi
    compute_psi
    compute_chi1
    compute_chi2
    compute_chi3
    compute_chi4
    compute_omega



NMR Observables
---------------
.. autosummary::
    :toctree: api/generated/

    chemical_shifts_shiftx2
    chemical_shifts_ppm
    chemical_shifts_spartaplus
    reindex_dataframe_by_atoms


Thermodynamic Quantities
------------------------
.. autosummary::
    :toctree: api/generated/

    dipole_moments
    static_dielectric
    isothermal_compressability_kappa_T
    thermal_expansion_alpha_P
    density

Order Parameters
----------------
.. autosummary::
    :toctree: api/generated/

    compute_nematic_order
