.. _analysis:

.. currentmodule:: mdtraj

Analysis Reference
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
    rmsf
    lprmsd
    Trajectory.superpose


Hydrogen Bonding
----------------
.. autosummary::
    :toctree: api/generated/

    baker_hubbard
    kabsch_sander
    wernet_nilsson

Pi-Stacking
-----------
.. autosummary::
    :toctree: api/generated/

    pi_stacking


Secondary Structure
-------------------
.. autosummary::
    :toctree: api/generated/

    compute_dssp


Shape Metrics
-------------------
.. autosummary::
    :toctree: api/generated/

    compute_gyration_tensor
    principal_moments
    asphericity
    acylindricity
    relative_shape_antisotropy



Surface Area, Radius of Gyration and Inertia
--------------------------------------------
.. autosummary::
    :toctree: api/generated/

    shrake_rupley
    compute_rg
    compute_inertia_tensor

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

    compute_J3_HN_C
    compute_J3_HN_CB
    compute_J3_HN_HA
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
    compute_directors
