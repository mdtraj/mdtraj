.. currentmodule:: mdtraj

Trajectory Analysis
===================

Trajectory analysis is the heart of MDTraj. These functions can be used to run
a variety of analyses on :class:`mdtraj.Trajectory` objects.
It's usually as simple as ::

    >>> import mdtraj as md
    >>> t = md.load('trajectory.pdb)
    >>> print(md.compute_phi(t))

Root-mean-square deviation (RMSD)
---------------------------------
.. autosummary::
    :toctree: generated/

    rmsd
    lprmsd
    Trajectory.superpose


Hydrogen Bonding
----------------
.. autosummary::
    :toctree: generated/

    baker_hubbard
    kabsch_sander
    wernet_nilsson


Surface Area and Radius of Gyration
-----------------------------------
.. autosummary::
    :toctree: generated/

    shrake_rupley
    compute_rg


Distances
---------
.. autosummary::
    :toctree: generated/

    compute_distances
    compute_displacements
    compute_contacts
    compute_drid
    compute_center_of_mass
    geometry.squareform

Bond Angles and Dihedrals
-------------------------
.. autosummary::
    :toctree: generated/

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
    :toctree: generated/

    chemical_shifts_shiftx2
    chemical_shifts_ppm
    chemical_shifts_spartaplus
    reindex_dataframe_by_atoms


