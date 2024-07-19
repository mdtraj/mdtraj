##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################


"""The mdtraj package contains tools for loading and saving molecular dynamics
trajectories in a variety of formats, including Gromacs XTC & TRR, CHARMM/NAMD
DCD, PDB, and HDF5.
"""

# silence cython related numpy warnings, see github.com/numpy/numpy/pull/432
import numpy as _  # noqa
from mdtraj.core import element
from mdtraj.core.topology import Amide, Aromatic, Double, Single, Topology, Triple
from mdtraj.core.trajectory import (
    Trajectory,
    iterload,
    join,
    load,
    load_frame,
    load_topology,
    open,
)
from mdtraj.formats.amberrst import load_ncrestrt, load_restrt
from mdtraj.formats.arc import load_arc
from mdtraj.formats.dcd import load_dcd
from mdtraj.formats.dtr import load_dtr, load_stk
from mdtraj.formats.hdf5 import load_hdf5
from mdtraj.formats.hoomdxml import load_hoomdxml
from mdtraj.formats.lammpstrj import load_lammpstrj
from mdtraj.formats.lh5 import load_lh5
from mdtraj.formats.mdcrd import load_mdcrd
from mdtraj.formats.mol2 import load_mol2
from mdtraj.formats.netcdf import load_netcdf
from mdtraj.formats.openmmxml import load_xml
from mdtraj.formats.pdb import load_pdb
from mdtraj.formats.prmtop import load_prmtop
from mdtraj.formats.psf import load_psf
from mdtraj.formats.registry import FormatRegistry
from mdtraj.formats.trr import load_trr
from mdtraj.formats.xtc import load_xtc
from mdtraj.formats.xyzfile import load_xyz
from mdtraj.geometry import (
    acylindricity,
    asphericity,
    baker_hubbard,
    compute_angles,
    compute_center_of_geometry,
    compute_center_of_mass,
    compute_chi1,
    compute_chi2,
    compute_chi3,
    compute_chi4,
    compute_chi5,
    compute_contacts,
    compute_dihedrals,
    compute_directors,
    compute_displacements,
    compute_distances,
    compute_distances_t,
    compute_drid,
    compute_dssp,
    compute_gyration_tensor,
    compute_inertia_tensor,
    compute_neighborlist,
    compute_neighbors,
    compute_nematic_order,
    compute_omega,
    compute_phi,
    compute_psi,
    compute_rdf,
    compute_rdf_t,
    compute_rg,
    density,
    dipole_moments,
    find_closest_contact,
    isothermal_compressability_kappa_T,
    kabsch_sander,
    principal_moments,
    relative_shape_antisotropy,
    shrake_rupley,
    static_dielectric,
    thermal_expansion_alpha_P,
    wernet_nilsson,
)
from mdtraj.nmr import (
    chemical_shifts_ppm,
    chemical_shifts_shiftx2,
    chemical_shifts_spartaplus,
    compute_chemical_shifts,
    compute_J3_HN_C,
    compute_J3_HN_CB,
    compute_J3_HN_HA,
    reindex_dataframe_by_atoms,
)

from . import reporters
from ._lprmsd import lprmsd
from ._rmsd import rmsd, rmsf

from . import _version

__version__ = _version.get_versions()["version"]

__all__ = (
    "reporters",
    "lprmsd",
    "rmsd",
    "rmsf",
    "element",
    "Amide",
    "Aromatic",
    "Double",
    "Single",
    "Topology",
    "Triple",
    "load_ncrestrt",
    "load_restrt",
    "load_arc",
    "load_dcd",
    "load_dtr",
    "load_stk",
    "load_hdf5",
    "load_hoomdxml",
    "load_lammpstrj",
    "load_lh5",
    "load_mdcrd",
    "load_mol2",
    "load_netcdf",
    "load_xml",
    "load_pdb",
    "load_prmtop",
    "load_psf",
    "load_trr",
    "load_xtc",
    "load_xyz",
    "FormatRegistry",
    "open",
    "load",
    "iterload",
    "load_frame",
    "load_topology",
    "join",
    "Trajectory",
    "baker_hubbard",
    "shrake_rupley",
    "kabsch_sander",
    "compute_distances",
    "compute_distances_t",
    "compute_displacements",
    "compute_angles",
    "compute_dihedrals",
    "compute_phi",
    "compute_psi",
    "compute_chi1",
    "compute_chi2",
    "compute_chi3",
    "compute_chi4",
    "compute_chi5",
    "compute_omega",
    "compute_rg",
    "compute_contacts",
    "compute_drid",
    "compute_center_of_mass",
    "compute_center_of_geometry",
    "wernet_nilsson",
    "compute_dssp",
    "compute_neighbors",
    "compute_neighborlist",
    "compute_rdf",
    "compute_rdf_t",
    "compute_nematic_order",
    "compute_inertia_tensor",
    "compute_gyration_tensor",
    "find_closest_contact",
    "compute_directors",
    "principal_moments",
    "asphericity",
    "acylindricity",
    "relative_shape_antisotropy",
    "dipole_moments",
    "static_dielectric",
    "isothermal_compressability_kappa_T",
    "thermal_expansion_alpha_P",
    "density",
    "compute_J3_HN_C",
    "compute_J3_HN_CB",
    "compute_J3_HN_HA",
    "chemical_shifts_ppm",
    "chemical_shifts_shiftx2",
    "chemical_shifts_spartaplus",
    "compute_chemical_shifts",
    "reindex_dataframe_by_atoms",
)


def capi():
    import os
    import sys

    module_path = sys.modules["mdtraj"].__path__[0]
    return {
        "lib_dir": os.path.join(module_path, "core", "lib"),
        "include_dir": os.path.join(module_path, "core", "lib"),
    }
