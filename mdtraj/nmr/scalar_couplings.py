#!/usr/bin/env python
# -*- coding: latin-1 -*-
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors: Robert McGibbon, TJ Lane, Osama El-Gabalawy
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
#############################################################################

"""
This file contains scripts for calculating scalar (J) Couplings from backbone dihedrals.
"""

##############################################################################
# Imports
##############################################################################

import numpy as np

from mdtraj.geometry import compute_phi

##############################################################################
# Globals
##############################################################################


J3_HN_HA_coefficients = {  # See full citations below in docstring references.
    "Ruterjans1999": dict(phi0=-60 * np.pi/180., A=7.90, B=-1.05, C=0.65),  # From Table 1. in paper.
    "Bax2007": dict(phi0=-60 * np.pi/180., A=8.4, B=-1.36, C=0.33),  # From Table 1. in paper
    "Bax1997": dict(phi0=-60 * np.pi/180., A=7.09, B=-1.42, C=1.55),  # From Table 2. in paper
    }

J3_HN_HA_uncertainties = {
    # Values in [Hz]
    "Ruterjans1999": 0.25,
    "Bax2007": 0.36,
    "Bax1997": 0.39
}

##############################################################################
# Functions
##############################################################################

def _J3_function(phi, A, B, C, phi0):
    """Return a scalar couplings with a given choice of karplus coefficients.  USES RADIANS!"""
    return A * np.cos(phi + phi0) ** 2. + B * np.cos(phi + phi0) + C


def compute_J3_HN_HA(traj, model="Bax2007"):
    """Calculate the scalar coupling between HN and H_alpha.

    This function does not take into account periodic boundary conditions (it
    will give spurious results if the three atoms which make up any angle jump
    across a PBC (are not "wholed"))

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory to compute J3_HN_HA for
    model : string, optional, default="Bax2007"
        Which scalar coupling model to use.  Must be one of Bax2007, Bax1999,
        or Ruterjans1999

    Returns
    -------
    indices : np.ndarray, shape=(n_phi, 4), dtype=int
        Atom indices (zero-based) of the phi dihedrals
    J : np.ndarray, shape=(n_phi, n_frames)
        Scalar couplings (J3_HN_HA, in [Hz]) of this trajectory.
        `J[k]` corresponds to the phi dihedral associated with
        atoms `indices[k]`

    Notes
    -----
    The coefficients are taken from the references below--please cite them.

    References
    ----------
    .. [1] Schmidt, J. M., Blümel, M., Löhr, F., & Rüterjans, H.
       "Self-consistent 3J coupling analysis for the joint calibration
       of Karplus coefficients and evaluation of torsion angles."
       J. Biomol. NMR, 14, 1 1-12 (1999)

    .. [2] Vögeli, B., Ying, J., Grishaev, A., & Bax, A.
       "Limits on variations in protein backbone dynamics from precise
       measurements of scalar couplings."
       J. Am. Chem. Soc., 129(30), 9377-9385 (2007)

    .. [3] Hu, J. S., & Bax, A.
       "Determination of ϕ and ξ1 Angles in Proteins from 13C-13C
       Three-Bond J Couplings Measured by Three-Dimensional Heteronuclear NMR.
       How Planar Is the Peptide Bond?."
       J. Am. Chem. Soc., 119(27), 6360-6368 (1997)

    """
    indices, phi = compute_phi(traj)

    if model not in J3_HN_HA_coefficients:
        raise(KeyError("model must be one of %s" % J3_HN_HA_coefficients.keys()))

    J = _J3_function(phi, **J3_HN_HA_coefficients[model])
    return indices, J
