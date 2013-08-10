# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MDTraj is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

##############################################################################
# Imports
##############################################################################

import copy
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _HAVE_OPT
if _HAVE_OPT:
    from mdtraj.geometry import C
    from mdtraj.utils.ffi import cpointer


def kabsch_sander(traj):
    """Compute the Kabsch-Sander hydrogen bond energy between each pair
    of residues in every frame.
    
    Hydrogen bonds are defined using an electrostatic definition, assuming
    partial charges of -0.42 e and +0.20 e to the carbonyl oxygen and amide
    hydrogen respectively, their opposites assigned to the carbonyl carbon
    and amide nitrogen. A hydrogen bond is identified if E in the following
    equation is less than -0.5 kcal/mol:
    
    E = 0.42 * 0.2 * 33.2 kcal/(mol*nm) * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
    
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    
    Returns
    -------
    matrices : list of scipy.sparse.csr_matrix
        The return value is a list of length equal to the number of frames
        in the trajectory. Each element is an n_residues x n_residues sparse
        matrix, where the existance of an entry at row `i`, column `j` with value
        `x` means that there exists a hydrogen bond between a backbone CO
        group at residue `i` with a backbone NH group at residue `j` whose
        Kabsch-Sander energy is less than -0.5 kcal/mol (the threshold for
        existance of the "bond"). The exact value of the energy is given by the
        value `x`.
        
    References
    ----------
    .. [1] Kabsch W, Sander C (1983). "Dictionary of protein secondary
    structure: pattern recognition of hydrogen-bonded and geometrical
    features". Biopolymers 22 (12): 2577-637. dio:10.1002/bip.360221211
    """
    if traj.topology is None:
        raise ValueError('kabsch_sander requires topology')
    if not _HAVE_OPT:
        raise RuntimeError('kabsch_sander implementation not available.')
    import scipy.sparse
    
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz',
                      shape=(None, None, 3), warn_on_cast=False)
    
    # putting -1 for one of the entries in the nhco indices list is a
    # signal that the residue does not have that atom type.
    get = lambda l: l[0] if len(l) > 0 else -1
    ca_indices, nhco_indices = [], []
    for residue in traj.topology.residues:
        ca_indices.append(get([a.index for a in residue.atoms if a.name == 'CA']))
        n = get([a.index for a in residue.atoms if a.name == 'N'])
        h = get([a.index for a in residue.atoms if a.name == 'H'])
        c = get([a.index for a in residue.atoms if a.name == 'C'])
        o = get([a.index for a in residue.atoms if a.name == 'O'])
        nhco_indices.append([n, h, c, o])

    nhco_indices = np.array(nhco_indices, np.int32)
    ca_indices = np.array(ca_indices, np.int32)
    n_residues = len(ca_indices)
    hbonds = np.empty((xyz.shape[0], n_residues, 2), np.int32)
    henergies = np.empty((xyz.shape[0], n_residues, 2), np.float32)
    hbonds.fill(-1)
    henergies.fill(np.nan)

    C.kabsch_sander(cpointer(xyz), cpointer(nhco_indices), cpointer(ca_indices),
                    xyz.shape[0], xyz.shape[1], n_residues,
                    cpointer(hbonds), cpointer(henergies))

    # The C code returns its info in a pretty inconvenient format.
    # Let's change it to a list of scipy CSR matrices.

    matrices = []
    hbonds_mask = (hbonds != -1)
    for i in range(xyz.shape[0]):
        # appologies for this cryptic code -- we need to deal with the low
        # level aspects of the csr matrix format.
        hbonds_frame = hbonds[i]
        mask = hbonds_mask[i]
        henergies_frame = henergies[i]
        
        indptr = np.zeros(n_residues + 1, np.int32)
        indptr[1:] = np.cumsum(mask.sum(axis=1))
        indices = hbonds_frame[mask].flatten()
        data = henergies_frame[mask].flatten()

        matrices.append(scipy.sparse.csr_matrix((data, indices, indptr), shape=(n_residues, n_residues)))
    
    return matrices

# import mdtraj as md
# import matplotlib.pyplot as pp
# t = md.load('MDTraj/testing/reference/2EQQ.pdb')
# #t = md.Trajectory(xyz=np.random.randn(10000,t.n_atoms, 3), topology=t.topology)
# 
# import time
# t0 = time.time()
# a = kabsch_sander(t)
# print a[0].sum(axis=0)
# print a[0].nonzero()
# print time.time() - t0
# #pp.matshow(kabsch_sander(t[0])[0].todense())
# #pp.show()