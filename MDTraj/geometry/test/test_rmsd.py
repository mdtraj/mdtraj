# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj.geometry import alignment, rmsd
from mdtraj import IRMSD
import numpy as np

xyz1 = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1.0]])

offset = 1.0 * np.random.normal(size=(3))
rotation = np.array([[1,0,0],[0,0,-1],[0,1,0]])
xyz2 = rotation.dot(xyz1.T).T + offset
xyz3 = xyz2 + np.random.normal(size=xyz1.shape)

def center(xyz_traj):
    for x in xyz_traj:
        x -= (x.astype('float64').mean(0))

def test_irmsd_against_numpy_qcp():
    rmsd_qcp = alignment.rmsd_qcp(xyz1, xyz3)
    
    n_atoms = len(xyz1)
    
    xyz1_traj = np.array([xyz1])
    xyz3_traj = np.array([xyz3])
    center(xyz1_traj)
    center(xyz3_traj)
    
    xyz1_irmsd, n_atoms_padded = rmsd.reshape_irmsd(xyz1_traj)
    xyz3_irmsd, n_atoms_padded = rmsd.reshape_irmsd(xyz3_traj)
    
    G1 = rmsd.calculate_G(xyz1_traj)
    G3 = rmsd.calculate_G(xyz3_traj)
    
    rmsd_IRMSD = IRMSD.rmsd_one_to_all(xyz1_irmsd, xyz3_irmsd, G1, G3, n_atoms, 0).astype('float64')[0]
    
    eq(rmsd_IRMSD, rmsd_qcp, decimal=5)
    

