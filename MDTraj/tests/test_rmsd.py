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


import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq
# from mdtraj.rmsd import RMSDCache, align_array
# from mdtraj._rmsd import getMultipleAlignDisplaceRMSDs_atom_major
from mdtraj.geometry.alignment import rmsd_qcp, compute_translation_and_rotation


# def test_axis_major():
#     t = md.load(get_fn('traj.h5'))
#     pt = md.rmsd_cache(t, major='axis')
#     calculated = pt.rmsds_to(pt, 10)
# 
#     reference = np.zeros(t.n_frames)
#     for i in range(t.n_frames):
#         reference[i] = rmsd_qcp(t.xyz[10], t.xyz[i])
# 
#     eq(calculated, reference, decimal=2)
# 
# 
# def test_atom_major():
#     t = md.load(get_fn('traj.h5'))
#     pt = md.rmsd_cache(t, major='atom')
#     calculated = pt.rmsds_to(pt, 10)
# 
#     reference = np.zeros(t.n_frames)
#     for i in range(t.n_frames):
#         reference[i] = rmsd_qcp(t.xyz[10], t.xyz[i])
# 
#     eq(calculated, reference, decimal=3)
# 

def test_trajectory_rmsd():
    t = md.load(get_fn('traj.h5'))[0:10]
    t.restrict_atoms(np.arange(16))
    calculated = t.rmsd(t, 0, parallel=False)
    
    pt = md.rmsd_cache(t, major='atom')
    reference = np.zeros(t.n_frames)
    for i in range(t.n_frames):
        reference[i] = rmsd_qcp(t.xyz[0], t.xyz[i])

    eq(calculated, reference, decimal=3)
    
# 
# def test_rmsd_to_self():
#     n_atoms = 16
#     np.random.seed(42)
#     conf_atom_1 = np.array(np.random.randn(1, n_atoms, 3), dtype=np.float32)
#     conf_atom_2 = np.copy(conf_atom_1)
# 
# 
#     r_atom_1 = RMSDCache(align_array(conf_atom_1, 'atom'), major='atom', n_atoms=n_atoms)
#     r_atom_2 = RMSDCache(align_array(conf_atom_2, 'atom'), major='atom', n_atoms=n_atoms)
# 
#     yield lambda: eq(float(r_atom_1.rmsds_to(r_atom_2, 0)[0]), 0.0, decimal=2)
# 
#     conf_axis_1 = np.copy(conf_atom_1[0].T.reshape(1, 3, n_atoms))
#     conf_axis_2 = np.copy(conf_axis_1)
# 
#     r_axis_1 = RMSDCache(align_array(conf_axis_1, 'axis'), major='axis', n_atoms=n_atoms)
#     r_axis_2 = RMSDCache(align_array(conf_axis_2, 'axis'), major='axis', n_atoms=n_atoms)
# 
#     yield lambda: eq(float(r_axis_1.rmsds_to(r_axis_2, 0)[0]), 0.0, decimal=2)
# 
# 
# def test_align_displace():
#     t = md.load(get_fn('traj.h5'))
#     t.center_coordinates()
#     pt = md.rmsd_cache(t, major='atom')
#     rmsd, rot = getMultipleAlignDisplaceRMSDs_atom_major(pt.cords, pt.cords, pt._traces(), pt._traces(), pt.cords, pt.cords, t.n_atoms, t.n_atoms, 0)
# 
#     for i in range(t.n_frames):
#         translation, rotation = compute_translation_and_rotation(t.xyz[0], t.xyz[i])
#         eq(rot[i], rotation)
#         eq(float(rmsd_qcp(t.xyz[0], t.xyz[i])), float(rmsd[i]), decimal=3)
