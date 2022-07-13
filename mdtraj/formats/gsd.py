##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Alexander Yang
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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division

import numpy as np
import os


from mdtraj.core.topology import Topology
from mdtraj.core.element import virtual_site
from mdtraj.utils import (ensure_type, cast_indices, 
        lengths_and_angles_to_tilt_factors)
from mdtraj.formats.registry import FormatRegistry
from mdtraj.utils.six import string_types

__all__ = ['load_gsd', 'write_gsd', 'load_gsd_topology']


@FormatRegistry.register_loader('.gsd')
def load_gsd(filename, top=None, start=None, n_frames=None, stride=None, 
        atom_indices=None, frame=None):
    """Load a GSD trajectory file.

    Parameters
    -----------
    filename : path-like
        Path of GSD trajectory file.
    top : {path-like, Trajectory, Topology}, None
        A pdb file, a trajectory, or a topology to supply topology information
        If None, topology information will be parsed from the GSD file
    start : int, None
        First frame to convert
    n_frames : int, None
        Number of frames after `start` to convert
    stride : int
        Read only every stride-th frame.   
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file.
    frame : int, optional
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    """
    from mdtraj.core.trajectory import Trajectory, _parse_topology
    import gsd.hoomd

    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_gsd. '
                        'you supplied %s'.format(type(filename)))

    if top is not None:
        topology = _parse_topology(top)
    else:
        topology = load_gsd_topology(filename)
    atom_indices = cast_indices(atom_indices)

    with gsd.hoomd.open(filename, 'rb') as f:
        if frame is not None:
            xyz, vectors, time = read_snapshot(frame, f[frame], 
                    topology, atom_indices=atom_indices)
            t = Trajectory(xyz=np.array(xyz), topology=topology, 
                    time=np.array([time]))
            t.unitcell_vectors = np.reshape(vectors, (-1,3,3))
            return t

        else:
            return hoomdtraj_to_traj(f, topology, start=start, n_frames=n_frames,
                    stride=stride, atom_indices=atom_indices)

def load_gsd_topology(filename, frame=0):
    """ Create an MDTraj.Topology from a GSD file 
    
    Parameters
    ----------
    filename : path-like
        Path of GSD trajectory file.
    frame : int, 0 
        Frame of GSD file to parse topology

    Returns
    -------
    top : mdtraj.Topology

    Notes
    -----
    GSD files support systems with variable topologies.
    For compatibility with MDTraj, only the topology from GSD frame 0 is
    used to construct the MDTraj topology.
    """
    import gsd.hoomd
    with gsd.hoomd.open(filename, 'rb') as gsdfile:
        top = Topology()
        generic_chain = top.add_chain()
        generic_residue = top.add_residue('A', generic_chain)
        all_particle_types = gsdfile[frame].particles.types
        for particle_type_id in gsdfile[frame].particles.typeid:
            top.add_atom(all_particle_types[particle_type_id], virtual_site,
                    generic_residue)

        for bond in gsdfile[frame].bonds.group:
            atom1, atom2 = bond[0], bond[1]
            top.add_bond(top.atom(atom1), top.atom(atom2))

    return top

def hoomdtraj_to_traj(f, topology, start=None, n_frames=None, 
        stride=None, atom_indices=None):
    """ Convert HOOMDTrajectory to MDtraj Trajectory 
   
    Parameters
    ----------
    f : gsd.hoomd.HOOMDTrajectory object
    topology : mdtraj.Topology
    start : int, None
        First frame to convert
    n_frames : int, None
        Number of frames after `start` to convert
    stride : int
        Read only every stride-th frame.   
    atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates
            from the file.

    Returns
    --------
    traj : mdtraj.Trajectory object
    
    """
    from mdtraj.core.trajectory import Trajectory
    if start is None:
        start = 0
    if n_frames is None:
        n_frames = len(f) - start
    if stride is None:
        stride = 1

    all_coords, all_times, all_vectors = [], [], []
    for i, snapshot in enumerate(f[start : start+n_frames : stride], start=start):
        xyz, box_vectors, time = read_snapshot(i, snapshot, topology, 
                atom_indices=atom_indices)
        all_coords.append(xyz)
        all_vectors.append(box_vectors)
        all_times.append(time)

    all_coords = np.array(all_coords)
    all_vectors = np.array(all_vectors)
    all_times = np.array(all_times)
    if len(all_coords) == 0:
        return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), 
                topology=topology)

    t = Trajectory(xyz=all_coords, topology=topology, time=all_times)
    t.unitcell_vectors = all_vectors
    return t

def read_snapshot(frame, snapshot, topology, atom_indices=None):
    """ Parse relevant information from a single HOOMD snapshot (frame) 
    
    Parameters
    ----------
    frame : int
        Frame index to read
    snapshot : gsd.hoomd.Snapshot
    topology : mdtraj.Topology
    atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates
            from the file.

    Returns
    --------
    xyz : np.ndarray, (n, 3)
    box_vectors : list, (3, 3)
    time : int
        Step of hoomd snapshot

    Notes
    -----
    If the GSD file has a variable topology compared to the
    given topology, an IOError will be raised.
    
    """

    _check_topology(frame, snapshot, topology)
    xyz = snapshot.particles.position

    if atom_indices is not None:
        xyz = xyz[atom_indices]
    lx, ly, lz, xy, xz, yz = snapshot.configuration.box
    box_vectors=[[lx, xy*ly, xz*lz],
                            [0.0, ly, yz*lz],
                            [0.0, 0.0,lz]]
    time = snapshot.configuration.step

    return xyz, box_vectors, time

def write_gsd(filename, xyz, top, cell_lengths=None, cell_angles=None):
    """Write one or more frames of data to a gsd file.

    Parameters
    ----------
    filename : path-like
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
        The cartesian coordinates of the atoms to write. 
    top : mdtraj.Topology
    cell_lengths : np.ndarray, dtype=np.double, shape=(n_frames, 3)
        The lengths (a,b,c) of the unit cell for each frame. 
    cell_angles : np.ndarray, dtype=np.double, shape=(n_frames, 3)
        The angles (\alpha, \beta, \gamma) defining the unit cell for
        each frame. (Units of degrees).
    """
    import gsd.hoomd
    xyz = ensure_type(xyz, np.float32, 3, 'xyz', can_be_none=False,
            shape=(None, None, 3), warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True)
    cell_lengths = ensure_type(cell_lengths, np.float32, 2, 'cell_lengths',
            can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True)
    if cell_angles is None:
        cell_angles = np.empty_like(cell_lengths)
        cell_angles.fill(90)
    cell_angles = ensure_type(cell_angles, np.float32, 2, 'cell_angles',
            can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
            add_newaxis_on_deficient_ndim=True)
    unique_types = list(set(a.name for a in top.atoms)) # Find unique atomtypes

    types = {a:i for i, a in enumerate(unique_types)} # Map str(atomtype) to index
    atomtype_ids = np.array([types[a.name] for a in top.atoms], dtype=np.int32)

    unique_bond_types = []
    bondtype_ids = []
    bond_groups = []
    if top.n_bonds > 0:
        unique_bond_types, bondtype_ids, bond_groups = _process_bonds(top)

    with gsd.hoomd.open(filename, 'wb') as hoomd_traj:
        for i, coords in enumerate(xyz):
            gsd_frame = gsd.hoomd.Snapshot()
            gsd_frame.particles.N = top.n_atoms
            gsd_frame.particles.position = coords 
            gsd_frame.particles.types = unique_types
            gsd_frame.particles.typeid = atomtype_ids
            gsd_frame.bonds.N = top.n_bonds
            gsd_frame.bonds.types = unique_bond_types
            gsd_frame.bonds.typeid = bondtype_ids
            gsd_frame.bonds.group = bond_groups

            gsd_frame.configuration.box = lengths_and_angles_to_tilt_factors(
                    cell_lengths[i][0], cell_lengths[i][1], cell_lengths[i][2],
                    cell_angles[i][0], cell_angles[i][1], cell_angles[i][2])
            hoomd_traj.append(gsd_frame)

def _process_bonds(top):
    """ Identify and relate unique bondtypes to an index """
    sorted_bond_types = [tuple(sorted([b[0].name, b[1].name])) for b in top.bonds]
    unique_bond_types = list(set("{}-{}".format(b[0], b[1]) 
            for b in sorted_bond_types))
    bondtypes = {a:i for i, a in enumerate(unique_bond_types)}
    bondtype_ids = [bondtypes["{}-{}".format(key[0], key[1])] 
        for key in sorted_bond_types]
    bond_groups = [(b[0].index, b[1].index) for b in top.bonds]

    return (unique_bond_types, bondtype_ids, bond_groups)

def _check_topology(frame, snapshot, topology):
    """ Verify snapshot topology matches given MDTraj topology 
    
    Notes
    -----
    Only looks for N particles and N bonds for speed purposes"""
    error_msg = ("GSD frame {} ".format(frame) + 
            "has inconsistent topology compared to given topology, " + 
            "this is unsupported in MDTraj.")
    if (snapshot.particles.N != topology.n_atoms or
            snapshot.bonds.N != topology.n_bonds):
        raise IOError(error_msg)
