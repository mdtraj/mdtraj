# Copyright 2012-present mdtraj developers
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
"""Compatibility loader for MSMBuilder2 "LH5" trajectory format.

To load .lh5 trajectory files produced by MSMBuilder2, you can import this
module and it will register a loader, making the format available via
``mdtraj.load``.

Examples
--------
>>> import mdtraj as md
>>> md.load('msmbuilder26_trajectoryfile.lh5')                # doctest: +SKIP
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/Library/Frameworks/EPD64.framework/Versions/7.3/lib/python2.7/site-packages/mdtraj-0.2-py2.7-macosx-10.5-x86_64.egg/mdtraj/trajectory.py", line 145, in load
    'with extensions in %s' % (filename, extension, _LoaderRegistry.keys()))
IOError: Sorry, no loader for filename=msmbuilder26_trajectoryfile.lh5 (extension=.lh5) was found. I can only load files with extensions in ['.h5', '.xtc', '.ncdf', '.pdb', '.xml', '.binpos', '.dcd', '.nc', '.trr']

>>> # but let's load the compatibility module. it will register itself
>>> import mdtraj.compatibility
>>> # and now we can load these legacy trajectories without issue
>>> md.load('msmbuilder26_trajectoryfile.lh5')                # doctest: +SKIP
<mdtraj.Trajectory with 501 frames, 22 atoms at 0x1019bcdd0>

Functions
---------
"""
##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import sys
import numpy as np

from mdtraj import Trajectory, Topology
import mdtraj.pdb.element
import mdtraj.trajectory
from mdtraj.utils import import_
from mdtraj.utils.six import PY3
from mdtraj.utils.six.moves import zip
if PY3:
    basestring = str

tables = import_('tables')

##############################################################################
# Globals
##############################################################################

__all__ = ['load_legacy_hdf']

##############################################################################
# Functions
##############################################################################


def load_legacy_hdf(filename, stride=1, frame=None, chunk=50000,
                     upconvert_int16=True):
    """Load a HDF5 file in the legacy MSMBuilder2 lh5 format

    Parameters
    ----------
    filename : str
        String filename of HDF Trajectory file.
    stride : int, default=1
        Only read every stride-th frame
    frame : {None, int}, default=None
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
    chunk : int, default=50000
       Size of the chunk to use for loading the file.
    upconvert_int16 : bool, default=True
        If the data is encoded in the file as int16, an automatic upconversion
        to float32 based on the Gromacs lossy encoding scheme is done.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """

    def _convert_from_lossy_integers(X, precision=1000):
        """Implementation of the lossy compression used in Gromacs XTC using
        the pytables library.  Convert 16 bit integers into 32 bit floats."""
        X2 = X.astype("float32")
        X2 /= float(precision)
        return X2

    if not isinstance(filename, basestring):
        raise TypeError('filename must be of type string for load_hdf5. '
            'you supplied %s' % type(filename))

    F = tables.File(filename, 'r')
    try:
        topology = _topology_from_arrays(F.root.AtomID[:], F.root.AtomNames[:],
                                         F.root.ChainID[:], F.root.ResidueID[:],
                                         F.root.ResidueNames[:])

        # find the xyz coordinates in the file
        if hasattr(F.root, 'XYZList'):
            xyz_node = F.root.XYZList
        else:
            raise ValueError('XYZ data not found in %s' % filename)

        if chunk < stride:
            raise ValueError('Chunk must be greater than or equal to the stride')

        # if they want only a single frame, quit early without
        # doing the chunking
        if frame is not None:

            xyz = xyz_node[int(frame)]
            if upconvert_int16 and xyz.dtype == np.int16:
                xyz = _convert_from_lossy_integers(xyz)

            trajectory = Trajectory(xyz=xyz, topology=topology, time=frame)
            return trajectory

        # adjust the stride/chunk
        if stride is not None:
            # Need to do this in order to make sure we stride correctly.
            # since we read in chunks, and then we need the strides
            # to line up
            while chunk % stride != 0:
                chunk -= 1
        else:
            stride = 1

        shape = xyz_node.shape
        begin_range_list = np.arange(0, shape[0], chunk)
        end_range_list = np.concatenate((begin_range_list[1:], [shape[0]]))

        def enum_chunks():
            for r0, r1 in zip(begin_range_list, end_range_list):
                xyz = np.array(xyz_node[r0: r1: stride])

                yield xyz

        xyz = np.concatenate([c for c in enum_chunks()])
        time = np.arange(len(F.root.XYZList))[::stride]

        if upconvert_int16 and xyz.dtype == np.int16:
            xyz = _convert_from_lossy_integers(xyz)
    except:
        raise
    finally:
        F.close()

    trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
    return trajectory


def _topology_from_arrays(AtomID, AtomNames, ChainID, ResidueID, ResidueNames):
    topology = Topology()

    # assert that the ChainID is just an array of empty strings, which appears
    # to be the case in our test systems for this legacy format
    assert np.all(chainid == '' for chainid in ChainID), 'Im not prepaed to parse multiple chains'
    chain0 = topology.add_chain()


    # register the residues
    registered_residues = {}
    for i in np.argsort(ResidueID):
        residue_name = ResidueNames[i]
        if not isinstance(residue_name, str):
            residue_name = residue_name.decode()
        if ResidueID[i] not in registered_residues:
            res = topology.add_residue(residue_name, chain0)
            registered_residues[ResidueID[i]] = res

    # register the atoms
    for i in np.argsort(AtomID):
        atom_name = AtomNames[i]
        if not isinstance(atom_name, str):
            atom_name = atom_name.decode()
        element_symbol = atom_name.lstrip('0123456789')[0]
        element = mdtraj.pdb.element.get_by_symbol(element_symbol)
        print(atom_name)
        topology.add_atom(atom_name, element,
                         registered_residues[ResidueID[i]])

    topology.create_standard_bonds()
    return topology


# register this reader with mdtraj!
mdtraj.trajectory._LoaderRegistry['.lh5'] = load_legacy_hdf
