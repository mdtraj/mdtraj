##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
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

"""
Code for merging and munging trajectories from FAH datasets.
"""
##############################################################################
# imports
##############################################################################

from __future__ import print_function, division
import os
import sys
import glob
import tarfile
from mdtraj.utils.contextmanagers import enter_temp_directory

__all__ = ['import_']

def keynat(string):
    '''A natural sort helper function for sort() and sorted()
    without using regular expression.

    >>> items = ('Z', 'a', '10', '1', '9')
    >>> sorted(items)
    ['1', '10', '9', 'Z', 'a']
    >>> sorted(items, key=keynat)
    ['1', '9', '10', 'Z', 'a']
    '''
    r = []
    for c in string:
        try:
            c = int(c)
            try:
                r[-1] = r[-1] * 10 + c
            except:
                r.append(c)
        except:
            r.append(c)
    return r

##############################################################################
# globals
##############################################################################


def concatenate_core17(path, top, output_filename):
    """Concatenate tar bzipped XTC files created by Folding@Home Core17.
    
    Parameters
    ----------
    path : str
        Path to directory containing "results-*.tar.bz2".  E.g. a single CLONE directory.
    top : mdtraj.Topology
        Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """
    
    from mdtraj.formats.hdf5 import HDF5TrajectoryFile
    import mdtraj as md
    
    glob_input = os.path.join(path, "results-*.tar.bz2")
    filenames = glob.glob(glob_input)
    filenames = sorted(filenames, key=keynat)
    
    if len(filenames) <= 0:
        return
    
    trj_file = HDF5TrajectoryFile(output_filename, mode='a')
    
    try:
        trj_file._create_earray(where='/', name='processed_filenames',atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass
    
    for filename in filenames:
        if filename in trj_file._handle.root.processed_filenames:
            print("Already processed %s" % filename)
            continue
        with enter_temp_directory():
            print("Processing %s" % filename)
            archive = tarfile.open(filename, mode='r:bz2')
            archive.extract("positions.xtc")
            trj = md.load("positions.xtc", top=top)
            
            for frame in trj:
                trj_file.write(coordinates=frame.xyz, cell_lengths=frame.unitcell_lengths, cell_angles=frame.unitcell_angles)
            
            trj_file._handle.root.processed_filenames.append([filename])
            

def concatenate_ocore(path, top, output_filename):
    """Concatenate XTC files created by Siegetank OCore.
    
    Parameters
    ----------
    path : str
        Path to stream directory containing frame directories /0, /1, /2
        etc.
    top : mdtraj.Topology
        Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """
    
    from mdtraj.formats.hdf5 import HDF5TrajectoryFile
    import mdtraj as md
    
    sorted_folders = sorted(os.listdir(path), key=lambda value: int(value))
    sorted_folders = [os.path.join(path, folder) for folder in sorted_folders]
    
    if len(sorted_folders) <= 0:
        return
    
    trj_file = HDF5TrajectoryFile(output_filename, mode='a')
    
    try:
        trj_file._create_earray(where='/', name='processed_folders',atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass
    
    for folder in sorted_folders:
        if folder in trj_file._handle.root.processed_folders:
            print("Already processed %s" % folder)
            continue
        print("Processing %s" % folder)
        xtc_filename = os.path.join(folder, "frames.xtc")
        trj = md.load(xtc_filename, top=top)
        
        for frame in trj:
            trj_file.write(coordinates=frame.xyz, cell_lengths=frame.unitcell_lengths, cell_angles=frame.unitcell_angles)
        
        trj_file._handle.root.processed_folders.append([folder])
            
