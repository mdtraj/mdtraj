##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
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
DCD, AMBER BINPOS, PDB, and HDF5.
"""
import sys
from types import ModuleType


all_by_module = {
    'mdtraj.core':              ['element'],
    'mdtraj._rmsd':             ['rmsd'],
    'mdtraj._lprmsd':           ['lprmsd'],
    'mdtraj.core.topology':     ['Topology'],
    'mdtraj.core.trajectory':   ['open', 'load', 'iterload', 'load_frame', 'Trajectory'],

    # Geometry
    'mdtraj.geometry.dihedral':
        ['compute_dihedrals', 'compute_phi', 'compute_psi', 'compute_omega',
        'compute_chi1', 'compute_chi2', 'compute_chi3','compute_chi4'],
    'mdtraj.geometry.distance':
        ['compute_distances', 'compute_displacements', 'compute_center_of_mass'],
    'mdtraj.geometry.angle':    ['compute_angles'],
    'mdtraj.geometry.contact':  ['compute_contacts'],
    'mdtraj.geometry.hbond':    ['wernet_nilsson', 'baker_hubbard', 'kabsch_sander'],
    'mdtraj.geometry.sasa':     ['shrake_rupley'],
    'mdtraj.geometry.drid':     ['compute_drid'],
    'mdtraj.geometry.dssp':     ['compute_dssp'],
    'mdtraj.geometry.rg':       ['compute_rg'],

    'mdtraj.nmr':
        ['compute_chemical_shifts', 'chemical_shifts_shiftx2', 'chemical_shifts_ppm',
         'chemical_shifts_spartaplus', 'reindex_dataframe_by_atoms',
         'compute_J3_HN_HA'],
    'mdtraj.reporters': ['HDF5Reporter', 'NetCDFReporter', 'DCDReporter'],

    # formats
    'mdtraj.formats.registry':  ['_FormatRegistry'],
    'mdtraj.formats.xtc':       ['load_xtc'],
    'mdtraj.formats.trr':       ['load_trr'],
    'mdtraj.formats.hdf5':      ['load_hdf5'],
    'mdtraj.formats.lh5':       ['load_lh5'],
    'mdtraj.formats.netcdf':    ['load_netcdf'],
    'mdtraj.formats.mdcrd':     ['load_mdcrd'],
    'mdtraj.formats.dcd':       ['load_dcd'],
    'mdtraj.formats.binpos':    ['load_binpos'],
    'mdtraj.formats.pdb':       ['load_pdb'],
    'mdtraj.formats.arc':       ['load_arc'],
    'mdtraj.formats.openmmxml': ['load_xml'],
    'mdtraj.formats.prmtop':    ['load_prmtop'],
    'mdtraj.formats.psf':       ['load_psf'],
    'mdtraj.formats.mol2':      ['load_mol2'],
    'mdtraj.formats.amberrst':  ['load_restrt', 'load_ncrestrt'],
}

# modules that should be imported when accessed as attributes of mdtraj
attribute_modules = frozenset(['version'])

def test(label='full', verbose=2):
    """Run tests for mdtraj using nose.

    Parameters
    ----------
    label : {'fast', 'full'}
        Identifies the tests to run. The fast tests take about 10 seconds,
        and the full test suite takes about two minutes (as of this writing).
    verbose : int, optional
        Verbosity value for test outputs, in the range 1-10. Default is 2.
    """
    import mdtraj
    from mdtraj.testing.nosetester import MDTrajTester
    tester = MDTrajTester(mdtraj)
    return tester.test(label=label, verbose=verbose, extra_argv=('--exe',))
# prevent nose from discovering this function, or otherwise when its run
# the test suite in an infinite loop
test.__test__ = False


object_origins = {}
for module, items in all_by_module.items():
    for item in items:
        object_origins[item] = module

class module(ModuleType):
    """Automatically import objects from the modules."""

    def __getattr__(self, name):
        if name in object_origins:
            module = __import__(object_origins[name], None, None, [name])
            for extra_name in all_by_module[module.__name__]:
                setattr(self, extra_name, getattr(module, extra_name))
            return getattr(module, name)
        elif name in attribute_modules:
            __import__('mdtraj.' + name)
        return ModuleType.__getattribute__(self, name)

    def __dir__(self):
        """Just show what we want to show."""
        result = list(new_module.__all__)
        result.extend(('__file__', '__path__', '__doc__', '__all__',
                       '__docformat__', '__name__', '__path__',
                       '__package__', '__version__'))
        return result

# keep a reference to this module so that it's not garbage collected
old_module = sys.modules[__name__]

# setup the new module and patch it into the dict of loaded modules
new_module = sys.modules[__name__] = module(__name__)
new_module.__dict__.update({
    '__file__':         __file__,
    '__package__':      'mdtraj',
    '__path__':         __path__,
    '__doc__':          __doc__,
    '__all__':          tuple(object_origins) + tuple(attribute_modules),
    '__docformat__':    'restructuredtext en',
    'test':             test,
})
