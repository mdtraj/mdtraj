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

import numpy as np
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif
from mdtraj.formats import mol2
from distutils.spawn import find_executable
from mdtraj.utils import enter_temp_directory
import tarfile
import pickle
import os
import numpy as np
import scipy.sparse
import subprocess
doc = DocStringFormatTester(mol2)


def test_load_mol2():
    trj = md.load(get_fn('imatinib.mol2'))
    ref_trj = md.load(get_fn('imatinib.pdb'))
    eq(trj.xyz, ref_trj.xyz)

    ref_top, ref_bonds = ref_trj.top.to_dataframe()
    top, bonds = trj.top.to_dataframe()
    eq(bonds, ref_bonds)

@skipif(find_executable('obabel') is None, 'You need obabel installed to run this test')
@skipif(os.environ.get("TRAVIS", None) == 'true', "Skip testing of entire FreeSolv database on Travis.")
def test_load_freesolv_gaffmol2_vs_sybylmol2_vs_obabelpdb():
    with enter_temp_directory():

        tar_filename = "freesolve_v0.3.tar.bz2"
        tar = tarfile.open(get_fn(tar_filename), mode="r:bz2")
        tar.extractall()
        tar.close()

        database = pickle.load(open("./v0.3/database.pickle"))

        for key in database:
            gaff_filename = "./v0.3/mol2files_gaff/%s.mol2" % key
            pdb_filename = "./v0.3/mol2files_gaff/%s.pdb" % key
            sybyl_filename = "./v0.3/mol2files_sybyl/%s.mol2" % key

            cmd = "obabel -imol2 %s -opdb > %s 2>/dev/null" % (sybyl_filename, pdb_filename)
            assert os.system(cmd) == 0

            t_pdb = md.load(pdb_filename)
            t_gaff = md.load(gaff_filename)
            t_sybyl = md.load(sybyl_filename)

            eq(t_pdb.n_atoms, t_gaff.n_atoms)
            eq(t_pdb.n_atoms, t_sybyl.n_atoms)

            eq(t_pdb.n_frames, t_gaff.n_frames)
            eq(t_pdb.n_frames, t_gaff.n_frames)

            eq(t_pdb.xyz, t_gaff.xyz, decimal=4)
            eq(t_pdb.xyz, t_sybyl.xyz, decimal=4)

            top_pdb, bonds_pdb = t_pdb.top.to_dataframe()
            top_gaff, bonds_gaff = t_gaff.top.to_dataframe()
            top_sybyl, bonds_sybyl = t_sybyl.top.to_dataframe()

            eq(top_sybyl.name.values, top_pdb.name.values)
            # eq(top_gaff.name.values, top_sybyl.name.values)  # THEY CAN HAVE DIFFERENT NAMES, so this isn't TRUE!

            def make_bonds_comparable(bond_array):
                """Create a bond connectivity matrix from a numpy array of atom pairs.  Avoids having to compare the order in which bonds are listed."""
                n_bonds = len(bond_array)
                data = np.ones(n_bonds)
                i = bond_array[:, 0]
                j = bond_array[:, 1]
                matrix = scipy.sparse.coo_matrix((data, (i, j)), shape=(t_pdb.n_atoms, t_pdb.n_atoms)).toarray()
                return matrix + matrix.T  # Symmetrize to account for (a ~ b) versus (b ~ a)

            bond_matrix_pdb = make_bonds_comparable(bonds_pdb)
            bond_matrix_gaff = make_bonds_comparable(bonds_gaff)
            bond_matrix_sybyl = make_bonds_comparable(bonds_sybyl)

            eq(bond_matrix_pdb, bond_matrix_gaff)
            eq(bond_matrix_pdb, bond_matrix_sybyl)

  # Third row from mol2 file copied below, used in testing.
  #       3 N1          8.5150   -0.1620    1.3310 n3        1 LIG     -0.732600

def test_mol2_dataframe():
    top, bonds = mol2.mol2_to_dataframes(get_fn("imatinib.mol2"))
    eq(top.name[2], "N1")
    eq(top.atype[2], "n3")
    eq(top.resName[2], "LIG")
    eq(float(top.charge[2]), -0.732600)

def test_mol2_warnings():
    trj = md.load_mol2(get_fn('lysozyme-ligand-tripos.mol2'))

def test_mol2_status_bits():
    trj = md.load_mol2(get_fn('status-bits.mol2'))
    eq(trj.topology.n_atoms, 18)
    eq(trj.topology.n_bonds, 18)
