import itertools
import os
from distutils.spawn import find_executable

import mdtraj as md
import numpy as np
import pytest
from mdtraj.geometry.dihedral import (compute_chi1, compute_chi2, compute_chi3,
                                      compute_chi4, compute_omega, compute_phi,
                                      compute_psi, compute_chi5)
from mdtraj.testing import eq


@pytest.fixture()
def traj(get_fn):
    return md.load(get_fn('2EQQ.pdb'))


@pytest.fixture()
def pymol():
    pymol = find_executable('pymol')
    if pymol is None:
        try_paths = ['/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL']
        for path in try_paths:
            if os.path.exists(path):
                pymol = path
                continue

    if pymol is None:
        raise pytest.skip("pymol executable not found")
    return pymol


def test_compute_phi(traj):
    rx, x = compute_phi(traj, opt=True)
    ry, y = compute_phi(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_compute_psi(traj):
    rx, x = compute_psi(traj, opt=True)
    ry, y = compute_psi(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_compute_chi(get_fn, tmpdir, pymol, traj):
    """Test that the proper chi angles are computed.

    The "reference data" is manually computed with pymol.
    """
    try:
        from jinja2 import Template
    except ImportError:
        pytest.skip("jinja2 is required for this test")

    pymol_script = Template('''
import numpy as np
from pymol import cmd, CmdException

CHI1_ATOMS = [["N", "CA", "CB", "CG"],
              ["N", "CA", "CB", "CG1"],
              ["N", "CA", "CB", "SG"],
              ["N", "CA", "CB", "OG"],
              ["N", "CA", "CB", "OG1"]]
CHI2_ATOMS = [["CA", "CB", "CG", "CD"],
            ["CA", "CB", "CG", "CD1"],
            ["CA", "CB", "CG1", "CD1"],
            ["CA", "CB", "CG", "OD1"],
            ["CA", "CB", "CG", "ND1"],
            ["CA", "CB", "CG", "SD"]]
CHI3_ATOMS = [["CB", "CG", "CD", "NE"],
            ["CB", "CG", "CD", "CE"],
            ["CB", "CG", "CD", "OE1"],
            ["CB", "CG", "SD", "CE"]]
CHI4_ATOMS = [["CG", "CD", "NE", "CZ"],
            ["CG", "CD", "CE", "NZ"]]

CHI5_ATOMS = ["CD", "NE", "CZ", "NH1"]

for n, chi_atoms in enumerate([CHI1_ATOMS, CHI2_ATOMS, CHI3_ATOMS, CHI4_ATOMS, CHI5_ATOMS]):
    indices, angles = [], []
    for r in range(1, cmd.count_atoms('name ca')+1):
        r_indices, r_angles = None, None
        for atoms in chi_atoms:
            try:
                selection = ['%d/%s' % (r, e) for e in atoms]
                r_indices = [cmd.id_atom(s) - 1 for s in selection]
                r_angles =  cmd.get_dihedral(*selection) * np.pi / 180.0
            except CmdException:
                pass
        if r_indices is not None:
            indices.append(r_indices)
            angles.append(r_angles)
    np.savetxt("{{indices_fn}}.%d" % (n+1), indices, '%d')
    np.savetxt("{{angles_fn}}.%d" % (n+1), angles)

cmd.quit()''')

    indices_fn = os.path.join(tmpdir, 'indices')
    angles_fn = os.path.join(tmpdir, 'angles')
    pymolscript_fn = os.path.join(tmpdir, 'pymolscript.py')

    with open(pymolscript_fn, 'w') as f:
        f.write(pymol_script.render({'indices_fn': indices_fn,
                                     'angles_fn': angles_fn}))

    os.system('%s %s -cr %s' % (pymol, get_fn('2EQQ.pdb'), pymolscript_fn))

    ref_indices = [np.loadtxt(indices_fn + '.%d' % i, dtype=int) for i in range(1, 5)]
    ref_angles = [np.loadtxt(angles_fn + '.%d' % i) for i in range(1, 5)]

    indices, angles = zip(*[compute_chi1(traj, opt=False), compute_chi2(traj, opt=False),
                            compute_chi3(traj, opt=False), compute_chi4(traj, opt=False),
                            compute_chi5(traj, opt=False)])
    indices_opt, angles_opt = zip(*[compute_chi1(traj, opt=True), compute_chi2(traj, opt=True),
                                    compute_chi3(traj, opt=True), compute_chi4(traj, opt=True),
                                    compute_chi5(traj, opt=True)])

    for x, y in zip(indices, ref_indices):
        eq(x, y)
    for x, y in zip(indices_opt, ref_indices):
        eq(x, y)
    for x, y in zip(angles, ref_angles):
        eq(x[0], y, decimal=4)
    for x, y in zip(angles_opt, ref_angles):
        eq(x[0], y, decimal=4)


def test_compute_omega(traj):
    rx, x = compute_omega(traj, opt=True)
    ry, y = compute_omega(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_shape_when_none(get_fn):
    t = md.load(get_fn('frame0.h5'))
    np.hstack((md.compute_phi(t)[1],
               md.compute_psi(t)[1],
               md.compute_chi1(t)[1],
               md.compute_chi2(t)[1],
               md.compute_chi3(t)[1],
               md.compute_chi1(t)[1],
               md.compute_omega(t)[1]))


def test_dihedral_1(pymol, tmpdir):
    xyz = '''MODEL        0
ATOM      1    A ACE     1       4.300  13.100   8.600  1.00  0.00
ATOM      2    B ACE     1       5.200  13.600   8.800  1.00  0.00
ATOM      3    C ACE     1       4.900  14.300   9.600  1.00  0.00
ATOM      4    D ACE     1       5.600  14.200   7.900  1.00  0.00
    '''
    script = '''
from pymol import cmd
with open('output.txt', 'w') as f:
    f.write('%f' % cmd.get_dihedral('1/A', '1/B', '1/C', '1/D'))
'''
    prevdir = os.path.abspath('.')
    try:
        os.chdir(tmpdir)
        with open('xyz.pdb', 'w') as f:
            f.write(xyz)
        with open('pymolscript.py', 'w') as f:
            f.write(script)

        os.system('%s %s -cr %s' % (pymol, 'xyz.pdb', 'pymolscript.py'))
        with open('output.txt') as f:
            pymol_value = np.deg2rad(float(f.read()))
        t = md.load('xyz.pdb')
    finally:
        os.chdir(prevdir)

    mdtraj_value = md.compute_dihedrals(t, [[0, 1, 2, 3]])[0, 0]

    np.testing.assert_array_almost_equal(pymol_value, mdtraj_value)


def test_dihedral_2chains(get_fn):
    # make sure that comput_phi is finding dihedrals from all of the chains
    # in a multi-chain topology
    t = md.load_pdb(get_fn('4OH9.pdb'))

    # remove the water
    water_indices = [a.index for a in t.top.atoms if a.residue.name != 'HOH']
    t.restrict_atoms(water_indices)

    # okay we've got two protein chains
    assert t.top.n_chains == 2
    phi_indices, angles = md.compute_phi(t)

    for chain in t.top.chains:
        chain_indices = [a.index for a in chain.atoms]

        # assert that at least one of the phi_indices involves atoms in this chain
        assert any(i in chain_indices for i in np.concatenate(phi_indices))


def test_generator():
    N_FRAMES = 2
    N_ATOMS = 5
    xyz = np.asarray(np.random.randn(N_FRAMES, N_ATOMS, 3), dtype=np.float32)
    ptraj = md.Trajectory(xyz=xyz, topology=None)

    quartets = np.array(list(itertools.combinations(range(N_ATOMS), 4)), dtype=np.int32)
    quartets2 = itertools.combinations(range(N_ATOMS), 4)
    a = md.compute_dihedrals(ptraj, quartets)
    b = md.compute_dihedrals(ptraj, quartets2)
    eq(a, b)
