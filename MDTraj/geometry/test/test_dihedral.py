import os
import sys
import shutil
from distutils.spawn import find_executable
import tempfile

import mdtraj as md
import numpy as np
from mdtraj.testing import get_fn, eq, SkipTest
from mdtraj.geometry.dihedral import (compute_chi1, compute_chi2, compute_chi3,
    compute_chi4, compute_omega, compute_phi, compute_psi)

traj = md.load(get_fn('2EQQ.pdb'))

def test_compute_phi():
    rx, x = compute_phi(traj, opt=True)
    ry, y = compute_phi(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_compute_psi():
    rx, x = compute_psi(traj, opt=True)
    ry, y = compute_psi(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_compute_chi():
    """Test that the proper chi angles are computed.

    The "reference data" is manually computed with pymol.
    """
    try:
        from jinja2 import Template
    except ImportError:
        raise SkipTest("jinja2 is required for this test")

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
            ["CA", "CB", "CG", "ND1"]]
CHI3_ATOMS = [["CB", "CG", "CD", "NE"],
            ["CB", "CG", "CD", "CE"],
            ["CB", "CG", "CD", "OE1"],
            ["CB", "CG", "SD", "CE"]]
CHI4_ATOMS = [["CG", "CD", "NE", "CZ"],
            ["CG", "CD", "CE", "NZ"]]

for n, chi_atoms in enumerate([CHI1_ATOMS, CHI2_ATOMS, CHI3_ATOMS, CHI4_ATOMS]):
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

    pymol = find_executable('pymol')
    if pymol is None:
        try_paths = ['/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL']
        for path in try_paths:
            if os.path.exists(path):
                pymol = path
                continue

    if pymol is None:
        raise SkipTest("pymol executable not found")

    try:
        tempdir = tempfile.mkdtemp()
        indices_fn = os.path.join(tempdir, 'indices')
        angles_fn = os.path.join(tempdir, 'angles')
        pymolscript_fn = os.path.join(tempdir, 'pymolscript.py')

        with open(pymolscript_fn, 'w') as f:
            f.write(pymol_script.render({'indices_fn': indices_fn,
                                         'angles_fn':  angles_fn}))

        os.system('%s %s -cr %s' % (pymol, get_fn('2EQQ.pdb'), pymolscript_fn))

        ref_indices = [np.loadtxt(indices_fn + '.%d' % i, dtype=int) for i in range(1,5)]
        ref_angles = [np.loadtxt(angles_fn + '.%d' % i) for i in range(1,5)]

    finally:
        shutil.rmtree(tempdir)

    indices, angles = zip(*[compute_chi1(traj, opt=False), compute_chi2(traj, opt=False),
                            compute_chi3(traj, opt=False), compute_chi4(traj, opt=False)])
    indices_opt, angles_opt = zip(*[compute_chi1(traj, opt=True), compute_chi2(traj, opt=True),
                                    compute_chi3(traj, opt=True), compute_chi4(traj, opt=True)])

    for x, y in zip(indices, ref_indices):
        eq(x, y)
    for x, y in zip(indices_opt, ref_indices):
        eq(x, y)
    for x, y in zip(angles, ref_angles):
        eq(x[0], y, decimal=4)
    for x, y in zip(angles_opt, ref_angles):
        eq(x[0], y, decimal=4)


def test_compute_omega():
    rx, x = compute_omega(traj, opt=True)
    ry, y = compute_omega(traj, opt=False)
    eq(rx, ry)
    eq(x, y)

def test_shape_when_none():
    t = md.load(get_fn('frame0.h5'))
    np.hstack((md.compute_phi(t)[1],
               md.compute_psi(t)[1],
               md.compute_chi1(t)[1],
               md.compute_chi2(t)[1],
               md.compute_chi3(t)[1],
               md.compute_chi1(t)[1],
               md.compute_omega(t)[1]))
