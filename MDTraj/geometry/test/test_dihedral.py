import mdtraj as md
from mdtraj.geometry.dihedral import (compute_chi_all, compute_omega,
                                      compute_phi, compute_psi)
from mdtraj.testing import get_fn, eq
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


def test_compute_chi_all():
    rx, x = compute_chi(traj, opt=True)
    ry, y = compute_chi(traj, opt=False)
    eq(rx, ry)
    eq(x, y)


def test_compute_omega():
    rx, x = compute_omega(traj, opt=True)
    ry, y = compute_omega(traj, opt=False)
    eq(rx, ry)
    eq(x, y)
