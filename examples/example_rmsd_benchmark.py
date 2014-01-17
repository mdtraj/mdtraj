"""Benchmarking RMSD speed
"""

import numpy as np
import mdtraj as md

# To benchmark the speed of the RMSD calculation, it's not really
# necessary to use 'real' coordinates, so let's just generate
# some random numbers from a normal distribution for the cartesian
# coordinates.

t = md.Trajectory(xyz=np.random.randn(1000, 100, 3), topology=None)
print t

# The Theobald QCP method requires centering the invidual conformations
# the origin. That's done on the fly when we call :func:`md.rmsd()`.

import time; start=time.time()
for i in range(100):
    md.rmsd(t, t, 0)
print 'md.rmsd(): %.2f rmsds / s' % ((t.n_frames * 100) / (time.time() - start))

# But for some applications like clustering, we want to run many
# rmsd() calculations per trajectory frame. Under these circumstances,
# the centering of the trajectories is going to be done many times, which
# leads to a slight slowdown. If we manually center the trajectory
# and then inform the rmsd() function that the centering has been
# precentered, we can achieve ~2x speedup, depending on your machine
# and the number of atoms.

t.center_coordinates()
start = time.time()
for i in range(100):
    md.rmsd(t, t, 0, precentered=True)
print 'md.rmsd(precentered=True): %.2f rmsds / s' % ((t.n_frames * 100) / (time.time() - start))

# Just for fun, let's compare this code to the straightforward
# numpy implementation of the same algorithm, which mdtraj has
# (mostly for testing) in the  mdtraj.geometry.alignment subpackage

from mdtraj.geometry.alignment import rmsd_qcp
start = time.time()
for k in range(t.n_frames):
    rmsd_qcp(t.xyz[0], t.xyz[k])
print 'pure numpy rmsd_qcp(): %.2f rmsds / s' % (t.n_frames / (time.time() - start))

# The :func:`md.rmsd()` code is *a lot* faster. If you go look at the :func:`rmsd_qcp`
# source code, you'll see that it's not because that code is particularly slow or
# unoptimized. It's about as good as you can do with numpy. The reason for the speed
# difference is that an inordinate amount of time was put into hand-optimizing
# an SSE3 implementation in C for the :func:`md.rmsd()` code.
