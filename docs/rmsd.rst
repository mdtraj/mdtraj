RMSD : Fast SSE2/SSE3-based code for RMSD computation
=====================================================

MDTraj includes the `IRMSD <https://github.com/simtk/irmsd>`_ library for computing the optimal root-mean-square-deviation between pairs of protein conformations. It is based on the Theobald QCP method, and because of an efficient matrix multiply routine, is several-fold faster than other RMSD packages using off-the-shelf naive or generic high-performance matrix multiplies. In particular:


* With one thread, IRMSD is 4x the speed of the original Theobald code and over 3x as fast as Theobald code modified to use GotoBLAS.
* IRMSD features automatic parallelization over multiple RMSD computations; with four threads, it is twice the speed of similarly parallelized Theobald code using GotoBLAS.
* IRMSD reaches the machine theoretical peak arithmetic performance when using multiple threads and small structures that fit into L1 cache.
* IRMSD reaches the machine limit of memory bandwidth when dealing with very large structures that do not fit into cache.

IRMSD also fixes a small numerical instability in the Theobald QCP method that manifests as a lack of convergence in approximately one in 10^9 RMSD computations.

The IRMSD code was written by Imran S. Haque. **If you use this code in computations that result in publication, please cite our paper: [upcoming]**

Example Usage
-------------

:: 

  >>> import mdtraj as md
  >>> import mdtraj.testing
  # load up a little trajectory that's distributed with mdtraj, just for this example
  >>> t = md.load(mdtraj.testing.get_fn('traj.h5'))
  >>> print t
  <mdtraj.Trajectory with 100 frames, 22 atoms>

It's alanine dipeptide. if you were doing this in your own code, you'd probably want to discard the coordinates for the hydrogens. In general, RMSD doesn't take into account the exchange symmetry of identical atoms, for example the hydrogens on a methyl group. ::
  
  >>> # prepare the trajectory for RMSD computations. this aligns the
  >>> # data to the appropriate btye-boundaries in memory, and adds extra     
  >>> # "padding atoms" at (0, 0, 0) to facilitate the use of fast SSE2/3 
  >>> # instructions.
  >>> c = mdtraj.rmsd_cache(t)
  
``rmsds_to`` computes the RMSD from each conformation in the cache to
one of the conformations in a different cache. Here, we compute the
distance from the first conformation to all of the others. ::

  >>> values = c.rmsds_to(c, 0)
  >>> print values.shape
  (100,)
  >>> print values[0:10]
  array([ 0.0001041 ,  0.00474321,  0.00918401,  0.01253726,  0.01449224,
          0.01535706,  0.01565618,  0.01561841,  0.01534577,  0.01479676], dtype=float32)


You might notice that the distance between the zero-th frame and itself is not exactly zero, its 0.0001041 nm. This is to be expected. The numerical precision is about one one-hundredth of an angstrom.

API
---

.. currentmodule :: mdtraj.rmsd

.. autofunction :: mdtraj.rmsd_cache

.. autoclass :: mdtraj.rmsd.RMSDCache

     

Low Level C-Code
----------------

The compiled code is written in C, and wrapped for python accessibility in cython. Both the C and Cython codes are compiled together into a single shared object file, accessible at ``mdtraj._rmsd``. Two public functions are declared.

.. automodule:: mdtraj._rmsd