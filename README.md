## MDTraj: a modern, open library for the analysis of molecular dynamics trajectories

[![Linux Build Status](https://travis-ci.org/mdtraj/mdtraj.svg?branch=master)](https://travis-ci.org/mdtraj/mdtraj)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/sqjgx3jh14vuxks5/branch/master?svg=true)](https://ci.appveyor.com/project/rmcgibbo/mdtraj/branch/master)
[![PyPI Version](https://badge.fury.io/py/mdtraj.svg)](https://pypi.python.org/pypi/mdtraj)
[![Anaconda-Server Version](https://anaconda.org/omnia/mdtraj/badges/version.svg)](https://anaconda.org/omnia/mdtraj)
[![Anaconda-Server Downloads](https://anaconda.org/omnia/mdtraj/badges/downloads.svg)](https://anaconda.org/omnia/mdtraj)
[![Research software impact](http://depsy.org/api/package/pypi/mdtraj/badge.svg)](http://depsy.org/package/python/mdtraj)

Read, write and analyze MD trajectories with only a few lines of Python code.

With MDTraj, you can

- Read and write from **every MD format imaginable** (`pdb`, `xtc`, `trr`, `dcd`, `binpos`, `netcdf`, `mdcrd`, `prmtop`, ...)
- Run **blazingly** fast RMSD calculations (4x the speed of the original Theobald QCP).
- Use tons of analysis functions like bonds/angles/dihedrals, hydrogen bonding identification, secondary structure assignment, NMR observables.
- Use a **lightweight API**, with a focus on **speed** and vectorized operations.

For details, see the website at [mdtraj.org](http://mdtraj.org). To get involved,
take a look at the [github issue tracker](https://github.com/mdtraj/mdtraj/issues)
and/or the user forums [discourse.mdtraj.org](http://discourse.mdtraj.org).

####  Citation [![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1016%2Fj.bpj.2015.08.015-blue.svg)](http://doi.org/10.1016/j.bpj.2015.08.015)

MDTraj is research software. If you make use of MDTraj in scientific publications, please cite it. The BibTeX reference is
```
@article{McGibbon2015MDTraj,
    author = {McGibbon, Robert T. and Beauchamp, Kyle A. and Harrigan, Matthew P. and Klein, Christoph and Swails, Jason M. and Hern{\'a}ndez, Carlos X.  and Schwantes, Christian R. and Wang, Lee-Ping and Lane, Thomas J. and Pande, Vijay S.},
    title = {MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories},
	journal = {Biophysical Journal},
    volume = {109},
    number = {8},
    pages = {1528 -- 1532},
    year = {2015},
	doi = {10.1016/j.bpj.2015.08.015}
}
```

#### License

GNU LGPL version 2.1, or at your option a later version of the license.
Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.
