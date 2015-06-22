## MDTraj: a modern, open library for the analysis of molecular dynamics trajectories

[![Linux Build Status](https://travis-ci.org/mdtraj/mdtraj.svg?branch=master)](https://travis-ci.org/mdtraj/mdtraj)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/sqjgx3jh14vuxks5/branch/master?svg=true)](https://ci.appveyor.com/project/rmcgibbo/mdtraj/branch/master)
[![PyPI Version](https://badge.fury.io/py/mdtraj.svg)](https://pypi.python.org/pypi/mdtraj)
[![Binstar Badge](https://binstar.org/omnia/mdtraj/badges/version.svg)](https://binstar.org/omnia/mdtraj)
[![Downloads](https://img.shields.io/pypi/dm/mdtraj.svg)](https://pypi.python.org/pypi/mdtraj)

Read, write and analyze MD trajectories with only a few lines of Python code.

With MDTraj, you can

- Read and write from **every MD format imaginable** (`pdb`, `xtc`, `trr`, `dcd`, `binpos`, `netcdf`, `mdcrd`, `prmtop`, ...)
- Run **blazingly** fast RMSD calculations (4x the speed of the original Theobald QCP).
- Use tons of analysis functions like bonds/angles/dihedrals, hydrogen bonding identification, secondary structure assignment, NMR observables.
- Use a **lightweight API**, with a focus on **speed** and vectorized operations.

For details, see the website at [mdtraj.org](http://mdtraj.org). To get involved,
take a look at the [github issue tracker](https://github.com/mdtraj/mdtraj/issues)
and/or the user forums [discourse.mdtraj.org](http://discourse.mdtraj.org).

####  Citation [![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1101%2F008896-blue.svg)](http://doi.org/10.1101/008896)

MDTraj is research software. If you make use of MDTraj in scientific publications, please cite it. The BibTeX reference is
```
@article {
	author = {McGibbon, Robert T. and Beauchamp, Kyle A. and Schwantes, Christian R. and Wang, Lee-Ping and Hern{\'a}ndez, Carlos X. and Harrigan, Matthew P. and Lane, Thomas J. and Swails, Jason M. and Pande, Vijay S.},
	title = {MDTraj: a modern, open library for the analysis of molecular dynamics trajectories},
	year = {2014},
	doi = {10.1101/008896},
	publisher = {Cold Spring Harbor Labs Journals},
	journal = {bioRxiv}
}
```

#### License

GNU LGPL version 2.1, or at your option a later version of the license.
Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.
