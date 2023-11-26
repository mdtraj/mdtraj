## Note on current project status

TLDR: MDTraj is currently undergoing a transition to new maintainers (the Folding@home consortium).
Please bear with us as we slowly respond to issues and clean up.

Hello!

For anybody wondering, the MDTraj repo is currently undergoing a transition to new maintainers;
The Folding@home consortium (@FoldingAtHome) will now be taking over support, maintainence
and management of MDTraj with @xuhuihuang leading repo management, and @mattwthompson, @apayne97,
and myself (@sukritsingh) helping to manage and maintain this repo more directly.

Existing supporters, contributors, etc. are all welcome to contribute as much as they wish.
The intention behind this transition is to ensure MDTraj receives continued support.

As we are slowly transitioning to this new support, please bear with us as we slowly respond
to issues, clean up any PRs, etc.

I have opened a [Discussions page](https://github.com/mdtraj/mdtraj/discussions) on github
for folks asking for help with code issues/not getting things working. Our hope is that the
issues page will be for discrete bugs, feature requests, related discussions etc., but this
is very flexible!

Best,

Sukrit Singh (@sukritsingh)

July 6th, 2023

--------------------------------

## MDTraj: an open-source library for analysis of molecular dynamics trajectories

[![Build Status](https://github.com/mdtraj/mdtraj/actions/workflows/main.yaml/badge.svg)](https://github.com/mdtraj/mdtraj/actions)
[![PyPI Version](https://badge.fury.io/py/mdtraj.svg)](https://pypi.python.org/pypi/mdtraj)
[![Anaconda-Server Version](https://anaconda.org/conda-forge/mdtraj/badges/version.svg)](https://anaconda.org/conda-forge/mdtraj)
[![Anaconda-Server Downloads](https://anaconda.org/conda-forge/mdtraj/badges/downloads.svg)](https://anaconda.org/conda-forge/mdtraj)
[![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1016%2Fj.bpj.2015.08.015-blue.svg)](http://doi.org/10.1016/j.bpj.2015.08.015)

Read, write and analyze MD trajectories with only a few lines of Python code.

With MDTraj, you can

- Read and write from **every MD format imaginable** (`pdb`, `xtc`, `trr`, `dcd`, `binpos`, `netcdf`, `mdcrd`, `prmtop`, `gsd`, ...)
- Run **blazingly** fast RMSD calculations (4x the speed of the original Theobald QCP).
- Use tons of analysis functions like bonds/angles/dihedrals, hydrogen bonding identification, secondary structure assignment, NMR observables.
- Use a **lightweight API**, with a focus on **speed** and vectorized operations.

For details, see the website at [mdtraj.org](http://mdtraj.org). To get involved,
take a look at the [github issue tracker](https://github.com/mdtraj/mdtraj/issues)
and/or the [gitter room](https://gitter.im/mdtraj/community).

####  Citation

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
