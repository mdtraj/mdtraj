
# MDTraj: an open-source library for analysis of molecular dynamics trajectories

[![Build Status](https://github.com/mdtraj/mdtraj/actions/workflows/main.yaml/badge.svg)](https://github.com/mdtraj/mdtraj/actions)
[![PyPI Version](https://badge.fury.io/py/mdtraj.svg)](https://pypi.python.org/pypi/mdtraj)
[![Anaconda-Server Version](https://anaconda.org/conda-forge/mdtraj/badges/version.svg)](https://anaconda.org/conda-forge/mdtraj)
[![Anaconda-Server Downloads](https://anaconda.org/conda-forge/mdtraj/badges/downloads.svg)](https://anaconda.org/conda-forge/mdtraj)


[![SPEC 0 — Minimum Supported Dependencies](https://img.shields.io/badge/SPEC-0-green?labelColor=%23004811&color=%235CA038)](https://scientific-python.org/specs/spec-0000/)
[![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1016%2Fj.bpj.2015.08.015-blue.svg)](http://doi.org/10.1016/j.bpj.2015.08.015)
[![Documentation Status](https://readthedocs.org/projects/mdtraj/badge/?version=latest)](https://mdtraj.readthedocs.io/en/latest/?badge=latest)

Read, write and analyze MD trajectories with only a few lines of Python code.

With MDTraj, you can

- Read and write from **every MD format imaginable** (`pdb`, `xtc`, `trr`, `dcd`, `netcdf`, `mdcrd`, `prmtop`, `gsd`, ...)
- Run **blazingly** fast RMSD calculations (4x the speed of the original Theobald QCP).
- Use tons of analysis functions like bonds/angles/dihedrals, hydrogen bonding identification, secondary structure assignment, NMR observables.
- Use a **lightweight API**, with a focus on **speed** and vectorized operations.

For details, see the [documentation](https://mdtraj.readthedocs.io/en/latest/). To get involved,
take a look at the [GitHub issue tracker](https://github.com/mdtraj/mdtraj/issues)
and/or the [Discussions](https://github.com/mdtraj/mdtraj/discussions) section.

##  Citation

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

## License

GNU LGPL version 2.1, or at your option a later version of the license.
Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

## Project Status
MDTraj is actively maintained by community members and representatives from Folding@Home, The Molecular Sciences Software Institute, and Open Force Field Initiative. We welcome new users and contributors!

If you are interested in contributing, please check out the [repository issue tracker](https://github.com/mdtraj/mdtraj/issues). For feature ideas or general questions, please post in the [repository discussions](https://github.com/mdtraj/mdtraj/discussions).
