RMSD : Fast SSE2/SSE3-based code for RMSD computation
=====================================================

This directory is not a python package (no `__init__.py`). It contains the
source code that gets compiled to produce the IRMSD implementation code, a
single python shared object file importable at `mdtraj._rmsd`. The object
oriented python wrapper for this code is  `mdtraj.rmsd`, and is written in
the `rmsd.py` file.

