##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors: Robert McGibbon
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
#############################################################################


##############################################################################
# Imports
##############################################################################

import os
import string
from distutils.spawn import find_executable as _find_executable
import pandas as pd

from mdtraj.utils import enter_temp_directory

##############################################################################
# Globals
##############################################################################

# Possible names for the external commands -- these are expected
# to be found in the PATH.
SHIFTX2 = ['shiftx2.py']
SPARTA_PLUS = ['sparta+', 'SPARTA+', 'SPARTA+.linux']
PPM = ['ppm_linux_64.exe']

__all__ = ['chemical_shifts_shiftx2', 'chemical_shifts_ppm', 'chemical_shifts_spartaplus']


def find_executable(names):
    for possible in names:
        result = _find_executable(possible)
        if result is not None:
            return result
    return None


##############################################################################
# Code
##############################################################################


def chemical_shifts_shiftx2(trj):
    """Predict chemical shifts of a trajectory using ShiftX2.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.

    Returns
    -------
    results : pandas DataFrame
        Dataframe containing results, with index consisting of (resSeq, atom_name) pairs
        and columns for each frame in trj.

    Notes
    -----
    You must have ShiftX2 available on your path; see (http://www.shiftx2.ca/).

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Beomsoo Han, Yifeng Liu, Simon Ginzinger, and David Wishart. "SHIFTX2: significantly improved protein chemical shift prediction." J. Biomol. NMR, 50, 1 43-57 (2011)
    """

    binary = find_executable(SHIFTX2)
    if binary is None:
        raise OSError('External command not found. Looked for %s in PATH. `chemical_shifts_shiftx2` requires the external program SHIFTX2, available at http://www.shiftx2.ca/' % ', '.join(SHIFTX2))

    results = []
    with enter_temp_directory():
        for i in range(trj.n_frames):
            trj[i].save("./trj%d.pdb" % i)
        cmd = "%s -b 'trj*.pdb'" % binary

        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your ShiftX2 installation or your input trajectory." % cmd))

        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d.pdb.cs" % i)
            d.rename(columns={"NUM": "resSeq", "RES": "resName", "ATOMNAME": "name"}, inplace=True)
            d["frame"] = trj.time[i]
            results.append(d)

    results = pd.concat(results)
    results = results.pivot_table(rows=["resSeq", "name"], cols="frame", values="SHIFT")
    return results


def chemical_shifts_ppm(trj):
    """Predict chemical shifts of a trajectory using ppm.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.

    Returns
    -------
    results : pandas.DataFrame
        Dataframe containing results, with index consisting of (resSeq, atom_name) pairs
        and columns for each frame in trj.

    Notes
    -----
    You must have ppm available on your path; see (http://spinportal.magnet.fsu.edu/ppm/ppm.html).

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Li, DW, and Bruschweiler, R. "PPM: a side-chain and backbone chemical shift predictor for the assessment of protein conformational ensembles." J Biomol NMR. 2012 Nov;54(3):257-65.
    """
    first_resSeq = trj.top.residue(0).resSeq

    binary = find_executable(PPM)
    if binary is None:
        raise OSError('External command not found. Looked for %s in PATH. `chemical_shifts_ppm` requires the external program PPM, available at http://spinportal.magnet.fsu.edu/ppm/ppm.html' % ', '.join(PPM))

    with enter_temp_directory():
        trj.save("./trj.pdb")
        cmd = "%s -pdb trj.pdb -mode detail" % binary

        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your PPM installation or your input trajectory." % cmd))

        d = pd.read_csv("./bb_details.dat", delim_whitespace=True)
        columns = ["resSeq", "resName", "name", "expt", "other"]

        d = pd.read_csv("./bb_details.dat", delim_whitespace=True, header=None).drop([0, 4], axis=1)
        d = d.rename(columns={1: "resSeq", 2: "resName", 3: "name"})
        d["resSeq"] += first_resSeq - 1  # Fix bug in PPM that reindexes to 1
        d = d.drop("resName", axis=1)
        d = d.set_index(["resSeq", "name"])
        d.columns = trj.time
        d.columns.name = "frame"

    return d


def _get_lines_to_skip(filename):
    """Determine the number of comment lines in a SPARTA+ output file."""
    format_string = """FORMAT %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"""
    handle = open(filename)
    for i, line in enumerate(handle):
        if line.find(format_string) != -1:
            return i + 2

    raise(Exception("No format string found in SPARTA+ file!"))


def chemical_shifts_spartaplus(trj):
    """Predict chemical shifts of a trajectory using SPARTA+.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.

    Returns
    -------
    results : pandas.DataFrame
        Dataframe containing results, with index consisting of (resSeq, atom_name) pairs
        and columns for each frame in trj.

    Notes
    -----
    You must have ppm available on your path; see (http://spin.niddk.nih.gov/bax/software/SPARTA+/).
    Also, the SPARTAP_DIR environment variable must be set so that SPARTA+
    knows where to find its database files.

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Shen, Y., and Bax, Ad. "SPARTA+: a modest improvement in empirical NMR chemical shift prediction by means of an artificial neural network." J. Biomol. NMR, 48, 13-22 (2010)
    """

    binary = find_executable(SPARTA_PLUS)
    if binary is None:
        raise OSError('External command not found. Looked for %s in PATH. `chemical_shifts_spartaplus` requires the external program SPARTA+, available at http://spin.niddk.nih.gov/bax/software/SPARTA+/' % ', '.join(SPARTA_PLUS))

    names = ["VARS", "resSeq", "resName", "name", "SS_SHIFT", "SHIFT", "RC_SHIFT", "HM_SHIFT", "EF_SHIFT", "SIGMA"]

    with enter_temp_directory():
        for i in range(trj.n_frames):
            trj[i].save("./trj%d.pdb" % i)

        cmd = "%s -in %s" % (binary, string.join(["trj%d.pdb" % i for i in range(trj.n_frames)]))

        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your SPARTA+ installation or your input trajectory." % cmd))

        lines_to_skip = _get_lines_to_skip("trj0_pred.tab")

        results = []
        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d_pred.tab" % i, skiprows=lines_to_skip, delim_whitespace=True, header=None, names=names)
            d["frame"] = trj.time[i]
            results.append(d)

    results = pd.concat(results)
    results = results.pivot_table(rows=["resSeq", "name"], cols="frame", values="SHIFT")

    return results
