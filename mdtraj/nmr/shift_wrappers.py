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


import os
import shutil
import subprocess

import numpy as np

from mdtraj.utils import enter_temp_directory, import_

##############################################################################
# Globals
##############################################################################

# Possible names for the external commands -- these are expected
# to be found in the PATH.
SHIFTX2 = ["shiftx2.py"]
SPARTA_PLUS = ["sparta+", "SPARTA+", "SPARTA+.linux"]
PPM = ["ppm_linux_64.exe"]

__all__ = [
    "chemical_shifts_shiftx2",
    "chemical_shifts_ppm",
    "chemical_shifts_spartaplus",
    "reindex_dataframe_by_atoms",
]


def find_executable(names):
    for possible in names:
        result = shutil.which(possible)
        if result is not None:
            return result
    return None


##############################################################################
# Code
##############################################################################


def compute_chemical_shifts(trj, model="shiftx2", **kwargs):
    """Predict chemical shifts of a trajectory using ShiftX2.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.
    model : str, optional, default="shiftx2"
        The program to use for calculating chemical shifts.  Must be one
        of shiftx2, ppm, or sparta+

    Returns
    -------
    results : pandas DataFrame
        Dataframe containing results, with index consisting of
        (resSeq, atom_name) pairs and columns for each frame in trj.

    Notes
    -----
    You must have the appropriate chemical soft programs installed
    and in your executable path.

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference, see docstrings for chemical_shifts_*
    for the various possible models.
    """
    if model == "shiftx2":
        return chemical_shifts_shiftx2(trj, **kwargs)
    elif model == "ppm":
        return chemical_shifts_ppm(trj, **kwargs)
    elif model == "sparta+":
        return chemical_shifts_spartaplus(trj, **kwargs)
    else:
        raise (ValueError("model must be one of shiftx2, ppm, or sparta+"))


def chemical_shifts_shiftx2(trj, pH=5.0, temperature=298.00):
    """Predict chemical shifts of a trajectory using ShiftX2.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.
    pH : float, optional, default=5.0
        pH value which gets passed to the ShiftX2 predictor.
    temperature : float, optional, default=298.00
        Temperature which gets passed to the ShiftX2 predictor.

    Returns
    -------
    results : pandas DataFrame
        Dataframe containing results, with index consisting of
        (resSeq, atom_name) pairs and columns for each frame in trj.

    Notes
    -----
    You must have ShiftX2 available on your path; see (http://www.shiftx2.ca/).

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Beomsoo Han, Yifeng Liu, Simon Ginzinger, and David Wishart.
       "SHIFTX2: significantly improved protein chemical shift
       prediction." J. Biomol. NMR, 50, 1 43-57 (2011)
    """
    pd = import_("pandas")
    binary = find_executable(SHIFTX2)
    if binary is None:
        raise OSError(
            "External command not found. Looked for {} in PATH. "
            "`chemical_shifts_shiftx2` requires the external program SHIFTX2, "
            "available at http://www.shiftx2.ca/".format(", ".join(SHIFTX2)),
        )

    results = []
    with enter_temp_directory():
        for i in range(trj.n_frames):
            fn = "./trj%d.pdb" % i
            trj[i].save(fn)
            subprocess.check_call(
                [
                    binary,
                    "-b",
                    fn,
                    "-p",
                    f"{pH:.1f}",
                    "-t",
                    f"{temperature:.2f}",
                ],
            )

            d = pd.read_csv("./trj%d.pdb.cs" % i)
            d.rename(
                columns={"NUM": "resSeq", "RES": "resName", "ATOMNAME": "name"},
                inplace=True,
            )
            d["frame"] = i
            results.append(d)

    return pd.concat(results).pivot_table(
        index=["resSeq", "name"],
        columns="frame",
        values="SHIFT",
    )


def chemical_shifts_ppm(trj):
    """Predict chemical shifts of a trajectory using ppm.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.

    Returns
    -------
    results : pandas.DataFrame
        Dataframe containing results, with index consisting of
        (resSeq, atom_name) pairs and columns for each frame in trj.

    Notes
    -----
    You must have ppm available on your path; see
    (http://spin.ccic.ohio-state.edu/index.php/download/index).

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Li, DW, and Bruschweiler, R. "PPM: a side-chain and backbone chemical
       shift predictor for the assessment of protein conformational ensembles."
       J Biomol NMR. 2012 Nov;54(3):257-65.
    """
    pd = import_("pandas")
    binary = find_executable(PPM)

    first_resSeq = trj.top.residue(0).resSeq

    if binary is None:
        raise OSError(
            f"External command not found. Looked for {PPM} in PATH. `chemical_shifts_ppm` "
            "requires the external program PPM, available at "
            "http://spin.ccic.ohio-state.edu/index.php/download/index",
        )

    with enter_temp_directory():
        trj.save("./trj.pdb")
        # -para old is on order to use newer ppm versions
        cmd = "%s -pdb trj.pdb -mode detail -para old" % binary

        return_flag = os.system(cmd)

        if return_flag != 0:
            raise (
                OSError(
                    f"Could not successfully execute command '{cmd}', check your PPM "
                    "installation or your input trajectory.",
                )
            )

        d = pd.read_table(
            "./bb_details.dat",
            index_col=False,
            header=None,
            sep=r"\s+",
        ).drop([3], axis=1)

        d = d.rename(columns={0: "resSeq", 1: "resName", 2: "name"})
        d["resSeq"] += first_resSeq - 1  # Fix bug in PPM that reindexes to 1
        d = d.drop("resName", axis=1)
        d = d.set_index(["resSeq", "name"])
        d.columns = np.arange(trj.n_frames)
        d.columns.name = "frame"

    return d


def _get_lines_to_skip(filename):
    """Determine the number of comment lines in a SPARTA+ output file."""
    format_string = """FORMAT %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"""
    handle = open(filename)
    for i, line in enumerate(handle):
        if line.find(format_string) != -1:
            return i + 2

    raise (Exception("No format string found in SPARTA+ file!"))


def chemical_shifts_spartaplus(trj, rename_HN=True):
    """Predict chemical shifts of a trajectory using SPARTA+.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.
    rename_HN : bool, optional, default=True
        SPARTA+ calls the amide proton "HN" instead of the standard "H".
        When True, this option renames the output as "H" to match the PDB
        and BMRB nomenclature.

    Returns
    -------
    results : pandas.DataFrame
        Dataframe containing results, with index consisting of
        (resSeq, atom_name) pairs and columns for each frame in trj.

    Notes
    -----
    You must have SPARTA+ available on your path; see
    (http://spin.niddk.nih.gov/bax/software/SPARTA+/). Also, the SPARTAP_DIR
    environment variable must be set so that SPARTA+ knows where to find
    its database files.

    Chemical shift prediction is for PROTEIN atoms; trajectory objects
    with ligands, solvent, ions, or other non-protein components may give
    UNKNOWN RESULTS.

    Please cite the appropriate reference below.

    References
    ----------
    .. [1] Shen, Y., and Bax, Ad. "SPARTA+: a modest improvement in empirical
       NMR chemical shift prediction by means of an artificial neural network."
       J. Biomol. NMR, 48, 13-22 (2010)
    """
    pd = import_("pandas")
    binary = find_executable(SPARTA_PLUS)
    if binary is None:
        raise OSError(
            f"External command not found. Looked for {SPARTA_PLUS} in PATH. "
            "`chemical_shifts_spartaplus` requires the external program SPARTA+, available at "
            "http://spin.niddk.nih.gov/bax/software/SPARTA+/",
        )

    names = [
        "resSeq",
        "resName",
        "name",
        "SS_SHIFT",
        "SHIFT",
        "RC_SHIFT",
        "HM_SHIFT",
        "EF_SHIFT",
        "SIGMA",
    ]

    with enter_temp_directory():
        for i in range(trj.n_frames):
            trj[i].save("./trj%d.pdb" % i)

        subprocess.check_call(
            [binary, "-in"] + [f"trj{i}.pdb" for i in range(trj.n_frames)] + ["-out", "trj0_pred.tab"],
        )

        lines_to_skip = _get_lines_to_skip("trj0_pred.tab")

        results = []
        for i in range(trj.n_frames):
            d = pd.read_table(
                "./trj%d_pred.tab" % i,
                names=names,
                header=None,
                sep=r"\s+",
                skiprows=lines_to_skip,
            )
            d["frame"] = i
            results.append(d)

    results = pd.concat(results)

    if rename_HN:
        results.name[results.name == "HN"] = "H"

    return results.pivot_table(
        index=["resSeq", "name"],
        columns="frame",
        values="SHIFT",
    )


def reindex_dataframe_by_atoms(trj, frame):
    """Reindex chemical shift output to use atom number (serial) indexing.

    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.
    frame : pandas.DataFrame
        Dataframe containing results, with index consisting of
        (resSeq, atom_name) pairs and columns for each frame in trj.

    Returns
    -------
    new_frame : pandas.DataFrame
        Dataframe containing results, with index consisting of atom
        indices (AKA the 'serial' entry in a PDB).  Columns correspond to
        each frame in trj.

    Notes
    -----
    Be aware that this function may DROP predictions if the atom naming
    is different between the input trajectory and the output of various
    chemical shift prediction tools.
    """

    top, bonds = trj.top.to_dataframe()
    top["serial"] = top.index
    top = top.set_index(["resSeq", "name"])

    new_frame = frame.copy()

    new_frame["serial"] = top.ix[new_frame.index].serial
    new_frame = new_frame.dropna().reset_index().set_index("serial").drop(["resSeq", "name"], axis=1)

    return new_frame
