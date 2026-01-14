import json
import os
import shutil
import subprocess

import numpy as np

from mdtraj.utils import enter_temp_directory, import_

# Possible names for the external command -- expected to be found in PATH if not provided
PEPSI_SAXS = ["Pepsi-SAXS", "pepsi_saxs", "pepsi-saxs"]


def find_executable(names):
    """Return the first executable found on PATH from a list of names, or None."""
    which = shutil.which
    for possible in names:
        result = which(possible)
        if result is not None:
            return result
    return None


def _weighted_scale_and_offset(x, y, sigma):
    """Weighted linear least squares fit of y ≈ c*x + b.

    Parameters
    ----------
    x : np.ndarray, shape=(n,)
        Predictor values (e.g., mean predicted I(q)).
    y : np.ndarray, shape=(n,)
        Observed values (experimental I(q)).
    sigma : np.ndarray, shape=(n,)
        Standard deviation (errors) of y.

    Returns
    -------
    c : float
        Best-fit scale factor.
    b : float
        Best-fit constant background offset.
    """
    w = 1.0 / (sigma * sigma)
    # Build X with two columns: x and ones, solve (X^T W X) beta = X^T W y
    X0 = x
    X1 = np.ones_like(x)
    # XT_W = np.array([np.sum(X0 * w), np.sum(X1 * w)])[:, None]
    # Compute the elements of X^T W X
    S00 = np.sum(w * X0 * X0)
    S01 = np.sum(w * X0 * X1)
    S11 = np.sum(w * X1 * X1)
    # Compute the elements of X^T W y
    T0 = np.sum(w * X0 * y)
    T1 = np.sum(w * X1 * y)
    # Solve 2x2 linear system
    A = np.array([[S00, S01], [S01, S11]])
    b_vec = np.array([T0, T1])
    c, b = np.linalg.solve(A, b_vec)
    return float(c), float(b)


def _weighted_chi2(x, y, sigma, c, b):
    """Weighted reduced chi^2 for y vs c*x + b."""
    r = (y - (c * x + b)) / sigma
    n = r.size
    return float(np.sum(r * r) / max(1, n - 1))


def compute_saxs_pepsi(trj, saxs_data, pepsi_path=None, num_frames=100, output_path=None, bulk_e_density=0.334):
    """Compute SAXS predictions for a trajectory using Pepsi-SAXS (two-pass ensemble fit).

    This routine performs two passes based on [2] in order to avoid overfitting of solvation
    parameters:
    1) Fit solvation parameters per-frame allowing Pepsi-SAXS to minimize chi^2.
       Aggregate the best-fit r0 (normalized by the default r0) and d_rho values across frames.
    2) Refit all frames with ensemble-averaged solvation parameters and scaleFactor=1
       (i.e., no per-frame scaling). Average the predicted I(q) across frames, and then
       compute an optimal global scale and constant background that minimize chi^2
       against the experimental data.

    Parameters
    ----------
    trj : mdtraj.Trajectory
        Trajectory object to fit.
    saxs_data : str
        Path to the experimental SAXS data file for Pepsi-SAXS. Must include errors.
    pepsi_path : str or None, optional, default=None
        Path to the Pepsi-SAXS executable. If None, the function will search PATH
        for common executable names.
    num_frames : int, optional, default=100
        Number of frames to sample from the trajectory. Frames are evenly spaced.
        If `num_frames` exceeds the number of frames in `trj`, it is clamped.
    output_path : str or None, optional, default=None
        If provided, use this directory to store intermediate files. When None,
        a temporary working directory is used and cleaned automatically.

    Returns
    -------
    results : pandas.DataFrame
        Dataframe with columns:
        - q: scattering vector values
        - exp_I: experimental I(q)
        - exp_sigma: experimental errors
        - fit_I: ensemble-averaged prediction after global scale and background fit
        Indexed from 0..n_q-1 with a simple RangeIndex.
    chi2 : float
        Final reduced chi^2 value for the global fit

    Notes
    -----
    - This function calls the external program Pepsi-SAXS. You must have it installed
      and accessible either via `pepsi_path` or on your system PATH.
    - Files are written to either a temporary directory (default) or to `output_path`.
      The function does not print progress.

    References
    ----------
    .. [1] Grudinin S, Garkavenko M, Kazennov A. Pepsi-SAXS: an adaptive method for rapid and accurate
    computation of small-angle X-ray scattering profiles. Acta Crystallogr D Struct Biol. 2017
    May 1;73(Pt 5):449-464. doi: 10.1107/S2059798317005745. Epub 2017 Apr 27. PMID: 28471369.
    .. [2] Ahmed, M.C., Crehuet, R., Lindorff-Larsen, K. (2020). Computing, Analyzing, and Comparing
    the Radius of Gyration and Hydrodynamic Radius in Conformational Ensembles of Intrinsically Disordered
    Proteins. In: Kragelund, B.B., Skriver, K. (eds) Intrinsically Disordered Proteins.
    Methods in Molecular Biology, vol 2141. Humana, New York, NY.
    https://doi.org/10.1007/978-1-0716-0524-0_21
    """
    pd = import_("pandas")

    if not isinstance(saxs_data, str):
        raise TypeError("saxs_data must be a path string to the experimental data file.")
    saxs_data = os.path.abspath(saxs_data)
    if not os.path.exists(saxs_data):
        raise OSError(f"Experimental SAXS data not found at path: {saxs_data}")

    if pepsi_path is None:
        pepsi_path = find_executable(PEPSI_SAXS)
        if pepsi_path is None:
            raise OSError(
                "External command not found. Looked for {} in PATH. "
                "Please provide `pepsi_path` or install Pepsi-SAXS.".format(", ".join(PEPSI_SAXS)),
            )
    pepsi_path = os.path.abspath(pepsi_path)

    # Select evenly spaced frames
    n_total = trj.n_frames
    n_use = int(min(num_frames, n_total))
    if n_use <= 0:
        raise ValueError("num_frames must be >= 1 and the trajectory must contain frames.")
    # Evenly spaced indices from [0 .. n_total-1], inclusive
    idx = np.linspace(0, n_total - 1, num=n_use, dtype=int)

    def _run_pepsi_saxs_two_pass():
        # Save selected frames as PDBs
        for i, k in enumerate(idx):
            trj[k].save(f"./frame_{i}.pdb")

        # First pass: fit per-frame with --json and allow Pepsi to vary solvation params
        for i in range(n_use):
            fn = f"frame_{i}.pdb"
            out = f"frame_{i}.fit"
            subprocess.run(
                [pepsi_path, fn, saxs_data, "-o", out, "-cst", "--json"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )

        # Find default r0 using a single constrained fit on frame 0
        subprocess.run(
            [
                pepsi_path,
                "frame_0.pdb",
                saxs_data,
                "-o",
                "find_r0.fit",
                "--r0_min_factor",
                "1.0",
                "--r0_max_factor",
                "1.0",
                "--r0_N",
                "1",
                "-cst",
                "--json",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        try:
            with open("find_r0.json") as fh:
                data = json.load(fh)
            default_r0 = data["Fitting to the experimental curve"]["Best r0 found"]["value"]
        except Exception as e:
            raise OSError("Failed to read default r0 from Pepsi-SAXS JSON output.") from e

        # Aggregate best-fit per-frame parameters
        r0_factors = []
        d_rho_vals = []
        for i in range(n_use):
            json_fn = f"frame_{i}.json"
            try:
                with open(json_fn) as fh:
                    jd = json.load(fh)
                d_rho = jd["Fitting to the experimental curve"]["Best d_rho found"]["value"]
                r0_best = jd["Fitting to the experimental curve"]["Best r0 found"]["value"]
                r0_factors.append(r0_best / default_r0)
                d_rho_vals.append(d_rho)
            except Exception:
                # Skip missing or malformed files; keep lightweight behavior
                continue

        if len(r0_factors) == 0 or len(d_rho_vals) == 0:
            raise OSError("No valid Pepsi-SAXS JSON results found to compute ensemble parameters.")

        # Ensemble parameters: average r0 factor and d_rho
        r0_factor = float(np.mean(r0_factors))
        d_rho_percent = 100.0 * float(np.mean(d_rho_vals)) / bulk_e_density

        # Second pass: refit all frames using ensemble parameters and scaleFactor=1
        for i in range(n_use):
            fn = f"frame_{i}.pdb"
            out = f"Ave_AA_frame{i}.fit"
            subprocess.run(
                [
                    pepsi_path,
                    fn,
                    saxs_data,
                    "-o",
                    out,
                    "--dro",
                    f"{d_rho_percent}",
                    "--r0_min_factor",
                    f"{r0_factor}",
                    "--r0_max_factor",
                    f"{r0_factor}",
                    "--r0_N",
                    "1",
                    "--scaleFactor",
                    "1",
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )

        # Read and average I(q) across frames; read q, experimental I, and errors from first file
        # Use minimal parsing for speed; Pepsi-SAXS .fit files typically have 5-line headers.
        first = pd.read_csv("Ave_AA_frame0.fit", header=5, sep=r"\s+", usecols=[0, 1, 2, 3])
        q = first.iloc[:, 0].to_numpy()
        exp_I = first.iloc[:, 1].to_numpy()
        exp_sigma = first.iloc[:, 2].to_numpy()
        fit_first = first.iloc[:, 3].to_numpy()

        # Accumulate mean I(q) without storing all frames
        # Start with first frame, then add others
        sum_fit = fit_first.astype(float, copy=True)
        for i in range(1, n_use):
            d = pd.read_csv(f"Ave_AA_frame{i}.fit", header=5, sep=r"\s+", usecols=[3])
            sum_fit += d.iloc[:, 0].to_numpy()

        mean_fit = sum_fit / float(n_use)

        # Fast weighted linear least squares to find global scale and background
        c, b = _weighted_scale_and_offset(mean_fit, exp_I, exp_sigma)
        fit_I = c * mean_fit + b
        chi2 = _weighted_chi2(mean_fit, exp_I, exp_sigma, c, b)

        # Assemble results dataframe
        results = pd.DataFrame(
            {"q": q, "exp_I": exp_I, "exp_sigma": exp_sigma, "fit_I": fit_I},
        )
        return results, chi2

    if output_path is None:
        with enter_temp_directory():
            return _run_pepsi_saxs_two_pass()
    else:
        # Use user-provided directory (create if needed)
        output_path = os.path.abspath(output_path)
        os.makedirs(output_path, exist_ok=True)
        cwd = os.getcwd()
        try:
            os.chdir(output_path)
            return _run_pepsi_saxs_two_pass()
        finally:
            os.chdir(cwd)
