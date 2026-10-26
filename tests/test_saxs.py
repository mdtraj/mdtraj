import numpy as np
import pytest

import mdtraj as md

# Import the SAXS API and executable discovery
from mdtraj.saxs import compute_saxs_pepsi
from mdtraj.saxs.saxs import PEPSI_SAXS, find_executable
from mdtraj.testing import eq

example_trajectory = "2DFC_traj.pdb"  # 10 frames of MD simulation of 2DFC
example_saxs_data = "2DFC_SAXS.dat"  # SAXS data from https://sastutorials.org/Proteins/Proteins.html


@pytest.mark.skipif(
    not find_executable(PEPSI_SAXS),
    reason="Pepsi-SAXS binary not found",
)
def test_pepsi_saxs_basic(get_fn):
    traj = md.load(get_fn(example_trajectory))
    exp_path = get_fn(example_saxs_data)
    results, chi2 = compute_saxs_pepsi(traj, saxs_data=exp_path, num_frames=10)

    # Basic shape and content checks
    print(results)
    eq(set(results.columns), {"q", "exp_I", "exp_sigma", "fit_I"})
    assert len(results) > 0, "No q-points parsed from SAXS file."

    # q should be increasing (there are also some repeated values in the experimental data)
    assert np.all(np.diff(results["q"].values) >= 0)

    # chi2 value of trajectory using our double pass scheme
    # Running Pepsi-SAXS --cst on just the 2DFC.pdb structure file
    # gives the comparable value of 1.13
    eq(float(chi2), 1.08, decimal=2)


@pytest.mark.skipif(
    not find_executable(PEPSI_SAXS),
    reason="Pepsi-SAXS binary not found",
)
def test_pepsi_saxs_num_frames_clamp(get_fn):
    traj = md.load(get_fn(example_trajectory))

    exp_path = get_fn(example_saxs_data)

    # Request more frames than available; function should clamp internally and succeed
    results, chi2 = compute_saxs_pepsi(traj, saxs_data=exp_path, num_frames=1000)

    print(results)
    eq(set(results.columns), {"q", "exp_I", "exp_sigma", "fit_I"})
    assert len(results) > 0
    assert np.isfinite(chi2)


@pytest.mark.skipif(
    not find_executable(PEPSI_SAXS),
    reason="Pepsi-SAXS binary not found",
)
def test_pepsi_saxs_output_directory(get_fn, tmp_path):
    # Use a multi-frame trajectory
    traj = md.load(get_fn(example_trajectory))
    exp_path = get_fn(example_saxs_data)

    # Provide an explicit output directory; ensure files are created there
    outdir = tmp_path / "saxs_output"
    results, chi2 = compute_saxs_pepsi(traj, saxs_data=exp_path, num_frames=10, output_path=str(outdir))

    print(results)
    assert outdir.exists()
    # Spot check a few expected intermediate files
    # Note: exact file names come from the implementation.
    assert (outdir / "frame_0.pdb").exists()
    assert (outdir / "Ave_AA_frame0.fit").exists()

    assert len(results) > 0
    assert np.isfinite(chi2)


def test_pepsi_saxs_bad_data_path_raises(get_fn):
    traj = md.load(get_fn(example_trajectory))

    # Non-existent data path should raise OSError
    with pytest.raises(OSError):
        compute_saxs_pepsi(traj, saxs_data="does_not_exist.dat", num_frames=1)


@pytest.mark.skipif(
    not find_executable(PEPSI_SAXS),
    reason="Pepsi-SAXS binary not found",
)
def test_pepsi_saxs_invalid_num_frames_raises(get_fn):
    traj = md.load(get_fn(example_trajectory))

    # num_frames <= 0 should raise ValueError
    with pytest.raises(ValueError):
        compute_saxs_pepsi(traj, saxs_data=get_fn(example_saxs_data), num_frames=0)
