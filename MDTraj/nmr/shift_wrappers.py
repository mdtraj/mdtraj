import string
import shutil
import pandas as pd
import os
import contextlib
import tempfile


@contextlib.contextmanager
def enter_temp_directory():
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield
    os.chdir(cwd)
    shutil.rmtree(temp_dir)


def chemical_shifts_shiftx2(trj):
    """Predict chemical shifts of a trajectory using ShiftX2.
    
    Parameters
    ----------
    trj : Trajectory
        Trajectory to predict shifts for.
        
    Returns
    -------
    
    results : Pandas DataFrame
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
    
    results = []
    with enter_temp_directory():
        for i in range(trj.n_frames):        
            trj[i].save("./trj%d.pdb" % i)
        cmd = "shiftx2.py -b 'trj*.pdb'"
        
        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your installation." % cmd))

        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d.pdb.cs" % i)
            d.rename(columns={"NUM":"resSeq", "RES":"resName", "ATOMNAME":"name"}, inplace=True)
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
    
    results : Pandas DataFrame
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
    
    with enter_temp_directory():
        trj.save("./trj.pdb")
        cmd = "ppm_linux_64.exe -pdb trj.pdb -mode detail"
        
        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your installation." % cmd))

        d = pd.read_csv("./bb_details.dat", delim_whitespace=True)
        columns = ["resSeq", "resName", "name", "expt", "other"]

        d = pd.read_csv("./bb_details.dat", delim_whitespace=True, header=None).drop([0, 4], axis=1)
        d = d.rename(columns={1:"resSeq", 2:"resName", 3:"name"})
        d["resSeq"] += first_resSeq - 1  # Fix bug in PPM that reindexes to 1
        d = d.drop("resName", axis=1)
        d = d.set_index(["resSeq", "name"])
        d.columns = trj.time
    
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
    
    results : Pandas DataFrame
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

    names = ["VARS", "resSeq", "resName", "name", "SS_SHIFT", "SHIFT", "RC_SHIFT", "HM_SHIFT", "EF_SHIFT", "SIGMA"]
    
    with enter_temp_directory():
        for i in range(trj.n_frames):
            trj[i].save("./trj%d.pdb" % i)
        
        cmd = "SPARTA+.linux -in %s" % string.join(["trj%d.pdb" % i for i in range(trj.n_frames)])

        return_flag = os.system(cmd)

        if return_flag != 0:
            raise(IOError("Could not successfully execute command '%s', check your installation." % cmd))

        lines_to_skip = _get_lines_to_skip("trj0_pred.tab")
        
        results = []
        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d_pred.tab" % i, skiprows=lines_to_skip, delim_whitespace=True, header=None, names=names)
            d["frame"] = trj.time[i]
            results.append(d)
        
    results = pd.concat(results)
    results = results.pivot_table(rows=["resSeq", "name"], cols="frame", values="SHIFT")
    
    return results
