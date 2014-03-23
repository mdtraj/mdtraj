import string
import shutil
import pandas as pd
import os
import contextlib
import tempfile


@contextlib.contextmanager
def enter_temp_directory():
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield
    os.chdir(cwd)
    shutil.rmtree(temp_dir)


def contains_only_protein(trj):
    return True


def chemical_shifts_shiftx2(trj):
    results = []
    with enter_temp_directory():
        for i in range(trj.n_frames):        
            trj[i].save("./trj%d.pdb" % i)
        cmd = "shiftx2.py -b 'trj*.pdb'"
        os.system(cmd)
        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d.pdb.cs" % i)
            d.rename(columns={"NUM":"resSeq", "RES":"resName", "ATOMNAME":"name"}, inplace=True)
            d["frame"] = trj.time[i]
            results.append(d)

    results = pd.concat(results)
    results = results.pivot_table(rows=["resSeq", "name"], cols="frame", values="SHIFT")
    return results


def chemical_shifts_ppm(trj):
    with enter_temp_directory():
        trj.save("./trj.pdb")
        cmd = "ppm_linux_64.exe -pdb trj.pdb -mode detail"
        os.system(cmd)
        d = pd.read_csv("./bb_details.dat", delim_whitespace=True)
        columns = ["resSeq", "resName", "name", "expt", "other"]

        d = pd.read_csv("./bb_details.dat", delim_whitespace=True, header=None).drop([0, 4], axis=1)
        d = d.rename(columns={1:"resSeq", 2:"resName", 3:"name"})
        d = d.drop("resName", axis=1)
        d = d.set_index(["resSeq", "name"])
        d.columns = trj.time
    
    return d


def get_lines_to_skip(filename):
    format_string = """FORMAT %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"""
    handle = open(filename)
    for i, line in enumerate(handle):
        if line.find(format_string) != -1:
            return i + 2
    
    raise(Exception("No format string found in SPARTA+ file!"))


def chemical_shifts_spartaplus(trj):
    #names = ["VARS", "RESID", "RESNAME", "ATOMNAME", "SS_SHIFT", "SHIFT", "RC_SHIFT", "HM_SHIFT", "EF_SHIFT", "SIGMA"]  # Names in actual file
    names = ["VARS", "resSeq", "resName", "name", "SS_SHIFT", "SHIFT", "RC_SHIFT", "HM_SHIFT", "EF_SHIFT", "SIGMA"]
    
    first_resSeq = trj.top.residue(0).resSeq
    
    with enter_temp_directory():
        for i in range(trj.n_frames):
            trj_i = trj[i].copy()
            trj[i].save("./trj%d.pdb" % i)
        
        cmd = "SPARTA+.linux -in %s" % string.join(["trj%d.pdb" % i for i in range(trj.n_frames)])
        os.system(cmd)
        lines_to_skip = get_lines_to_skip("trj0_pred.tab")
        
        results = []
        for i in range(trj.n_frames):
            d = pd.read_csv("./trj%d_pred.tab" % i, skiprows=lines_to_skip, delim_whitespace=True, header=None, names=names)
            d["frame"] = trj.time[i]
            results.append(d)
        
    results = pd.concat(results)
    results = results.pivot_table(rows=["resSeq", "name"], cols="frame", values="SHIFT")
    
    return results
