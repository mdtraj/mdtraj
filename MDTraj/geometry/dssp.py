import numpy as np
import sys; sys.path.insert(0, '/Users/rmcgibbo/projects/mdtraj/build/lib.macosx-10.5-x86_64-2.7/')

import mdtraj as md
from mdtraj.utils.six import PY2
from mdtraj.utils import ensure_type
from mdtraj.geometry.hbond import _prep_kabsch_sander_arrays
from mdtraj.geometry import _geometry

__all__ = ['compute_dssp']

def compute_dssp(traj):
    """compute_drid(traj)

    Compute Dictionary of protein secondary structure (DSSP) secondary
    structure assignments for each frame of this trajectory
    
    Parameters
    ----------
    traj : md.Trajectory
        A trajectory

    Returns
    -------
    assignments : list of str
        The return value is a list of strings, of length n_frames. Each string
        is of length n_residues, and contains the secondary-structure code of
        each residue.

    Notes
    -----
    The DSSP assignment codes are:

       - 'H' : Alpha helix
       - 'B' : Residue in isolated beta-bridge
       - 'E' : Extended strand, participates in beta ladder
       - 'G' : 3-helix (3/10 helix)
       - 'I' : 5 helix (pi helix)
       - 'T' : hydrogen bonded turn
       - 'S' : bend
       - ' ' : Loops and irregular elements

    References
    ----------
    .. [1] Kabsch W, Sander C (1983). "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features". Biopolymers 22 (12): 2577-637. dio:10.1002/bip.360221211
    """
    if traj.topology is None:
        raise ValueError('kabsch_sander requires topology')
    if not _geometry._processor_supports_sse41():
        raise RuntimeError('This CPU does not support the required instruction set (SSE4.1)')
    
    xyz, nco_indices, ca_indices, proline_indices = _prep_kabsch_sander_arrays(traj)
    chain_ids = np.array([r.chain.index for r in traj.top.residues], dtype=np.int32)

    return _geometry._dssp(xyz, nco_indices, ca_indices, proline_indices, chain_ids)
