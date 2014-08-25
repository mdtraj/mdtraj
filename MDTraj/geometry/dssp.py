import numpy as np
import sys; sys.path.insert(0, '/Users/rmcgibbo/projects/mdtraj/build/lib.macosx-10.5-x86_64-2.7/')

import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry.hbond import _prep_kabsch_sander_arrays
from mdtraj.geometry import _geometry


def compute_dssp(traj):
    """compute_drid(traj)

    Compute Dictionary of protein secondary structure (DSSP) secondary
    structure assignments for each frame of this trajectory
    
    Parameters
    ----------
    traj : md.Trajectory
        A trajectory

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
   
    _geometry._dssp(xyz, nco_indices, ca_indices, proline_indices, chain_ids)
   
    
if __name__ == '__main__':
    from mdtraj.testing import get_fn
    t = md.load(get_fn('1bpi.pdb'))[0]
    compute_dssp(t)
    