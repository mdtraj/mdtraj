##############################################################################
# Imports
##############################################################################

import numpy as np

import mdtraj as md
from mdtraj.utils.six import PY2
from mdtraj.utils import ensure_type
from mdtraj.geometry.hbond import _prep_kabsch_sander_arrays
from mdtraj.geometry import _geometry
if PY2:
    from string import maketrans
else:
    maketrans = str.maketrans

##############################################################################
# GLOBALS
##############################################################################

SIMPLIFIED_CODE_TRANSLATION = maketrans('HGIEBTS ', 'HHHEECCC')
__all__ = ['compute_dssp']


##############################################################################
# CODE
##############################################################################


def compute_dssp(traj, simplified=True):
    """Compute Dictionary of protein secondary structure (DSSP) secondary structure assignments
    
    Parameters
    ----------
    traj : md.Trajectory
        A trajectory

    Returns
    -------
    assignments : np.ndarray, shape=(n_frames, n_residues), dtype=S1
        The assignments is a 2D array of character codes (see below), giving
        the secondary structure of each residue in each frame.
    simplified  : bool, default=True.
        Use the simplified 3-category assignment scheme. Otherwise the original
        8-category scheme is used.

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

    The simplified DSSP codes are:

       - 'H' : Helix. Either of the 'H', 'G', or 'I' codes.
       - 'E' : Strand. Either of the 'E', or 'B' codes.
       - 'C' : Coil. Either of the 'T', 'S' or ' ' codes.

    Our implementation is based on DSSP-2.2.0, written by Maarten L. Hekkelman
    and distributed under the Boost Software license.

    References
    ----------
    .. [1] Kabsch W, Sander C (1983). "Dictionary of protein secondary
       structure: pattern recognition of hydrogen-bonded and geometrical
       features". Biopolymers 22 (12): 2577-637. dio:10.1002/bip.360221211
    """
    if traj.topology is None:
        raise ValueError('kabsch_sander requires topology')
    if not _geometry._processor_supports_sse41():
        raise RuntimeError('This CPU does not support the required instruction set (SSE4.1)')
    
    xyz, nco_indices, ca_indices, proline_indices = _prep_kabsch_sander_arrays(traj)
    chain_ids = np.array([r.chain.index for r in traj.top.residues], dtype=np.int32)

    value = _geometry._dssp(xyz, nco_indices, ca_indices, proline_indices, chain_ids)
    if simplified:
        value = value.translate(SIMPLIFIED_CODE_TRANSLATION)

    n_frames = xyz.shape[0]
    n_residues = nco_indices.shape[0]
    as2darray = np.fromstring(value, dtype=np.dtype('S1')).reshape(n_frames, n_residues)
    return as2darray