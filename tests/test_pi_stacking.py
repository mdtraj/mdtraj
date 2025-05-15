##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
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
##############################################################################

import mdtraj as md
from mdtraj.geometry.pi_stacking import pi_stacking


def test_pi_stacking(get_fn):
    """Test that pi-stacking is being detected correctly."""
    t: md.Trajectory = md.load(get_fn("4BFQ.pdb"))
    lig_aromatic_atms_to_test = [
        ("C16", "C17", "C18", "C19", "C20", "C21"),  # pi-stacking
        ("N4", "C9", "C15", "C16", "C21", "C8"),  # pi-stacking
        ("N5", "C10", "C11", "C12", "C13", "C14"),  # NOT pi-stacking
    ]
    lig_grps = []
    for lig_grp in lig_aromatic_atms_to_test:
        lig_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 5 and resSeq 301 and (name {lig_grp[0]} or name {' or name '.join(lig_grp[1:])})",
                )
            ),
        )
    protein_grps = []
    # TYR186
    protein_stacking_atms = [("CE2", "CD2", "CG", "CD1", "CE1", "CZ")]
    for protein_grp in protein_stacking_atms:
        protein_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 0 and resSeq 186 and (name {protein_grp[0]} or name {' or name '.join(protein_grp[1:])})",  # noqa
                )
            ),
        )
    stacking_interactions = pi_stacking(t, lig_grps, protein_grps)
    # Since it returns per-frame, just get the 'first'
    stacking_interactions = stacking_interactions[0]
    assert len(stacking_interactions) == 2
    assert (lig_grps[0], protein_grps[0]) in stacking_interactions
    assert (lig_grps[1], protein_grps[0]) in stacking_interactions
    assert (lig_grps[2], protein_grps[0]) not in stacking_interactions


def test_t_stacking_more_relaxed_radius(get_fn):
    """Test that t-stacking version of pi-stacking is being detected correctly."""
    t: md.Trajectory = md.load(get_fn("4BFQ.pdb"))
    lig_aromatic_atms_to_test = [
        ("C2", "C3", "C4", "N1", "C6", "N6"),  # t-stacking
        ("N5", "C10", "C11", "C12", "C13", "C14"),  # NOT pi-stacking
    ]
    lig_grps = []
    for lig_grp in lig_aromatic_atms_to_test:
        lig_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 5 and resSeq 302 and (name {lig_grp[0]} or name {' or name '.join(lig_grp[1:])})",
                )
            ),
        )
    protein_grps = []
    # TYR91
    protein_stacking_atms = [("CE2", "CD2", "CG", "CD1", "CE1", "CZ")]
    for protein_grp in protein_stacking_atms:
        protein_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 0 and resSeq 91 and (name {protein_grp[0]} or name {' or name '.join(protein_grp[1:])})",
                )
            ),
        )
    # More permissive edge intersection radius for t-stack
    stacking_interactions = pi_stacking(t, lig_grps, protein_grps, edge_intersection_radius=0.21)
    stacking_interactions = stacking_interactions[0]
    assert len(stacking_interactions) == 1
    assert (lig_grps[0], protein_grps[0]) in stacking_interactions
    assert (lig_grps[1], protein_grps[0]) not in stacking_interactions


def test_t_stacking_default_settings(get_fn):
    """Test that t-stacking version of pi-stacking is being detected correctly."""

    t: md.Trajectory = md.load(get_fn("6A22.pdb"))
    lig_aromatic_atms_to_test = [
        ("S10", "C9", "C8", "N7", "C6"),  # t-stacking
        ("C14", "C15", "C16", "C17", "C18", "C19"),  # NOT pi-stacking
    ]
    lig_grps = []
    for lig_grp in lig_aromatic_atms_to_test:
        lig_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 9 and resSeq 9000 and (name {lig_grp[0]} or name {' or name '.join(lig_grp[1:])})",
                )
            ),
        )
    protein_grps = []
    # PHR118/378
    protein_stacking_atms = [("CE2", "CD2", "CG", "CD1", "CE1", "CZ")]
    for protein_grp in protein_stacking_atms:
        protein_grps.append(
            tuple(
                int(idx)
                for idx in t.top.select(
                    f"chainid 2 and resSeq 378 and (name {protein_grp[0]} or name {' or name '.join(protein_grp[1:])})",  # noqa
                )
            ),
        )
    stacking_interactions = pi_stacking(t, lig_grps, protein_grps)
    stacking_interactions = stacking_interactions[0]
    assert len(stacking_interactions) == 1
    assert (lig_grps[0], protein_grps[0]) in stacking_interactions
    assert (lig_grps[1], protein_grps[0]) not in stacking_interactions
