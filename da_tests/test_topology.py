import pytest
import mdtraj as md


def test_delete_atom_by_index_deprecated(get_fn):
    top = md.load(get_fn("ala_ala_ala.pdb")).topology
    with pytest.deprecated_call(): 
        top.delete_atom_by_index(0)