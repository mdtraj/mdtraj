from mdtraj.pdb import element
from mdtraj.testing import assert_raises

def test_immutable():
    with assert_raises(AttributeError):
        element.hydrogen.mass = 1
    with assert_raises(AttributeError):
        element.radium.symbol = 'sdfsdfsdf'
    with assert_raises(AttributeError):
        element.iron.name = 'sdfsdf'

    assert element.hydrogen.mass == 1.007947
    assert element.radium.symbol == 'Ra'
    assert element.iron.name == 'iron'
