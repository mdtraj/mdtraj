from mdtraj.pdb import element
from mdtraj.testing import assert_raises

def test_immutable():
    def f():
        element.hydrogen.mass = 1
    def g():
        element.radium.symbol = 'sdfsdfsdf'
    def h():
        element.iron.name = 'sdfsdf'

    assert_raises(AttributeError, f)
    assert_raises(AttributeError, g)
    assert_raises(AttributeError, h)
    assert element.hydrogen.mass == 1.007947
    assert element.radium.symbol == 'Ra'
    assert element.iron.name == 'iron'
