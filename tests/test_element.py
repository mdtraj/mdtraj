import mdtraj as md
import pytest
import pickle
from mdtraj import element
from mdtraj.testing import eq


def test_immutable():
    def f():
        element.hydrogen.mass = 1

    def g():
        element.radium.symbol = 'sdfsdfsdf'

    def h():
        element.iron.name = 'sdfsdf'

    pytest.raises(AttributeError, f)
    pytest.raises(AttributeError, g)
    pytest.raises(AttributeError, h)
    assert element.hydrogen.mass == 1.007947
    assert element.radium.symbol == 'Ra'
    assert element.iron.name == 'iron'


def test_element_0(get_fn):
    t = md.load(get_fn('bpti.pdb'))

    a = t.top.atom(15)
    H = element.Element.getBySymbol('H')

    eq(a.element, element.hydrogen)


def test_element_pickle():
    """Test that every Element object can pickle and de-pickle"""
    for el in dir(element):
        if isinstance(el, element.Element):
            assert el == pickle.loads(pickle.dumps(el))
