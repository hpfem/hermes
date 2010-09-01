from numpy import array

from hermes1d import Mesh

eps = 1e-12

def test_mesh1():
    m = Mesh(-1, 1, 5, 2)
    pts, p = m.get_mesh_data()
    assert (abs(pts - array([-1, -0.6, -0.2, 0.2, 0.6, 1])) < eps).all()
    assert (abs(p - array([2, 2, 2, 2, 2])) < eps).all()

    m = Mesh(-1, 1, 4, 1)
    pts, p = m.get_mesh_data()
    assert (abs(pts - array([-1, -0.5, 0, 0.5, 1])) < eps).all()
    assert (abs(p - array([1, 1, 1, 1])) < eps).all()

def test_mesh2():
    m = Mesh(-1, 1, 5, 2, 1)
    pts, p = m.get_mesh_data()
    assert (abs(pts - array([-1, -0.6, -0.2, 0.2, 0.6, 1])) < eps).all()
    assert (abs(p - array([2, 2, 2, 2, 2])) < eps).all()

def test_mesh3():
    m = Mesh([0, 1, 3, 4], [2, 2, 2])
    pts, p = m.get_mesh_data()
    assert (abs(pts - array([0, 1, 3, 4])) < eps).all()
    assert (abs(p - array([2, 2, 2])) < eps).all()

    m = Mesh([-1, -0.5, 0, 0.5, 1], [2, 2, 2, 2])
    pts, p = m.get_mesh_data()
    assert (abs(pts - array([-1, -0.5, 0, 0.5, 1])) < eps).all()
    assert (abs(p - array([2, 2, 2, 2])) < eps).all()
