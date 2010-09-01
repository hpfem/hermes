from hermes1d.fekete.fekete import Mesh1D, Function

def test1():
    m = Mesh1D((-5, -4, 3, 10), (1, 5, 1))

def test2():
    eps = 1e-12
    func = lambda x: x**2
    f = Function(func, Mesh1D((-5, -4, 3, 10), (2, 5, 2)))
    for x in [-5, -4.5, -4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3, 4, 5, 6, 7, 10]:
        assert abs(f(x) - func(x)) < eps

    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 2)))
    for x in [-5, -4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3, 4, 5, 6, 7, 10]:
        assert abs(f(x) - func(x)) < eps
    x = -4.9
    assert abs(f(x) - func(x)) > 0.08
    x = -4.5
    assert abs(f(x) - func(x)) > 0.24

    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-5, -4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3, 10]:
        assert abs(f(x) - func(x)) < eps
    x = -4.9
    assert abs(f(x) - func(x)) > 0.08
    x = -4.5
    assert abs(f(x) - func(x)) > 0.24
    x = 4
    assert abs(f(x) - func(x)) > 5.9
    x = 5
    assert abs(f(x) - func(x)) > 9.9
    x = 6
    assert abs(f(x) - func(x)) > 11.9
    x = 7
    assert abs(f(x) - func(x)) > 11.9
    x = 8
    assert abs(f(x) - func(x)) > 9.9
    x = 9
    assert abs(f(x) - func(x)) > 5.9

def test3():
    eps = 1e-12
    func = lambda x: x**2
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3]:
        assert abs(f(x) - func(x)) < eps

    func = lambda x: x**3
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3]:
        assert abs(f(x) - func(x)) < eps

    func = lambda x: x**4
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3]:
        assert abs(f(x) - func(x)) < eps

    func = lambda x: x**5
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    for x in [-4, -3, -2, -1, 0, 0.01, 1e-5, 1, 2, 3]:
        assert abs(f(x) - func(x)) < eps

    func = lambda x: x**6
    f = Function(func, Mesh1D((-5, -4, 3, 10), (1, 5, 1)))
    x = -1
    assert abs(f(x) - func(x)) > 61.9
    x = 0
    assert abs(f(x) - func(x)) > 61.9
    x = 1
    assert abs(f(x) - func(x)) > 61.6
    x = 2
    assert abs(f(x) - func(x)) > 28.9

def test4():
    eps = 1e-12
    func = lambda x: x**2
    orig_mesh = Mesh1D((-5, -4, 3, 10), (1, 5, 1))
    mesh1     = Mesh1D((-5, -4, 3, 10), (1, 1, 1))
    f = Function(func, orig_mesh)
    g = f.project_onto(mesh1)
    h = Function(func, mesh1)
    assert g == Function(func, mesh1)
    assert h == h.project_onto(orig_mesh)

def test5():
    eps = 1e-12
    func = lambda x: x**2
    mesh1 = Mesh1D((-5, -4, 3, 10), (2, 5, 2))
    mesh2 = Mesh1D((-5, -4, 3, 10), (2, 2, 2))
    mesh3 = Mesh1D((-5, -4, 3, 10), (2, 2, 1))
    mesh4 = Mesh1D((-5, 10), (2,))
    mesh5 = Mesh1D((-5, 10), (3,))
    mesh6 = Mesh1D((-5, 10), (1,))
    f = Function(func, mesh1)
    g = Function(func, mesh2)
    h = Function(func, mesh3)
    l = Function(func, mesh4)

    assert f == g
    assert g == f
    assert f == l
    assert g == l
    assert f != h
    assert h != f
    assert g != h
    assert h != g

    assert f == Function(lambda x: x**2, mesh1)
    assert f != Function(lambda x: x**3, mesh1)
    assert f == Function(lambda x: x**2, mesh2)
    assert f == Function(lambda x: x**2, mesh4)
    assert f == Function(lambda x: x**2, mesh5)
    assert f != Function(lambda x: x**2, mesh6)

def test6():
    mesh1 = Mesh1D((-5, -4, 3, 10), (2, 5, 2))
    mesh2 = Mesh1D((-5, -4, 3, 10), (2, 2, 2))
    mesh3 = Mesh1D((-5, -4, 3, 10), (2, 2, 1))
    mesh4 = Mesh1D((-5, 10), (2,))
    mesh5 = Mesh1D((-5, 10), (3,))
    mesh6 = Mesh1D((-5, 10), (1,))
    mesh7 = Mesh1D((-5, 10), (1,))
    mesh8 = Mesh1D((-5, 0, 10), (1, 4))

    assert mesh1 == mesh1
    assert not (mesh1 != mesh1)
    assert mesh1 != mesh2
    assert mesh1 != mesh3
    assert mesh1 != mesh4
    assert mesh1 != mesh5
    assert mesh1 != mesh6
    assert mesh6 == mesh7
    assert mesh1.union(mesh1) == mesh1

    assert mesh1.union(mesh2) == mesh1
    assert mesh2.union(mesh1) == mesh1

    assert mesh1.union(mesh3) == mesh1
    assert mesh3.union(mesh1) == mesh1

    assert mesh1.union(mesh4) == mesh1
    assert mesh4.union(mesh1) == mesh1

    assert mesh1.union(mesh5) == Mesh1D((-5, -4, 3, 10), (3, 5, 3))
    assert mesh5.union(mesh1) == Mesh1D((-5, -4, 3, 10), (3, 5, 3))

    assert mesh1.union(mesh6) == mesh1
    assert mesh6.union(mesh1) == mesh1

    assert mesh1.union(mesh8) == Mesh1D((-5, -4, 0, 3, 10), (2, 5, 5, 4))
    assert mesh8.union(mesh1) == Mesh1D((-5, -4, 0, 3, 10), (2, 5, 5, 4))

def test7():
    mesh1 = Mesh1D((-5, -4, 3, 10), (2, 5, 2))
    mesh2 = Mesh1D((-5, -4, 3, 10), (2, 2, 2))
    mesh3 = Mesh1D((-5, -4, 3, 10), (2, 2, 1))
    mesh4 = Mesh1D((-5, 10), (2,))
    mesh5 = Mesh1D((-5, 10), (3,))
    mesh6 = Mesh1D((-5, 10), (1,))
    mesh8 = Mesh1D((-5, 0, 10), (1, 4))

    assert mesh1.restrict_to_interval(-5, 10) == mesh1
    assert mesh1.restrict_to_interval(-4.5, 10) == Mesh1D((-4.5, -4, 3, 10),
            (2, 5, 2))
    assert mesh1.restrict_to_interval(-4, 10) != mesh1
    assert mesh1.restrict_to_interval(-4, 10) == Mesh1D((-4, 3, 10), (5, 2))
    assert mesh1.restrict_to_interval(-3.5, 10) == Mesh1D((-3.5, 3, 10), (5, 2))
    assert mesh1.restrict_to_interval(3, 10) == Mesh1D((3, 10), (2,))
    assert mesh1.restrict_to_interval(3.5, 10) == Mesh1D((3.5, 10), (2,))

    assert mesh2.restrict_to_interval(-5, 10) == mesh2
    assert mesh2.restrict_to_interval(-4.5, 10) == Mesh1D((-4.5, -4, 3, 10),
            (2, 2, 2))
    assert mesh2.restrict_to_interval(-4, 10) != mesh2
    assert mesh2.restrict_to_interval(-4, 10) == Mesh1D((-4, 3, 10), (2, 2))
    assert mesh2.restrict_to_interval(-3.5, 10) == Mesh1D((-3.5, 3, 10), (2, 2))
    assert mesh2.restrict_to_interval(3, 10) == Mesh1D((3, 10), (2,))
    assert mesh2.restrict_to_interval(3.5, 10) == Mesh1D((3.5, 10), (2,))

    assert mesh3.restrict_to_interval(-5, 10) == mesh3
    assert mesh3.restrict_to_interval(-4.5, 10) == Mesh1D((-4.5, -4, 3, 10),
            (2, 2, 1))
    assert mesh3.restrict_to_interval(-4, 10) != mesh3
    assert mesh3.restrict_to_interval(-4, 10) == Mesh1D((-4, 3, 10), (2, 1))
    assert mesh3.restrict_to_interval(-3.5, 10) == Mesh1D((-3.5, 3, 10), (2, 1))
    assert mesh3.restrict_to_interval(3, 10) == Mesh1D((3, 10), (1,))
    assert mesh3.restrict_to_interval(3.5, 10) == Mesh1D((3.5, 10), (1,))
