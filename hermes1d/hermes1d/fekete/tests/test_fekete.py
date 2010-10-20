from math import sin, cos, log, pi, sqrt, e, log, exp

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

def test8():
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
    zero = Function(lambda x: 0., Mesh1D((-5, 10), (1,)))
    assert zero.l2_norm() < eps
    assert Function(lambda x: 0., mesh1) == zero
    assert Function(lambda x: 0., mesh2) == zero
    assert Function(lambda x: 0., mesh3) == zero
    assert Function(lambda x: 0., mesh4) == zero
    assert Function(lambda x: 0., mesh5) == zero
    assert Function(lambda x: 0., mesh6) == zero

    assert f - g == zero
    assert (f-g).l2_norm() < eps
    assert g - f == zero
    assert (g - f).l2_norm() < eps
    assert f - l == zero
    assert (f - l).l2_norm() < eps
    assert g - l == zero
    assert (g - l).l2_norm() < eps
    assert f - h != zero
    assert (f - h).l2_norm() > eps
    assert h - f != zero
    assert (h - f).l2_norm() > eps
    assert g - h != zero
    assert (g - h).l2_norm() > eps
    assert h - g != zero
    assert (h - g).l2_norm() > eps

    assert f - Function(lambda x: x**2, mesh1) == zero
    assert f - Function(lambda x: x**3, mesh1) != zero
    assert f - Function(lambda x: x**2, mesh2) == zero
    assert f - Function(lambda x: x**2, mesh4) == zero
    assert f - Function(lambda x: x**2, mesh5) == zero
    assert f - Function(lambda x: x**2, mesh6) != zero

def test_l2_h1_1():
    eps = 1e-12
    eps_low = 1e-10
    func = lambda x: sin(x)
    mesh1 = Mesh1D((0, pi), (20,))
    mesh2 = Mesh1D((0, pi/2, pi), (20, 20))
    mesh3 = Mesh1D((0, pi/2), (20,))
    f = Function(func, mesh1)
    g = Function(func, mesh2)
    h = Function(func, mesh3)
    assert abs(f.l2_norm(method="Fekete")-sqrt(pi/2)) < eps
    assert abs(g.l2_norm(method="Fekete")-sqrt(pi/2)) < eps
    assert abs(h.l2_norm(method="Fekete")-sqrt(pi/4)) < eps

    assert abs(f.l2_norm(method="FE")-sqrt(pi/2)) < eps
    assert abs(g.l2_norm(method="FE")-sqrt(pi/2)) < eps
    assert abs(h.l2_norm(method="FE")-sqrt(pi/4)) < eps

    assert abs(f.h1_norm()-sqrt(pi)) < eps_low
    assert abs(g.h1_norm()-sqrt(pi)) < eps_low
    assert abs(h.h1_norm()-sqrt(pi/2)) < eps_low

    func = lambda x: cos(x)
    mesh1 = Mesh1D((0, pi/4, pi/2, 3*pi/4, pi), (20, 20, 20, 20))
    f = Function(func, mesh1)
    assert abs(f.l2_norm(method="Fekete")-sqrt(pi/2)) < eps
    assert abs(f.l2_norm(method="FE")-sqrt(pi/2)) < eps

    assert abs(f.h1_norm()-sqrt(pi)) < eps

def test_l2_h1_2():
    eps = 1e-9
    func = lambda x: log(x)
    mesh1 = Mesh1D((1, 1.5, 2, 2.5, e), (20, 20, 20, 20))
    f = Function(func, mesh1)
    l2_norm_exact = sqrt(e-2)
    h1_norm_exact = sqrt(e-1-exp(-1))
    assert abs(f.l2_norm(method="Fekete")-l2_norm_exact) < eps
    assert abs(f.l2_norm(method="FE")-l2_norm_exact) < eps
    assert abs(f.h1_norm()-h1_norm_exact) < eps

def test_power():
    eps = 1e-12
    func = lambda x: x
    mesh1 = Mesh1D((0, 1), (1,))
    mesh2 = Mesh1D((0, 1), (2,))
    mesh3 = Mesh1D((0, 1), (3,))
    f = Function(func, mesh1)
    assert abs(f.l2_norm() - sqrt(1./3)) < eps
    assert f**2 != Function(lambda x: x**2, mesh1)
    assert f**2 == Function(lambda x: x**2, mesh2)
    assert f**2 == Function(lambda x: x**2, mesh3)

    func = lambda x: x
    mesh1 = Mesh1D((5, 6), (1,))
    mesh2 = Mesh1D((5, 6), (2,))
    mesh3 = Mesh1D((5, 6), (3,))
    f = Function(func, mesh1)
    assert f**2 != Function(lambda x: x**2, mesh1)
    assert f**2 == Function(lambda x: x**2, mesh2)
    assert f**2 == Function(lambda x: x**2, mesh3)

    func = lambda x: x**3+x
    mesh1 = Mesh1D((5, 6), (3,))
    mesh2 = Mesh1D((5, 6), (5,))
    mesh3 = Mesh1D((5, 6), (6,))
    mesh4 = Mesh1D((5, 6), (9,))
    f = Function(func, mesh1)
    assert f**2 != Function(lambda x: x**6+2*x**4+x**2, mesh1)
    assert f**2 != Function(lambda x: x**6+2*x**4+x**2, mesh2)
    assert f**2 == Function(lambda x: x**6+2*x**4+x**2, mesh3)
    assert f**3 == Function(lambda x: x**3 + 3*x**5 + 3*x**7 + x**9, mesh4)

    func = lambda x: x**3+x
    mesh1 = Mesh1D((5, 5.1, 6), (3, 3))
    mesh2 = Mesh1D((5, 5.1, 6), (5, 5))
    mesh3 = Mesh1D((5, 5.1, 6), (6, 6))
    mesh4 = Mesh1D((5, 5.1, 6), (9, 9))
    mesh5 = Mesh1D((5, 5.1, 6), (9, 7))
    mesh6 = Mesh1D((5, 5.1, 6), (6, 5))
    mesh7 = Mesh1D((5, 6), (6,))
    mesh8 = Mesh1D((5, 6), (9,))
    f = Function(func, mesh1)
    assert f**2 != Function(lambda x: x**6+2*x**4+x**2, mesh1)
    assert f**2 != Function(lambda x: x**6+2*x**4+x**2, mesh2)
    assert f**2 != Function(lambda x: x**6+2*x**4+x**2, mesh6)
    assert f**2 == Function(lambda x: x**6+2*x**4+x**2, mesh3)
    assert f**3 == Function(lambda x: x**3 + 3*x**5 + 3*x**7 + x**9, mesh4)
    assert f**3 != Function(lambda x: x**3 + 3*x**5 + 3*x**7 + x**9, mesh5)

    assert f**2 == Function(lambda x: x**6+2*x**4+x**2, mesh7)
    assert f**3 == Function(lambda x: x**3 + 3*x**5 + 3*x**7 + x**9, mesh8)
