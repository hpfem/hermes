from math import sin, cos, pi, atan2

from numpy import array, zeros, dot, eye
from numpy.linalg import inv, norm

from hermes2d._numerical_flux import (matrix_R, matrix_R_inv, matrix_D_minus,
        flux_riemann, flux_riemann_invert, numerical_flux, A_x as matrix_A_x,
        A_z as matrix_A_z, f_x as vector_f_x, f_z as vector_f_z, R as R_const,
        c_v)

eps = 1e-10

def R(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_R(i, j, *w)
    return A

def R_inv(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_R_inv(i, j, *w)
    return A

def D_minus(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_D_minus(i, j, *w)
    return A

def A_x(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_A_x(i, j, *w)
    return A

def A_z(w):
    A = zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j] = matrix_A_z(i, j, *w)
    return A

def f_x(w):
    A = zeros((4,))
    for i in range(4):
        A[i] = vector_f_x(i, *w)
    return A

def f_z(w):
    A = zeros((4,))
    for i in range(4):
        A[i] = vector_f_z(i, *w)
    return A

def A_minus(w):
    return dot(R(w), dot(D_minus(w), R_inv(w)))

def f_riemann(w_l, w_r):
    return f_x(w_l) + dot(A_minus(w_r), w_r) - dot(A_minus(w_l), w_l)

def T_rot(beta):
    # this is the 3D rotation matrix (2.2.51 in Pavel's master thesis)
    # in 2D without the middle column and row, alpha = 0
    alpha = 0
    return array([
        [1, 0, 0, 0],
        [0, cos(alpha)*cos(beta), sin(beta), 0],
        [0, -cos(alpha)*sin(beta), cos(beta), 0],
        [0, 0, 0, 1]
        ])

def calc_p(w):
    w0, w1, w3, w4 = w
    p = R_const/c_v * (w4 - (w1**2 + w3**2)/(2*w0))
    return p

def flux(w, n):
    return dot(A_x(w), w) * n[0] + dot(A_z(w), w) * n[1]

def tangent_w(w, n):
    _p = calc_p(w)
    vel = w[1:3]/w[0]
    vel = vel - dot(vel, n)*n
    w0 = w[0]
    w1, w3 = w[0]*vel
    w4 = w[3]
    w4 = _p * c_v / R_const + (w1*w1+w3*w3) / (2*w0);
    w = array([w0, w1, w3, w4])
    return w

# ---------------------------------------

def test_inv():
    w = array([1.1, -10, 13, 700.1])
    assert (abs(dot(R(w), R_inv(w))-eye(4)) < eps).all()
    assert (abs(R_inv(w) - inv(R(w))) < eps).all()
    w = array([1.1, 10, 13, 700.1])
    assert (abs(dot(R(w), R_inv(w))-eye(4)) < eps).all()
    assert (abs(R_inv(w) - inv(R(w))) < eps).all()
    w = array([3.1, 10, 13, 800.1])
    assert (abs(dot(R(w), R_inv(w))-eye(4)) < eps).all()
    assert (abs(R_inv(w) - inv(R(w))) < eps).all()

def test_flux():
    w_l = array([1.1, -10, 13, 700.1])
    w_r = array([1.1, -10, 13, 800.1])
    assert (abs(f_riemann(w_l, w_l) - f_x(w_l)) < eps).all()
    assert (abs(f_riemann(w_r, w_r) - f_x(w_r)) < eps).all()

    alpha = 0
    m = T_rot(alpha)
    w_r = dot(m, w_r)
    w_l = dot(m, w_l)
    m = T_rot(pi)
    flux1 = f_riemann(w_l, w_r)
    flux2 = -dot(inv(m), f_riemann(dot(m, w_r), dot(m, w_l)))
    flux3 = flux_riemann(w_l, w_r)
    flux4 = -dot(inv(m), flux_riemann(dot(m, w_r), dot(m, w_l)))
    flux5 = -flux_riemann_invert(w_l, w_r)
    assert (abs(flux1 - flux2) < eps).all()
    assert (abs(flux1 - flux3) < eps).all()
    assert (abs(flux1 - flux4) < eps).all()
    assert (abs(flux1 - flux5) < eps).all()

    w_l = array([1.1, -10, 13, 700.1])
    w_r = array([1.1, -10, 13, 800.1])
    alpha = 0.3
    m = T_rot(alpha)
    w_r = dot(m, w_r)
    w_l = dot(m, w_l)
    m = T_rot(pi)
    flux1 = f_riemann(w_l, w_r)
    flux2 = -dot(inv(m), f_riemann(dot(m, w_r), dot(m, w_l)))
    flux3 = flux_riemann(w_l, w_r)
    flux4 = -dot(inv(m), flux_riemann(dot(m, w_r), dot(m, w_l)))
    flux5 = -flux_riemann_invert(w_l, w_r)
    assert (abs(flux1 - flux2) < eps).all()
    assert (abs(flux1 - flux3) < eps).all()
    assert (abs(flux1 - flux4) < eps).all()
    assert (abs(flux1 - flux5) < eps).all()

def test_flux_rot():
    def testit(w, n):
        n = n/norm(n)
        alpha = atan2(n[1], n[0])
        f = array([0, calc_p(w), 0, 0])
        w = tangent_w(w, n)
        flux2 = dot(T_rot(alpha), flux(w, n))
        assert (abs(flux2 - f) < eps).all()

    w = array([1.8, 1.3, 1.2, 800.])
    n = array([1, -1])
    testit(w, n)

    w = array([1.8, 1.3, 1.2, 800.])
    n = array([1, 1])
    testit(w, n)

    w = array([1.8, 1.4, 1.2, 800.])
    n = array([1, 1])
    testit(w, n)

def test_numerical_flux():
    def testit(w_l, w_r):
        w_l = array(w_l, dtype="double")
        w_r = array(w_r, dtype="double")
        def testn(n):
            # consistency:
            f1 = numerical_flux(w_l, w_l, n)
            f2 = flux(w_l, n)
            assert (abs(f1 - f2) < eps).all()

            # conservativity:
            f1 = numerical_flux(w_l, w_r, n)
            f2 = -numerical_flux(w_r, w_l, -n)
            assert (abs(f1 - f2) < eps).all()

            f1 = numerical_flux(w_l, w_r, n)
            f2 = numerical_flux(w_r, w_l, -n)
            assert not (abs(f1 - f2) < eps).all()
        test_normals = [
            array([1, 0]),
            array([0, 1]),
            array([-1, 0]),
            array([0, -1]),
            array([1, 1]),
            array([1, -1]),
            array([-1, 1]),
            array([-1, -1]),
            ]
        for n in test_normals:
            n = n/float(norm(n))
            testn(n)


    test_states = [
        array([1, 1, 0, 10]),
        array([1, 0, 1, 10]),
        array([1, 1, 1, 10]),
        array([1, -1, 0, 10]),
        array([1, 0, -1, 10]),
        array([1, -1, -1, 10]),
        array([1, 1, -1, 10]),
        array([1, -1, 1, 10]),

        array([1.5, 1.1, 0, 10.3]),
        array([1.5, 0, 1.2, 10.3]),
        array([1.5, 1.1, 1.2, 10.3]),
        array([1.5, -1.1, 0, 10.3]),
        array([1.5, 0, -1.2, 10.3]),
        array([1.5, -1.1, -1.2, 10.3]),
        array([1.5, 1.1, -1.2, 10.3]),
        array([1.5, -1.1, 1.2, 10.3]),

        array([1.1, -10, 13, 700.1]),
        array([1.1, -10, 13, 800.1]),
        ]
    for x in test_states:
        for y in test_states:
            testit(x, y)
