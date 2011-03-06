from math import exp, e
from numpy import pi, arange, array, sin, cos
from numpy.linalg import solve

from hermes1d.h1d_wrapper.h1d_wrapper import (assemble_projection_matrix_rhs,
        Mesh, FESolution, Function as Function2)
from hermes1d.fekete.fekete import Function, Mesh1D
from hermes_common.matrix import CSCMatrix, AVector

class FunctionSin(Function2):

    def eval_f(self, x):
        return sin(x)

    def eval_dfdx(self, x):
        return cos(x)

f_sin = FunctionSin()


def test_l2_h1_proj_run():
    """
    Test that the projections run.

    It doesn't test if it's correct.
    """
    pts = arange(0, 2*pi, 1)
    orders = [3]*(len(pts)-1)
    m = Mesh(pts, orders)
    n_dof = m.assign_dofs()
    A = CSCMatrix(n_dof)
    rhs = AVector(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CSCMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_h1 = FESolution(m, x).to_discrete_function()
    sol_l2.plot(False)
    sol_h1.plot(False)

def test_l2_h1_proj1():
    """
    Tests the correctness of the projections.
    """
    pts = arange(0, 2*pi, 3)
    orders = [2]*(len(pts)-1)
    m = Mesh(pts, orders)

    pts = array(list(arange(0, pts[-1], 0.1)) + [pts[-1]])
    orders = [6]*(len(pts)-1)
    f_exact = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CSCMatrix(n_dof)
    rhs = AVector(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CSCMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert (sol_l2 - f_exact).l2_norm() < 0.07
    assert (sol_h1 - f_exact).l2_norm() < 0.07

def test_l2_h1_proj2():
    """
    Tests the correctness of the projections.
    """
    pts = arange(0, 2*pi, 3)
    orders = [4]*(len(pts)-1)
    m = Mesh(pts, orders)

    pts = array(list(arange(0, pts[-1], 0.1)) + [pts[-1]])
    orders = [6]*(len(pts)-1)
    f_exact = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CSCMatrix(n_dof)
    rhs = AVector(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="L2")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CSCMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f_sin, projection_type="H1")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert (sol_l2 - f_exact).l2_norm() < 0.002
    assert (sol_h1 - f_exact).l2_norm() < 0.002

def test_l2_h1_proj3():
    """
    Tests conversion to FE basis.
    """
    pts = arange(0, 2*pi, 0.1)
    orders = [2]*(len(pts)-1)
    m = Mesh(pts, orders)

    f = Function(lambda x: sin(x), Mesh1D(pts, orders))

    n_dof = m.assign_dofs()
    A = CSCMatrix(n_dof)
    rhs = AVector(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f, projection_type="L2")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_l2 = FESolution(m, x).to_discrete_function()
    A = CSCMatrix(n_dof)
    assemble_projection_matrix_rhs(m, A, rhs, f, projection_type="H1")
    x = solve(A.to_scipy_csc().todense(), rhs.to_numpy())
    sol_h1 = FESolution(m, x).to_discrete_function()
    assert sol_l2 == f
    assert sol_h1 == f

def test_l2_h1_proj4():
    """
    Tests conversion to FE basis.
    """
    pts = arange(0, 2*pi, 0.4)
    orders = [2]*(len(pts)-1)
    m = Mesh1D(pts, orders)

    f = Function(lambda x: sin(x), m)
    assert f.project_onto(m, proj_type="Fekete") == f
    assert f.project_onto(m, proj_type="L2") == f
    assert f.project_onto(m, proj_type="H1") == f

    orders = [3]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    assert f.project_onto(m, proj_type="Fekete") == f
    assert f.project_onto(m, proj_type="L2") == f
    assert f.project_onto(m, proj_type="H1") == f

    orders = [4]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    assert f.project_onto(m, proj_type="Fekete") == f
    assert f.project_onto(m, proj_type="L2") == f
    assert f.project_onto(m, proj_type="H1") == f

    pts = arange(0, 2*pi, 3)
    orders = [2]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    pts = array(list(arange(0, pts[-1], 0.1)) + [pts[-1]])
    orders = [6]*(len(pts)-1)
    f_exact = Function(lambda x: sin(x), Mesh1D(pts, orders))

    sol_l2 = f_exact.project_onto(m, proj_type="L2")
    sol_h1 = f_exact.project_onto(m, proj_type="H1")
    assert (sol_l2 - f_exact).l2_norm() < 0.07
    assert (sol_h1 - f_exact).l2_norm() < 0.07

def test_l2_h1_proj5():
    """
    Tests exact projections.

    The exact results were generated using:

    from sympy import sin, cos, integrate, var, pi, exp, log, E, S
    var("x")
    f = exp(x)
    S(1)/2 * integrate(f, (x, -1, 1)) + S(3)/2*x*integrate(x*f, (x, -1, 1))

    """
    f_exact = lambda x: exp(x)
    f_exact_l2 = lambda x: e/2 - exp(-1)/2 + 3*x*exp(-1)
    # TODO: The constant term here has to be checked:
    f_exact_h1 = lambda x: +3*e/8 + 3*exp(-1)/8 + 3*x*(e + exp(-1))/8
    pts = [-1, -0.5, 0, 0.5, 1]
    orders = [20]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    f = Function(f_exact, m)

    pts = [-1, 1]
    orders = [1]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    f_proj_l2_exact = Function(f_exact_l2, m)
    f_proj_l2 = f.project_onto(m, proj_type="L2")
    f_proj_h1_exact = Function(f_exact_h1, m)
    f_proj_h1 = f.project_onto(m, proj_type="H1")
    eps_l2 = 1e-3
    eps_h1 = 0.03
    assert (f_proj_l2 - f_proj_l2_exact).l2_norm() < eps_l2
    assert (f_proj_h1 - f_proj_h1_exact).l2_norm() < eps_h1

    # Make sure that if we exchange the L2 and H1 solutions, then the test
    # fails:
    assert (f_proj_l2 - f_proj_h1_exact).l2_norm() > max(eps_l2, eps_h1)
    assert (f_proj_h1 - f_proj_l2_exact).l2_norm() > max(eps_l2, eps_h1)

def test_l2_h1_proj6():
    """
    Tests exact projections.

    Slightly more complicated example.
    """
    f_exact = lambda x: sin(x)*exp(x)
    f_exact_l2 = lambda x: -e*cos(1)/4 + e*sin(1)/4 + exp(-1)*sin(1)/4 + 3*x*(e*sin(1)/2 - exp(-1)*sin(1)/2 - cos(1)*exp(-1))/2 + cos(1)*exp(-1)/4
    pts = [-1, -0.5, 0, 0.5, 1]
    orders = [20]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    f = Function(f_exact, m)

    pts = [-1, 1]
    orders = [1]*(len(pts)-1)
    m = Mesh1D(pts, orders)
    f_proj_l2 = Function(f_exact_l2, m)
    assert (f.project_onto(m, proj_type="L2") - f_proj_l2).l2_norm() < 0.03
