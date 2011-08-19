#! /usr/bin/env python

# This example solves a general second-order linear equation with non-constant
# coefficients, and shows how integration orders in linear and bilinear forms
# can be defined manually.

# PDE: -d/dx(a_11(x,y)du/dx) - d/dx(a_12(x,y)du/dy) - d/dy(a_21(x,y)du/dx) - d/dy(a_22(x,y)du/dy)
#      + a_1(x,y)du/dx + a_21(x,y)du/dy + a_0(x,y)u = rhs(x,y)
#
# Domain: arbitrary
#
# BC:  Dirichlet for boundary marker 1: u = g_D(x,y)
#      Natural for any other boundary marker:   (a_11(x,y)*nu_1 + a_21(x,y)*nu_2) * dudx
#                                             + (a_12(x,y)*nu_1 + s_22(x,y)*nu_2) * dudy = g_N(x,y)

# Import modules
from hermes2d import Mesh, MeshView, OrderView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView, VonMisesFilter, \
        OrderView

from hermes2d.examples.c07 import set_bc, set_forms
from hermes2d.examples import get_07_mesh

# The following parameters can be changed:
P_INIT = 2             # Initial polynomial degree of all mesh elements.
INIT_REF_NUM = 4       # Number of initial uniform refinements

# Load the mesh
mesh = Mesh()
mesh.load(get_07_mesh())

# Perform initial mesh refinements.
for i in range(INIT_REF_NUM):
    mesh.refine_all_elements()

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf)

# Initialize views
sview = ScalarView("Coarse solution", 0, 100, 798, 700)
oview = OrderView("Polynomial orders", 800, 100, 798, 700)

# Initialize the linear system.
ls = LinSystem(wf)
ls.set_spaces(space)

# Assemble and solve the matrix problem
sln = Solution()
ls.assemble()
ls.solve_system(sln)

# View the solution
sln.plot()

# View the mesh
mesh.plot(space=space)
