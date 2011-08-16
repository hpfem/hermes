#! /usr/bin/env python

# This example explains how to use Newton boundary conditions. Again,
# a Filter is used to visualize the solution gradient.
#
# PDE: Laplace equation -Laplace u = 0 (this equation describes, among
# many other things, also stationary heat transfer in a homogeneous linear
# material).
#
# BC: u = T1 ... fixed temperature on Gamma_3 (Dirichlet)
#     du/dn = 0 ... insulated wall on Gamma_2 and Gamma_4 (Neumann)
#     du/dn = H*(u - T0) ... heat flux on Gamma_1 (Newton)
#
# Note that the last BC can be written in the form  du/dn - H*u = -h*T0.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView

from hermes2d.examples.c06 import set_bc, set_forms
from hermes2d.examples import get_example_mesh

# The following parameters can be changed:

UNIFORM_REF_LEVEL = 2;   # Number of initial uniform mesh refinements.
CORNER_REF_LEVEL = 12;   # Number of mesh refinements towards the re-entrant corner.
P_INIT = 6;              # Uniform polynomial degree of all mesh elements.

# Boundary markers
NEWTON_BDY = 1

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())

# Perform initial mesh refinements.
for i in range(UNIFORM_REF_LEVEL):
    mesh.refine_all_elements()
mesh.refine_towards_vertex(3, CORNER_REF_LEVEL)

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf)

# Initialize the linear system.
ls = LinSystem(wf)
ls.set_spaces(space)

# Assemble and solve the matrix problem
sln = Solution()
ls.assemble()
ls.solve_system(sln)

# Visualize the approximation
sln.plot()

# Visualize the mesh
mesh.plot(space=space)
