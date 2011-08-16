#! /usr/bin/env python

# This example shows how to define Neumann boundary conditions. In addition,
# you will see how a Filter is used to visualize gradient of the solution
#
# PDE: Poisson equation -Laplace u = f, where f = CONST_F
#
# BC: u = 0 on Gamma_4 (edges meeting at the re-entrant corner)
#     du/dn = CONST_GAMMA_1 on Gamma_1 (bottom edge)
#     du/dn = CONST_GAMMA_2 on Gamma_2 (top edge, circular arc, and right-most edge)
#     du/dn = CONST_GAMMA_3 on Gamma_3 (left-most edge)
#
# You can play with the parameters below. For most choices of the four constants,
# the solution has a singular (infinite) gradient at the re-entrant corner.
# Therefore we visualize not only the solution but also its gradient.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, Solution, ScalarView, WeakForm, DummySolver

from hermes2d.examples.c05 import set_bc, set_forms
from hermes2d.examples import get_example_mesh

P_INIT = 4                           # initial polynomial degree in all elements
CORNER_REF_LEVEL = 12                # number of mesh refinements towards the re-entrant corner

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())

# Perform initial mesh refinements.
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
