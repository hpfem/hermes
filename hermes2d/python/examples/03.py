#! /usr/bin/env python

# This example shows how to solve a first simple PDE:
#   - load the mesh,
#   - perform initial refinements
#   - create a H1 space over the mesh
#   - define weak formulation
#   - initialize matrix solver
#   - assemble and solve the matrix system
#   - visualize the solution
#
# PDE: Poisson equation -Laplace u = CONST_F with homogeneous (zero)
#      Dirichlet boundary conditions.
#
# You can change the constant right-hand side CONST_F, the
# initial polynomial degree P_INIT, and play with various initial
# mesh refinements at the beginning.

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        WeakForm, Solution, ScalarView, LinSystem, DummySolver
from hermes2d.forms import set_forms
from hermes2d.examples.c03 import set_bc
from hermes2d.examples import get_example_mesh

P_INIT = 5                # Uniform polynomial degree of mesh elements.

# Problem parameters.
CONST_F = 2.0

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())

# Sample "manual" mesh refinement
mesh.refine_all_elements()

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm(1)
set_forms(wf)

# Initialize the linear system
ls = LinSystem(wf)
ls.set_spaces(space)

# Assemble and solve the matrix problem.
sln = Solution()
ls.assemble()
ls.solve_system(sln)

# Visualize the solution
sln.plot()

# Visualize the mesh
mesh.plot(space=space)
