#! /usr/bin/env python

# This example illustrates how to use non-homogeneous(nonzero)
# Dirichlet boundary conditions.
#
# PDE: Poisson equation -Laplace u = CONST_F, where CONST_F is
# a constant right-hand side. It is not difficult to see that
# the function u(x,y) = (-CONST_F/4)*(x^2 + y^2) satisfies the
# above PDE. Since also the Dirichlet boundary conditions
# are chosen to match u(x,y), this function is the exact
# solution.
#
# Note that since the exact solution is a quadratic polynomial,
# Hermes will compute it exactly if all mesh elements are quadratic
# or higher (then the exact solution lies in the finite element space).
# If some elements in the mesh are linear, Hermes will only find
# an approximation.

# Import modules
from hermes2d import (Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space,
        LinSystem, Solution, ScalarView, WeakForm, DummySolver)

from hermes2d.examples.c04 import set_bc
from hermes2d.examples import get_example_mesh
from hermes2d.forms import set_forms

# Below you can play with the parameters CONST_F, P_INIT, and UNIFORM_REF_LEVEL.
INIT_REF_NUM = 2         # number of initial uniform mesh refinements
P_INIT = 2               # initial polynomial degree in all elements

# Load the mesh file
mesh = Mesh()
mesh.load(get_example_mesh())

# Perform initial mesh refinements
for i in range(INIT_REF_NUM):
    mesh.refine_all_elements()

# Create an H1 space with default shapeset
space = H1Space(mesh, P_INIT)
set_bc(space)

# Initialize the weak formulation
wf = WeakForm()
set_forms(wf)

# Initialize the linear system
ls = LinSystem(wf)
ls.set_spaces(space)

# Assemble and solve the matrix problem
sln = Solution()
ls.assemble()
ls.solve_system(sln)

# Visualize the solution
sln.plot()

# Visualize the mesh
mesh.plot(space=space)
