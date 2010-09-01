#! /usr/bin/env python

# This example explains how to create two spaces over a mesh and use them
# to solve a simple problem of linear elasticity. At the end, VonMises
# filter is used to visualize the stress.
#
# PDE: Lame equations of linear elasticity
#
# BC: du_1/dn = f_0 on Gamma_3 and du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5
#     du_2/dn = f_1 on Gamma_3 and du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5
#     u_1 = 0 and u_2 = 0 on Gamma_1

# Import modules
from hermes2d import Mesh, MeshView, H1Shapeset, PrecalcShapeset, H1Space, \
        LinSystem, WeakForm, DummySolver, Solution, ScalarView, VonMisesFilter

from hermes2d.examples.c08 import set_bc, set_forms
from hermes2d.examples import get_sample_mesh

# The following parameter can be changed:
P_INIT = 4

# Load the mesh file
mesh = Mesh()
mesh.load(get_sample_mesh())

# Perform uniform mesh refinement
mesh.refine_all_elements()

# Create the x- and y- displacement space using the default H1 shapeset
xdisp = H1Space(mesh, P_INIT)
ydisp = H1Space(mesh, P_INIT)
set_bc(xdisp, ydisp)

# Initialize the weak formulation
wf = WeakForm(2)
set_forms(wf)

# Initialize the linear system.
ls = LinSystem(wf)
ls.set_spaces(xdisp, ydisp)

# Assemble and solve the matrix problem
xsln = Solution()
ysln = Solution()
ls.assemble()
ls.solve_system(xsln, ysln, lib="scipy")

# Visualize the solution
view = ScalarView("Von Mises stress [Pa]", 50, 50, 1200, 600)
E = float(200e9)
nu = 0.3
l = (E * nu) / ((1 + nu) * (1 - 2*nu))
mu = E / (2*(1 + nu))
stress = VonMisesFilter(xsln, ysln, mu, l)
view.show(stress)

# Visualize the mesh
mesh.plot(space=xdisp)
