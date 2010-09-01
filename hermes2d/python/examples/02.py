#! /usr/bin/env python

# This example demonstrates how to create a space over a mesh
# and how to visualize the finite element basis functions
# using the BaseView class. The variable P_INIT defines the
# initial degree of all mesh elements. Currently, it is set
# to 1. After visualizing the basis functions using the hints
# printed in the terminal window, change the value of P_INIT
# to 2, 3, etc. to also see higher-order shape functions.
#
# You can use this example to visualize all shape functions
# on the reference square and reference triangle domains,
# just load the corresponding mesh at the beginning of the file.

# Import modules
from hermes2d import Mesh, H1Shapeset, PrecalcShapeset, H1Space, \
        BaseView

from hermes2d.forms import set_forms
from hermes2d.examples import get_example_mesh

P_INIT = 3

# Load the mesh file
domain_mesh = get_example_mesh()    # Original L-shape domain
mesh = Mesh()
mesh.load(domain_mesh)

# Refine all elements (optional)
mesh.refine_all_elements()

# Create a shapeset and an H1 space
space = H1Space(mesh)

# Assign element orders and initialize the space
space.set_uniform_order(P_INIT)    # Set uniform polynomial order
                                   # P_INIT to all mesh elements.

# View the basis functions
bview = BaseView()
bview.show(space)
bview.wait()
