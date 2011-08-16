#! /usr/bin/env python

# This example shows how to load a mesh, perform various types
# of "manual"  element refinements.

# Import modules
from hermes2d import Mesh, MeshView
from hermes2d.examples import get_example_mesh

# Load the mesh file
domain_mesh = get_example_mesh()
mesh = Mesh()
mesh.load(domain_mesh)

# Perform some sample initial refinements
mesh.refine_all_elements();           # Refines all elements.
#mesh.refine_towards_vertex(3, 4);    # Refines mesh towards vertex #3 (4x).
#mesh.refine_towards_boundary(2, 4);  # Refines all elements along boundary 2 (4x).

# Display the mesh
mesh.plot()
