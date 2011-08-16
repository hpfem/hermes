#! /usr/bin/env python

"""
Displays a mesh file.

Usage:

./show_mesh.py file.mesh
"""

from hermes2d import Mesh, MeshView
import sys

filename = sys.argv[1]

mesh = Mesh()
mesh.load(filename)

mview = MeshView(filename, 100, 100, 500, 500)
mview.show(mesh, lib="mpl", method="orders")
