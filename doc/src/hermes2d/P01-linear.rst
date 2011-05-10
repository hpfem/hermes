Tutorial Part I (Linear Problems)
=================================

The first tutorial chapter begins with explaining the mesh
data format and showing how to load, refine, and visualize meshes. 
Then you will learn how to set up a finite element space and 
solve a first simple problem - a Poisson equation with zero Dirichlet 
boundary conditions. After that we show how to prescribe more
general boundary conditions. We also explain how Hermes
handles numerical quadrature since this is both very important 
for higher-order finite element methods and very different 
from standard low-order FEM codes. At the end of this 
chapter you will learn how to solve systems of equations and 
axisymmetric 3D problems. 

.. toctree::
   :maxdepth: 2

   P01-linear/01-mesh   
   P01-linear/02-space
   P01-linear/03-poisson
   P01-linear/essential_and_natural_bc
   P01-linear/04-bc-dirichlet
   P01-linear/05-bc-neumann
   P01-linear/06-bc-newton
   P01-linear/quadrature
   P01-linear/07-general
   P01-linear/08-system
   P01-linear/filters
   P01-linear/09-axisym








