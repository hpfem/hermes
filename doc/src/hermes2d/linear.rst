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
time-dependent problems. 

.. toctree::
   :maxdepth: 2

   linear/mesh   
   linear/space
   linear/poisson
   linear/bc
   linear/quadrature
   linear/general
   linear/system
   linear/filters








