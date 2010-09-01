// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

// Set component "comp" of the solution to be a constant "val" everywhere
// Note: This function does not touch Dirichlet boundary 
// conditions, those must be set to "val" separately.
void set_vertex_dofs_constant(Mesh* mesh, double val, int comp=0, int sln=0);

// Multiply (all components) of the solution at all points by 'val'.
// Caution: This does not work when Dirichlet conditions 
// are present - the lifts must be multiplied separately.
void multiply_dofs_with_constant(Mesh* mesh, double val, int sln=0);

// Copies all solution coefficients for component "comp" from 
// solution "sln_src" to target solution "sln_trg"
void copy_dofs(int sln_src, int sln_trg, Mesh* mesh, int comp=0);
