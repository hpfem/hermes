// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "discrete.h"

// Set component "comp" of the solution to be a constant "val" everywhere
// Note: This function does not touch Dirichlet boundary 
// conditions, those must be set to "val" separately.
void set_vertex_dofs_constant(Mesh* mesh, double val, int comp, int sln)
{
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    e->coeffs[sln][comp][0] = val;
    e->coeffs[sln][comp][1] = val;
  }
  delete I;
}

// Multiply (all components) of the solution at all points by 'val'.
// Caution: This does not work when Dirichlet conditions 
// are present - the lifts must be multiplied separately.
void multiply_dofs_with_constant(Mesh* mesh, double val, int sln)
{
  int n_dof = mesh->get_n_dof();
  double *y = new double[n_dof];
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) e->copy_coeffs_to_vector(y, sln);
  for (int i = 0; i < n_dof; i++) y[i] *= val;
  I->reset();	
  while ((e = I->next_active_element()) != NULL)
    e->get_coeffs_from_vector(y, sln);  
  delete I; 
  delete [] y;
}

// Copies all solution coefficients for component "comp" from 
// solution "sln_src" to target solution "sln_trg"
void copy_dofs(int sln_src, int sln_trg, Mesh* mesh, int comp=0) 
{
  if(sln_src < 0 || sln_src > MAX_SLN_NUM) error("wrong solution index in copy_dofs().");
  if(sln_trg < 0 || sln_trg > MAX_SLN_NUM) error("wrong solution index in copy_dofs().");
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    e->copy_dofs(sln_src, sln_trg);
  }
}
