// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "common.h"
#include "matrix.h"
#include "error.h"
#include "callstack.h"

//  Solvers
#include "solver/solver.h"
#include "solver/umfpack_solver.h"
#include "solver/superlu.h"
#include "solver/amesos.h"
#include "solver/petsc.h"
#include "solver/mumps.h"
#include "solver/nox.h"
#include "solver/aztecoo.h"

#define HERMES_TINY 1.0e-20


// ludcmp, lubksb - LU decomposition and back-substitution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void ludcmp(double **a, int n, int *indx, double *d)
{
  _F_
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double *vv = new double[n];
  MEM_CHECK(vv);

  *d = 1.0;
  for (i = 0; i < n; i++) {
    big=0.0;
    for (j = 0; j < n; j++) if ((temp = fabs(a[i][j])) > big) big = temp;
    if (big == 0.0) EXIT("Singular matrix in routine LUDCMP!");
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i]*fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = HERMES_TINY;
    if (j != n-1) {
      dum = 1.0 / (a[j][j]);
      for (i = j+1; i < n; i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
}

// choldc, cholsl - Cholesky decomposition and solution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void choldc(double **a, int n, double p[])
{
  _F_
  int i, j, k;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      double sum = a[i][j];
      k = i;
      while (--k >= 0) sum -= a[i][k] * a[j][k];
      if (i == j) {
	if (sum <= 0.0) EXIT("CHOLDC failed!");
	else p[i] = sqrt(sum);
      }
      else a[j][i] = sum / p[i];
    }
  }
}

// Simple dot product.
double vec_dot(double *r, double *s, int n_dof)
{
  double result = 0;
  for (int i=0; i < n_dof; i++) result += r[i]*s[i];
  return result;
}

template<typename Scalar>
double vec_dot(Vector<Scalar> *r, Vector<Scalar> *s, int n_dof)
{
  double result = 0;
#ifndef HERMES_COMMON_COMPLEX
  for (int i=0; i < n_dof; i++) result += r->get(i)*s->get(i);
#endif
  return result;
}

// SparseMatrix<Scalar> ////////////////////////////////////////////////////////////////////////////////////

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

template class SparseMatrix<double>;
template class SparseMatrix<std::complex<double>>;

template<typename Scalar>
SparseMatrix<Scalar>::SparseMatrix()
{
  _F_
  size = 0;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
{
  _F_
  this->size = size;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
SparseMatrix<Scalar>::~SparseMatrix()
{
  _F_
  delete [] pages;
}

template<typename Scalar>
void SparseMatrix<Scalar>::prealloc(unsigned int n)
{
  _F_
  this->size = n;

  pages = new Page *[n];
  MEM_CHECK(pages);
  memset(pages, 0, n * sizeof(Page *));
}

template<typename Scalar>
void SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
{
  _F_
  if (pages[col] == NULL || pages[col]->count >= PAGE_SIZE) {
    Page *new_page = new Page;
    MEM_CHECK(new_page);
    new_page->count = 0;
    new_page->next = pages[col];
    pages[col] = new_page;
  }
  pages[col]->idx[pages[col]->count++] = row;
}

template<typename Scalar>
int SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
{
  _F_
  // gather all pages in the buffer, deleting them along the way
  int *end = buffer;
  while (page != NULL) {
    memcpy(end, page->idx, sizeof(int) * page->count);
    end += page->count;
    Page *tmp = page;
    page = page->next;
    delete tmp;
  }

  // sort the indices and remove duplicities
  qsort_int(buffer, end - buffer);
  int *q = buffer;
  for (int *p = buffer, last = -1; p < end; p++) if (*p != last) *q++ = last = *p;

  return q - buffer;
}

template<typename Scalar>
int SparseMatrix<Scalar>::get_num_indices()
{
  _F_
  int total = 0;
  for (unsigned int i = 0; i < size; i++)
    for (Page *page = pages[i]; page != NULL; page = page->next)
      total += page->count;

  return total;
}

template<typename Scalar>
SparseMatrix<Scalar>* create_matrix(MatrixSolverType matrix_solver)
{
  _F_
  switch (matrix_solver) 
  {
    case SOLVER_AMESOS:
    case SOLVER_AZTECOO:
      {
        return new EpetraMatrix<Scalar>;
        break;
      }
    case SOLVER_MUMPS: 
      {
        return new MumpsMatrix<Scalar>;
        break;
      }
    case SOLVER_PETSC: 
      {
        return new PetscMatrix<Scalar>;
        break;
      }
    case SOLVER_UMFPACK: 
      {
        return new UMFPackMatrix<Scalar>;
        break;
      }
    case SOLVER_SUPERLU: 
    {
      return new SuperLUMatrix<Scalar>;
      break;
    }
    default: 
      error("Unknown matrix solver requested.");
  }
  return NULL;
}

template<typename Scalar>
Solver<Scalar>* create_linear_solver(MatrixSolverType matrix_solver, Matrix<Scalar>* matrix, Vector<Scalar>* rhs)
{
  _F_
  Vector<Scalar>* rhs_dummy = NULL;
  switch (matrix_solver) 
  {
    case SOLVER_AZTECOO:
      {
        info("Using AztecOO."); 
        if (rhs != NULL) return new AztecOOSolver<Scalar>(static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs));
        else return new AztecOOSolver<Scalar>(static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs_dummy));
        break;
      }
    case SOLVER_AMESOS:
      {
        info("Using Amesos.");         
        if (rhs != NULL) return new AmesosSolver<Scalar>("Amesos_Klu", static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs));
        else return new AmesosSolver<Scalar>("Amesos_Klu", static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs_dummy));
        break;
      }
    case SOLVER_MUMPS: 
      {
        info("Using Mumps.");         
        if (rhs != NULL) return new MumpsSolver<Scalar>(static_cast<MumpsMatrix<Scalar>*>(matrix), static_cast<MumpsVector<Scalar>*>(rhs)); 
        else return new MumpsSolver<Scalar>(static_cast<MumpsMatrix<Scalar>*>(matrix), static_cast<MumpsVector<Scalar>*>(rhs_dummy)); 
        break;
      }
    case SOLVER_PETSC: 
      {
        info("Using PETSc.");        
        if (rhs != NULL) return new PetscLinearSolver<Scalar>(static_cast<PetscMatrix<Scalar>*>(matrix), static_cast<PetscVector<Scalar>*>(rhs)); 
        else return new PetscLinearSolver<Scalar>(static_cast<PetscMatrix<Scalar>*>(matrix), static_cast<PetscVector<Scalar>*>(rhs_dummy)); 
        break;
      }
    case SOLVER_UMFPACK: 
      {
        info("Using UMFPack.");
        if (rhs != NULL) return new UMFPackLinearSolver<Scalar>(static_cast<UMFPackMatrix<Scalar>*>(matrix), static_cast<UMFPackVector<Scalar>*>(rhs)); 
        else return new UMFPackLinearSolver<Scalar>(static_cast<UMFPackMatrix<Scalar>*>(matrix), static_cast<UMFPackVector<Scalar>*>(rhs_dummy));  
        break;
      }
    case SOLVER_SUPERLU: 
    {
      info("Using SuperLU.");       
      if (rhs != NULL) return new SuperLUSolver<Scalar>(static_cast<SuperLUMatrix<Scalar>*>(matrix), static_cast<SuperLUVector<Scalar>*>(rhs)); 
      else return new SuperLUSolver<Scalar>(static_cast<SuperLUMatrix<Scalar>*>(matrix), static_cast<SuperLUVector<Scalar>*>(rhs_dummy)); 
      break;
    }
    default: 
      error("Unknown matrix solver requested.");
  }
  return NULL;
}

template<typename Scalar>
Vector<Scalar>* create_vector(MatrixSolverType matrix_solver)
{
  _F_
  switch (matrix_solver) 
  {
    case SOLVER_AMESOS:
    case SOLVER_AZTECOO:
      {
        return new EpetraVector<Scalar>;
        break;
      }
    case SOLVER_MUMPS: 
      {
        return new MumpsVector<Scalar>;
        break;
      }
    case SOLVER_PETSC: 
      {
        return new PetscVector<Scalar>;
        break;
      }
    case SOLVER_UMFPACK: 
      {
        return new UMFPackVector<Scalar>;
        break;
      }
    case SOLVER_SUPERLU: 
    {
      return new SuperLUVector<Scalar>;
      break;
    }
    default: 
      error("Unknown matrix solver requested.");
  }
  return NULL;
}