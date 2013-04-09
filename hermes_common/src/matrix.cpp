// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file matrix.cpp
\brief Basic matrix classes and operations.
*/
#include "common.h"
#include "matrix.h"
#include "callstack.h"

#include "solvers/linear_matrix_solver.h"
#include "solvers/umfpack_solver.h"
#include "solvers/superlu_solver.h"
#include "solvers/amesos_solver.h"
#include "solvers/petsc_solver.h"
#include "solvers/mumps_solver.h"
#include "solvers/aztecoo_solver.h"
#include "qsort.h"
#include "api.h"

void Hermes::Algebra::DenseMatrixOperations::ludcmp(double **a, int n, int *indx, double *d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double *vv = new double[n];

  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big = 0.0;
    for (j = 0; j < n; j++)
      if((temp = fabs(a[i][j])) > big)
        big = temp;
    if(big == 0.0)
    {
      delete [] vv;
      throw Exceptions::Exception("Singular matrix in routine LUDCMP!");
    }
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if((dum = vv[i]*fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if(j != imax)
    {
      for (k = 0; k < n; k++)
      {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0.0) a[j][j] = 1.0e-20;
    if(j != n-1)
    {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
}

void Hermes::Algebra::DenseMatrixOperations::choldc(double **a, int n, double p[])
{
  int i, j, k;
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      double sum = a[i][j];
      k = i;
      while (--k >= 0) sum -= a[i][k] * a[j][k];
      if(i == j)
      {
        if(sum <= 0.0)
          throw Exceptions::Exception("CHOLDC failed!");
        else p[i] = sqrt(sum);
      }
      else a[j][i] = sum / p[i];
    }
  }
}

template<typename Scalar>
void Hermes::Algebra::Matrix<Scalar>::set_row_zero(unsigned int n)
{
  throw Hermes::Exceptions::MethodNotOverridenException("Matrix<Scalar>::set");
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::SparseMatrix()
{
  this->size = 0;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
{
  this->size = size;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::~SparseMatrix()
{
  if(pages)
  {
    for (unsigned int i = 0; i < this->size; i++)
      if(pages[i])
        delete pages[i];
    delete [] pages;
  }
}

template<typename Scalar>
void Hermes::Algebra::SparseMatrix<Scalar>::prealloc(unsigned int n)
{
  this->size = n;

  pages = new Page *[n];
  memset(pages, 0, n * sizeof(Page *));
}

template<typename Scalar>
void Hermes::Algebra::SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
{
  if(pages[col] == NULL || pages[col]->count >= PAGE_SIZE)
  {
    Page *new_page = new Page;
    new_page->count = 0;
    new_page->next = pages[col];
    pages[col] = new_page;
  }
  pages[col]->idx[pages[col]->count++] = row;
}

template<typename Scalar>
int Hermes::Algebra::SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
{
  // gather all pages in the buffer, deleting them along the way
  int *end = buffer;
  while (page != NULL)
  {
    memcpy(end, page->idx, sizeof(int) * page->count);
    end += page->count;
    Page *tmp = page;
    page = page->next;
    delete tmp;
  }

  // sort the indices and remove duplicities
  qsort_int(buffer, end - buffer);
  int *q = buffer;
  for (int *p = buffer, last = -1; p < end; p++) if(*p != last) *q++= last = *p;

  return q - buffer;
}

template<typename Scalar>
int Hermes::Algebra::SparseMatrix<Scalar>::get_num_indices()
{
  int total = 0;
  for (unsigned int i = 0; i < this->size; i++)
    for (Page *page = pages[i]; page != NULL; page = page->next)
      total += page->count;

  return total;
}

template<typename Scalar>
SparseMatrix<Scalar>* Hermes::Algebra::create_matrix()
{
  switch (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
  {
  case Hermes::SOLVER_AMESOS:
    {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
      return new EpetraMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_AZTECOO:
    {
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
      return new EpetraMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_MUMPS:
    {
#ifdef WITH_MUMPS
      return new MumpsMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_PETSC:
    {
#ifdef WITH_PETSC
      return new PetscMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_UMFPACK:
    {
#ifdef WITH_UMFPACK
      return new UMFPackMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_SUPERLU:
    {
#ifdef WITH_SUPERLU
      return new SuperLUMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
      break;
    }
  default:
    throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
  }
  return NULL;
}

template<typename Scalar>
Vector<Scalar>* Hermes::Algebra::create_vector()
{
  switch (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
  {
  case Hermes::SOLVER_AMESOS:
    {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
      return new EpetraVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_AZTECOO:
    {
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
      return new EpetraVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_MUMPS:
    {
#ifdef WITH_MUMPS
      return new MumpsVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_PETSC:
    {
#ifdef WITH_PETSC
      return new PetscVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_UMFPACK:
    {
#ifdef WITH_UMFPACK
      return new UMFPackVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_SUPERLU:
    {
#ifdef WITH_SUPERLU
      return new SuperLUVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
      break;
    }
  default:
    throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
  }
  return NULL;
}

template class Hermes::Algebra::SparseMatrix<double>;
template class Hermes::Algebra::SparseMatrix<std::complex<double> >;

template HERMES_API Vector<double>* Hermes::Algebra::create_vector();
template HERMES_API SparseMatrix<double>*  Hermes::Algebra::create_matrix();

template HERMES_API Vector<std::complex<double> >* Hermes::Algebra::create_vector();
template HERMES_API SparseMatrix<std::complex<double> >*  Hermes::Algebra::create_matrix();