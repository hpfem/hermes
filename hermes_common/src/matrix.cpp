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
#include "solvers/paralution_solver.h"
#include "qsort.h"
#include "api.h"

namespace Hermes
{
  namespace Algebra
  {
    void DenseMatrixOperations::ludcmp(double **a, int n, int *indx, double *d)
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

    void DenseMatrixOperations::choldc(double **a, int n, double p[])
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
    void Matrix<Scalar>::set_row_zero(unsigned int n)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("Matrix<Scalar>::set");
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix()
    {
      this->size = 0;
      pages = NULL;

      row_storage = false;
      col_storage = false;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
    {
      this->size = size;
      pages = NULL;

      row_storage = false;
      col_storage = false;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::~SparseMatrix()
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
    void SparseMatrix<Scalar>::prealloc(unsigned int n)
    {
      this->size = n;

      pages = new Page *[n];
      memset(pages, 0, n * sizeof(Page *));
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
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
    int SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
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
      for (int *p = buffer, last = -1; p < end; p++)
        if(*p != last)
          *q++ = last = *p;

      return q - buffer;
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_indices()
    {
      int total = 0;
      for (unsigned int i = 0; i < this->size; i++)
        for (Page *page = pages[i]; page != NULL; page = page->next)
          total += page->count;

      return total;
    }


    template<typename Scalar>
    Vector<Scalar>::Vector() : size(0)
    {
    }

    template<typename Scalar>
    Vector<Scalar>::Vector(unsigned int size) : size(size)
    {
    }

    template<typename Scalar>
    void Vector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->set(i, vec->get(i));
    }

    template<typename Scalar>
    void Vector<Scalar>::set_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->length(); i++) this->set(i, vec[i]);
    }

    template<typename Scalar>
    void Vector<Scalar>::add_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec->get(i));
    }

    template<typename Scalar>
    void Vector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec[i]);
    }

    template HERMES_API void Vector<double>::set_vector(Vector<double>* vec);
    template HERMES_API void Vector<double>::set_vector(double* vec);
    template HERMES_API void Vector<double>::add_vector(Vector<double>* vec);
    template HERMES_API void Vector<double>::add_vector(double* vec);

    template HERMES_API void Vector<std::complex<double> >::set_vector(Vector<std::complex<double> >* vec);
    template HERMES_API void Vector<std::complex<double> >::set_vector(std::complex<double> * vec);
    template HERMES_API void Vector<std::complex<double> >::add_vector(Vector<std::complex<double> >* vec);
    template HERMES_API void Vector<std::complex<double> >::add_vector(std::complex<double> * vec);


    template<>
    HERMES_API SparseMatrix<double>* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new UMFPackMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SuperLUMatrix<double>;
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

    template<>
    HERMES_API Vector<double>* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsVector<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new UMFPackVector<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionVector<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SuperLUVector<double>;
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

    template<>
    HERMES_API SparseMatrix<std::complex<double> >* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new UMFPackMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SuperLUMatrix<std::complex<double> >;
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

    template<>
    HERMES_API Vector<std::complex<double> >* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");

#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new UMFPackVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SuperLUVector<std::complex<double> >;
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

    template class SparseMatrix<double>;
    template class SparseMatrix<std::complex<double> >;

    template class Vector<double>;
    template class Vector<std::complex<double> >;
  }
}