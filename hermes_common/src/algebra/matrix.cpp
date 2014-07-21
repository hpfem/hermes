// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://www.hpfem.org/.
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
#include "util/memory_handling.h"

#include "solvers/linear_matrix_solver.h"
#include "solvers/interfaces/umfpack_solver.h"
#include "solvers/interfaces/superlu_solver.h"
#include "solvers/interfaces/amesos_solver.h"
#include "solvers/interfaces/petsc_solver.h"
#include "solvers/interfaces/mumps_solver.h"
#include "solvers/interfaces/aztecoo_solver.h"
#include "solvers/interfaces/paralution_solver.h"
#include "qsort.h"
#include "api.h"

namespace Hermes
{
  namespace Algebra
  {
    template<typename Scalar>
    Matrix<Scalar>::Matrix(unsigned int size) : size(size)
    {
    }

    template<typename Scalar>
    void Matrix<Scalar>::set_row_zero(unsigned int n)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("Matrix<Scalar>::set_row_zero");
    }

    template<typename Scalar>
    void Matrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i < this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void Matrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const
    {
      if (!vector_out_initialized)
        vector_out = malloc_with_check<Scalar>(this->size);
      for (unsigned int i = 0; i < this->size; i++)
      {
        vector_out[i] = Scalar(0.);
        for (unsigned int j = 0; j < this->size; j++)
          vector_out[i] += this->get(i, j) * vector_in[j];
      }
    }

    template<typename Scalar>
    void Matrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("multiply_with_Scalar()");
    }

    template<typename Scalar>
    unsigned int Matrix<Scalar>::get_size() const
    {
      return this->size;
    };

    template<>
    void Matrix<double>::add(unsigned int m, unsigned int n, double *mat, int *rows, int *cols, const int size)
    {
      for (unsigned int i = 0; i < m; i++)
      {
        for (unsigned int j = 0; j < n; j++)
        {
          double entry = mat[i * size + j];
          if (entry > HermesEpsilon || entry < -HermesEpsilon)
          {
            if (rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
              add(rows[i], cols[j], entry);
          }
        }
      }
    }

    template<>
    void Matrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double>* mat, int *rows, int *cols, const int size)
    {
      for (unsigned int i = 0; i < m; i++)
      {
        for (unsigned int j = 0; j < n; j++)
        {
          std::complex<double> entry = mat[i * size + j];
          if (entry.real() > HermesEpsilon || entry.real() < -HermesEpsilon || entry.imag() > HermesEpsilon || entry.imag() < -HermesEpsilon)
          {
            if (rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
              add(rows[i], cols[j], mat[i * size + j]);
          }
        }
      }
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix() : Matrix<Scalar>()
    {
      pages = nullptr;
      next_pages = nullptr;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
    {
      this->size = size;
      pages = nullptr;
      next_pages = nullptr;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::~SparseMatrix()
    {
      this->free();
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::free()
    {
      if (pages)
      {
        free_with_check(pages);
      }

      if (next_pages)
      {
        for (unsigned int i = 0; i < this->size; i++)
          delete next_pages[i];
        free_with_check(next_pages);
      }
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::finish()
    {
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_row_entries(unsigned int row) const
    {
      return -1;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::extract_row_copy(unsigned int row, unsigned int len,
      unsigned int &n_entries, double *vals,
      unsigned int *idxs) const
    {
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_col_entries(unsigned int col) const
    {
      return -1;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::extract_col_copy(unsigned int col, unsigned int len,
      unsigned int &n_entries, double *vals,
      unsigned int *idxs) const
    {
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_sparse_matrix(SparseMatrix<Scalar>* mat)
    {
      add_as_block(0, 0, mat);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      int ndof = mat->get_size();
      if (this->get_size() != (unsigned int)num_stages * ndof)
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in SparseMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
        this->add_as_block(ndof*i, ndof*i, mat);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, SparseMatrix<Scalar>* mat)
    {
      if ((this->get_size() < offset_i + mat->get_size()) || (this->get_size() < offset_j + mat->get_size()))
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in SparseMatrix<Scalar>::add_as_block()");
      unsigned int block_size = mat->get_size();
      for (unsigned int r = 0; r < block_size; r++)
      {
        for (unsigned int c = 0; c < block_size; c++)
        {
          this->add(offset_i + r, offset_j + c, mat->get(r, c));
        }
      }
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* SparseMatrix<Scalar>::duplicate() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("SparseMatrix* duplicate()");
    }

    template<typename Scalar>
    unsigned int SparseMatrix<Scalar>::get_nnz() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("get_nnz()");
      return 0;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::prealloc(unsigned int n)
    {
      this->size = n;
      pages = malloc_with_check<SparseMatrix<Scalar>, Page>(n, this);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      if (pages[col].count >= PAGE_SIZE)
      {
        Page* final_page = &(pages[col]);
        while (final_page->next != nullptr && final_page->count >= PAGE_SIZE)
          final_page = final_page->next;

        if (final_page->next == nullptr && final_page->count >= PAGE_SIZE)
        {
          final_page->next = new Page(true);
          final_page = final_page->next;
        }
        final_page->idx[final_page->count++] = row;
      }
      else
        pages[col].idx[pages[col].count++] = row;
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
    {
      // gather all pages in the buffer, deleting them along the way
      int *end = buffer;
      while (page != nullptr)
      {
        memcpy(end, page->idx, sizeof(int)* page->count);
        end += page->count;
        Page *tmp = page;
        page = page->next;
        if (tmp->dyn_stored)
          delete tmp;
      }

      // sort the indices and remove duplicities
      qsort_int(buffer, end - buffer);
      int *q = buffer;
      for (int *p = buffer, last = -1; p < end; p++)
      if (*p != last)
        *q++ = last = *p;

      return q - buffer;
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_indices()
    {
      int total = 0;
      for (unsigned int i = 0; i < this->size; i++)
      for (Page *page = &pages[i]; page != nullptr; page = page->next)
        total += page->count;

      return total;
    }

    template<>
    HERMES_API SparseMatrix<double>* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
      {
                                    return new CSCMatrix<double>;
      }

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
                                   if (use_direct_solver)
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
                                 if (use_direct_solver)
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
                                   return new CSCMatrix<double>;
#else
                                   throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
      {
                                          if (use_direct_solver)
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
                                   return new CSCMatrix<double>;
#else
                                   throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
                                   break;
      }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return nullptr;
    }

    template<>
    HERMES_API SparseMatrix<std::complex<double> >* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
      {
                                    return new CSCMatrix<std::complex<double> >;
      }
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
                                   if (use_direct_solver)
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
                                 if (use_direct_solver)
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
                                   return new CSCMatrix<std::complex<double> >;
#else
                                   throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
                                   break;
      }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
      {
                                          if (use_direct_solver)
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
                                   return new CSCMatrix<std::complex<double> >;
#else
                                   throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
                                   break;
      }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return nullptr;
    }

    template class Matrix<double>;
    template class Matrix<std::complex<double> >;

    template class SparseMatrix<double>;
    template class SparseMatrix<std::complex<double> >;
  }
}