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
      if(!vector_out_initialized)
        vector_out = new Scalar[this->size];
      for(int i = 0; i < this->size; i++)
      {
        vector_out[i] = Scalar(0.);
        for(int j = 0; j < this->size; j++)
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

    template<typename Scalar>
    void Matrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      for (unsigned int i = 0; i < m; i++)       // rows
        for (unsigned int j = 0; j < n; j++)     // cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
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
    SparseMatrix<Scalar>::SparseMatrix(const SparseMatrix<Scalar>& mat)
    {
      this->size = mat.get_size();

      if (mat.pages)
      {
        this->pages = new Page *[this->size];
        memset(this->pages, 0, this->size * sizeof(Page *));

        for (int col = 0; col < this->size; col++)
        {
          Page *page = mat.pages[col];
          Page *new_page = new Page;
          new_page->count = page->count;
          memcpy(new_page->idx, page->idx, sizeof(int) * page->count);
          new_page->next = this->pages[col];
          this->pages[col] = new_page; 
        }
      }
      else
      {
        this->pages = NULL;
      }

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
      if(this->get_size() != (unsigned int) num_stages * ndof)
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in SparseMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
        this->add_as_block(ndof*i, ndof*i, mat);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, SparseMatrix<Scalar>* mat)
    {
      if((this->get_size() < offset_i + mat->get_size() )||(this->get_size() < offset_j + mat->get_size() ))
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

    template class Matrix<double>;
    template class Matrix<std::complex<double> >;

    template class SparseMatrix<double>;
    template class SparseMatrix<std::complex<double> >;
  }
}
