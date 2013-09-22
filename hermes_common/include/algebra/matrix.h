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
/*! \file matrix.h
\brief Basic matrix classes and operations.
*/
#ifndef __HERMES_COMMON_MATRIX_H
#define __HERMES_COMMON_MATRIX_H

#include "algebra_utilities.h"
#include "algebra_mixins.h"
#include "mixins.h"

namespace Hermes
{
  /// \brief Namespace containing classes for vector / matrix operations.
  namespace Algebra
  {
    /// \brief General (abstract) matrix representation in Hermes.
    template<typename Scalar>
    class HERMES_API Matrix : public Hermes::Mixins::Loggable, public Algebra::Mixins::MatrixRhsImportExport<Scalar>
    {
    public:
      /// constructor of matrix
      /// @param[in] size size of matrix
      Matrix(unsigned int size = 0);
      virtual ~Matrix() {};

      /// allocate the memory for stiffness matrix and right-hand side
      virtual void alloc() = 0;

      /// free the memory associated with stiffness matrix and right-hand side
      virtual void free() = 0;

      /// Get the value from a position
      /// @return the value from the specified position
      /// @param[in] m - the number of row
      /// @param[in] n - the number of column
      virtual Scalar get(unsigned int m, unsigned int n) const = 0;

      /// Zero the matrix.
      virtual void zero() = 0;

      /// set the stiffness matrix
      ///
      /// @param[in] m    - the row where to set
      /// @param[in] n    - the column where to set
      /// @param[in] v    - value
      virtual void set_row_zero(unsigned int n);

      /// update the stiffness matrix
      ///
      /// @param[in] m    - the row where to update
      /// @param[in] n    - the column where to update
      /// @param[in] v    - value
      virtual void add(unsigned int m, unsigned int n, Scalar v) = 0;

      /// update the stiffness matrix
      ///
      /// @param[in] m         - number of rows of given block
      /// @param[in] n         - number of columns of given block
      /// @param[in] mat    - block of values
      /// @param[in] rows      - array with row indexes
      /// @param[in] cols      - array with column indexes
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);

      /// Add a number to each diagonal entry.
      virtual void add_to_diagonal(Scalar v);

      /// Multiply with a vector.
      virtual void multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized = false) const;

      /// Multiply with a Scalar.
      virtual void multiply_with_Scalar(Scalar value);

      /// Get size of matrix
      /// @return size of matrix
      virtual unsigned int get_size() const;

    protected:
      unsigned int size;  ///< matrix size
    };

    /// \brief General (abstract) sparse matrix representation in Hermes.
    template<typename Scalar>
    class HERMES_API SparseMatrix : public Matrix<Scalar> {
    public:
      SparseMatrix();
      /// Constructor of sparse matrix
      /// @param[in] size size of matrix
      SparseMatrix(unsigned int size);
      SparseMatrix(const SparseMatrix<Scalar>& mat);
      virtual ~SparseMatrix();

      /// prepare memory
      ///
      /// @param[in] n - number of unknowns
      virtual void prealloc(unsigned int n);

      /// add indices of nonzero matrix element
      ///
      /// @param[in] row  - row index
      /// @param[in] col  - column index
      virtual void pre_add_ij(unsigned int row, unsigned int col);

      /// Finish manipulation with matrix (called before solving)
      virtual void finish();

      /// Add matrix
      /// @param mat matrix to add
      virtual void add_sparse_matrix(SparseMatrix<Scalar>* mat);
      
      /// Add matrix to diagonal
      /// Matrices must be the same type of solver
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);

      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, SparseMatrix<Scalar>* mat);
      
      /// Return the number of entries in a specified row
      ///
      /// @param[in] row - index of the row
      /// @return - the number of entries in the row 'row'
      virtual int get_num_row_entries(unsigned int row) const;

      /// Extract the copy of a row
      ///
      /// @param[in] row - global row to extract
      /// @param[in] len - length of 'vals' and 'idxs' arrays.
      /// @param[out] n_entries - number of nonzero entries extracted.
      /// @param[out] vals - extracted values for this row.
      /// @param[out] idxs - extracted global column indices for the corresponding values.
      virtual void extract_row_copy(unsigned int row, unsigned int len,
        unsigned int &n_entries, double *vals, unsigned int *idxs) const;

      /// Return the number of entries in a specified column
      ///
      /// @param[in] col - index of the column
      /// @return - the number of entries in the column 'col'
      virtual int get_num_col_entries(unsigned int col) const;

      /// Extract the copy of a column
      ///
      /// @param[in] col - global column to extract
      /// @param[in] len - length of 'vals' and 'idxs' arrays.
      /// @param[out] n_entries - number of nonzero entries extracted.
      /// @param[out] vals - extracted values for this column.
      /// @param[out] idxs - extracted global row indices for the corresponding values.
      virtual void extract_col_copy(unsigned int col, unsigned int len,
        unsigned int &n_entries, double *vals, unsigned int *idxs) const;

      /// Duplicate sparse matrix (including allocation).
      virtual SparseMatrix<Scalar>* duplicate() const;

      /// Get fill-in.
      virtual double get_fill_in() const = 0;

      unsigned row_storage:1; ///< \todo document
      unsigned col_storage:1; ///< \todo document

      /// get number of nonzero numbers in matrix
      /// @return number of nonzero numbers in matrix
      virtual unsigned int get_nnz() const;

    protected:
      /// Size of page (max number of indices stored in one page).
      static const int PAGE_SIZE = 62;

      /// Structure for storing indices in sparse matrix
      struct Page {
        /// number of indices stored
        int count;
        /// buffer for storring indices
        int idx[PAGE_SIZE];
        /// pointer to next page
        Page *next;
      };

      /// array of pages with indices array. Each field of arra contains pages for one column
      Page **pages;

      /// gather all pages in the buffer, delete them, sort buffer and remove duplicities
      /// @param[in] page first page with indices
      /// @param[out] buffer buffer to which indices will be copied
      /// @param[in] max maximum indices to be stored (probably)
      /// \todo max parameter does nothing (not implemented)
      /// @return number of indices
      int sort_and_store_indices(Page *page, int *buffer, int *max);
      /// get number of indices in all pages
      /// @return number of indices
      int get_num_indices();

      /// mem stat
      int mem_size;
    };

    /// \brief Function returning a matrix according to the users's choice.
    /// @return created matrix
    template<typename Scalar> HERMES_API
      SparseMatrix<Scalar>*  create_matrix(bool use_direct_solver = false);
  }
}
#endif
