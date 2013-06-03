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
/*! \file cs_matrix.h
\brief Basic cs (Compressed sparse) matrix classes and operations.
*/
#ifndef __HERMES_COMMON_CS_MATRIX_H
#define __HERMES_COMMON_CS_MATRIX_H

#include "matrix.h"

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class HERMES_API CSIterator;
    template <typename Scalar> class HERMES_API CSCIterator;
    template <typename Scalar> class HERMES_API CSRIterator;
  }

  /// \brief Namespace containing classes for vector / matrix operations.
  namespace Algebra
  {
    /// \brief General CS Matrix class.
    /// Either row- or column- specific (see subclassses).
    template <typename Scalar>
    class HERMES_API CSMatrix : public SparseMatrix<Scalar>
    {
    public:
      /// Creates matrix in CS format using size, nnz, and the three arrays.
      /// @param[in] size size of matrix (num of rows and columns)
      /// @param[in] nnz number of nonzero values
      /// @param[in] ap index to ap/ax, where each column /row starts (size is matrix size + 1)
      /// @param[in] ai row / column indices
      /// @param[in] ax values
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);

      /// \brief Default constructor.
      CSMatrix();
      /// \brief Constructor with specific size
      /// Calls alloc.
      /// @param[in] size size of matrix (number of rows and columns)
      CSMatrix(unsigned int size);
      virtual ~CSMatrix();

      /// Main addition method.
      /// Virtual - the method body is 1:1 for CSCMatrix, inverted for CSR.
      virtual void add(unsigned int Ai_data_index, unsigned int Ai_index, Scalar v);
      /// Main get method.
      /// Virtual - the method body is 1:1 for CSCMatrix, inverted for CSR.
      virtual Scalar get(unsigned int Ai_data_index, unsigned int Ai_index);
      
      virtual void add_to_diagonal(Scalar v);
      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(CSMatrix<Scalar>* mat);
      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, CSMatrix<Scalar>* mat);
      /// Add matrix to diagonal.
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_to_diagonal_blocks(int num_stages, CSMatrix<Scalar>* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);

      /// Utility method.
      virtual void alloc();
      /// Utility method.
      virtual void free();
      /// Utility method.
      virtual void zero();
      /// Utility method.
      virtual void set_row_zero(unsigned int n);

      /// Utility method.
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf") = 0;
      /// Utility method.
      virtual unsigned int get_matrix_size() const;
      /// Utility method.
      virtual unsigned int get_nnz() const;
      /// Utility method.
      virtual double get_fill_in() const;

      // Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);

      // Duplicates a matrix (including allocation).
      virtual CSMatrix* duplicate();
      // Exposes pointers to the CS arrays.
      /// @return pointer to #Ap
      int *get_Ap();
      // Exposes pointers to the CS arrays.
      /// @return pointer to #Ai
      int *get_Ai();
      // Exposes pointers to the CS arrays.
      /// @return pointer to #Ax
      Scalar *get_Ax();

    protected:
      // UMFPack specific data structures for storing the system matrix (CSC format).
      /// Matrix entries (column-wise).
      Scalar *Ax;
      /// Row / Column indices of values in Ax.
      int *Ai;
      /// Index to Ax/Ai, where each column / row starts.
      int *Ap;
      /// Number of non-zero entries ( =  Ap[size]).
      unsigned int nnz;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /// \brief General CSC Matrix class.
    /// (can be used in umfpack, in that case use the
    /// UMFPackMatrix subclass, or with EigenSolver, or anything else).
    template <typename Scalar>
    class HERMES_API CSCMatrix : public CSMatrix<Scalar>
    {
    public:
      /// \brief Default constructor.
      CSCMatrix();

      /// \brief Constructor with specific size
      /// Calls alloc.
      /// @param[in] size size of matrix (number of rows and columns)
      CSCMatrix(unsigned int size);

      virtual ~CSCMatrix();

      virtual Scalar get(unsigned int m, unsigned int n);

      virtual void add(unsigned int m, unsigned int n, Scalar v);

      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");

      friend class Hermes::Solvers::CSCIterator<Scalar>;

      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(CSMatrix<Scalar>* mat);

      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, CSMatrix<Scalar>* mat);

      // Duplicates a matrix (including allocation).
      virtual CSMatrix<Scalar>* duplicate();
    };

    /// \brief General CSR Matrix class.
    /// (can be used in umfpack, in that case use the
    /// UMFPackMatrix subclass, or with EigenSolver, or anything else).
    template <typename Scalar>
    class HERMES_API CSRMatrix : public CSMatrix<Scalar>
    {
    public:
      /// \brief Default constructor.
      CSRMatrix();

      /// \brief Constructor with specific size
      /// Calls alloc.
      /// @param[in] size size of matrix (number of rows and columns)
      CSRMatrix(unsigned int size);

      virtual ~CSRMatrix();

      virtual Scalar get(unsigned int m, unsigned int n);

      virtual void add(unsigned int m, unsigned int n, Scalar v);

      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");

      friend class Hermes::Solvers::CSCIterator<Scalar>;

      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(CSMatrix<Scalar>* mat);

      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, CSMatrix<Scalar>* mat);

      // Duplicates a matrix (including allocation).
      virtual CSMatrix<Scalar>* duplicate();

      /// Important - normal SparseMatrix has the pages structure suitable for CSC matrix, so we need
      /// to override the structure creation here.
      /// add indices of nonzero matrix element
      ///
      /// @param[in] row  - row index
      /// @param[in] col  - column index
      virtual void pre_add_ij(unsigned int row, unsigned int col);
    };
  }
  
  namespace Solvers
  {
    /// \brief CS matrix iterator. \todo document members
    template <typename Scalar>
    class CSIterator
    {
    protected:
      CSIterator(Hermes::Algebra::CSMatrix<Scalar>* mat);
      bool init();
      virtual void get_current_position(int& i, int& j, Scalar& val) = 0;
      virtual bool move_to_position(int i, int j) = 0;
      bool move_ptr();
      void add_to_current_position(Scalar val);

      int size;
      int nnz;
      int* Ai;
      int* Ap;
      Scalar* Ax;
      int Ai_pos;
      int Ap_pos;

      friend class Hermes::Algebra::CSMatrix<Scalar>;
    };

    /// \brief CSC matrix iterator. \todo document members
    template <typename Scalar>
    class CSCIterator : public CSIterator<Scalar>
    {
    protected:
      CSCIterator(Hermes::Algebra::CSCMatrix<Scalar>* mat);
      void get_current_position(int& i, int& j, Scalar& val);
      bool move_to_position(int i, int j);

      friend class Hermes::Algebra::CSCMatrix<Scalar>;
    };

    /// \brief CSC matrix iterator. \todo document members
    template <typename Scalar>
    class CSRIterator : public CSIterator<Scalar>
    {
    protected:
      CSRIterator(Hermes::Algebra::CSRMatrix<Scalar>* mat);
      void get_current_position(int& i, int& j, Scalar& val);
      bool move_to_position(int i, int j);

      friend class Hermes::Algebra::CSCMatrix<Scalar>;
    };
  }
}
#endif
