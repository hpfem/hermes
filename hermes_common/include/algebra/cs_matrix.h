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

#include "algebra/matrix.h"

namespace Hermes
{
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

      /// Switches CSR / CSC arrays.
      void switch_orientation();

      /// Main addition method.
      /// Virtual - the method body is 1:1 for CSCMatrix, inverted for CSR.
      virtual void add(unsigned int Ai_data_index, unsigned int Ai_index, Scalar v);
      
      /// Main get method.
      /// Virtual - the method body is 1:1 for CSCMatrix, inverted for CSR.
      virtual Scalar get(unsigned int Ai_data_index, unsigned int Ai_index) const;

      /// Utility method.
      virtual void alloc();
      /// Utility method.
      virtual void free();
      /// Utility method.
      virtual void zero();
      /// Utility method.
      virtual void set_row_zero(unsigned int n);

      /// Matrix export method.
      /// Utility version
      /// \See Matrix<Scalar>::export_to_file.
      void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf", bool invert_storage = false);
      
      /// Reading matrix
      /// Utility version
      /// \See Matrix<Scalar>::import_from_file.
      void import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt, bool invert_storage = false);
      
      /// Utility method.
      virtual unsigned int get_nnz() const;
      /// Utility method.
      virtual double get_fill_in() const;

      /// Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);

      /// Exposes pointers to the CS arrays.
      /// @return pointer to #Ap
      int *get_Ap() const;
      /// Exposes pointers to the CS arrays.
      /// @return pointer to #Ai
      int *get_Ai() const;
      /// Exposes pointers to the CS arrays.
      /// @return pointer to #Ax
      Scalar *get_Ax() const;

    protected:
      /// UMFPack specific data structures for storing the system matrix (CSC format).
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
    /// CSCMatrix subclass, or with EigenSolver, or anything else).
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

      virtual Scalar get(unsigned int m, unsigned int n) const;

      virtual void add(unsigned int m, unsigned int n, Scalar v);

      void multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const;

      void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");
      void import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt);

      /// Duplicates a matrix (including allocation).
      virtual CSMatrix<Scalar>* duplicate() const;
    };

    /// \brief General CSR Matrix class.
    /// (can be used in umfpack, in that case use the
    /// CSCMatrix subclass, or with EigenSolver, or anything else).
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

      virtual Scalar get(unsigned int m, unsigned int n) const;

      virtual void add(unsigned int m, unsigned int n, Scalar v);

      void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");
      void import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt);

      /// Duplicates a matrix (including allocation).
      virtual SparseMatrix<Scalar>* duplicate() const;

      /// Important - normal SparseMatrix has the pages structure suitable for CSC matrix, so we need
      /// to override the structure creation here.
      /// add indices of nonzero matrix element
      ///
      /// @param[in] row  - row index
      /// @param[in] col  - column index
      virtual void pre_add_ij(unsigned int row, unsigned int col);
    };
  }
}
#endif
