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

#include "common.h"
#include "vector.h"
#include "exceptions.h"
#include "mixins.h"

namespace Hermes
{
  enum MatrixSolverType
  {
    SOLVER_UMFPACK = 0,
    SOLVER_PETSC,
    SOLVER_MUMPS,
    SOLVER_SUPERLU,
    SOLVER_AMESOS,
    SOLVER_AZTECOO
  };

  /// \brief Namespace containing classes for vector / matrix operations.
  namespace Algebra
  {
    /// Contains operation on dense matrices.
    namespace DenseMatrixOperations
    {
      /// Creates a new (full) matrix with m rows and n columns with entries of the type T.
      /// The entries can be accessed by matrix[i][j]. To delete the matrix, just
      /// do "delete matrix".
      template<typename T>
      T **new_matrix(unsigned int m, unsigned int n = 0)
      {
        if(!n) n = m;
        T **vec = (T **) new char[sizeof(T *) * m + sizeof(T) * m * n];
        memset(vec, 0, sizeof(T *) * m + sizeof(T) * m * n);
        T *row = (T *) (vec + m);
        for (unsigned int i = 0; i < m; i++, row += n) vec[i] = row;
        return vec;
      }

      template<typename T>
      T **new_matrix_malloc(unsigned int m, unsigned int n = 0)
      {
        if(!n) n = m;
        T **vec = (T **) malloc(sizeof(T *) * m + sizeof(T) * m * n);
        memset(vec, 0, sizeof(T *) * m + sizeof(T) * m * n);
        T *row = (T *) (vec + m);
        for (unsigned int i = 0; i < m; i++, row += n) vec[i] = row;
        return vec;
      }

      /// Copies a matrix. Both matrices has to be equal to or larger than provided sizes.
      /// Size compatibility check is not done.
      template<typename T>
      void copy_matrix(T** dest, T** src, unsigned int m, unsigned int n = 0)
      {
        if(n == 0) n = m;
        for(unsigned int i = 0; i < m; i++)
        {
          memcpy(dest[i], src[i], n*sizeof(T));
        }
      }

      /// \brief Saves a dense matrix to a octave file format.
      /// \param[in] matrix_name A name of a matrix in Octave. It can be used to create an output filename.
      /// \param[in] matrix A pointer to an array of pointers to a rows of the matrix. Such a structure can be generated using new_matrix() or it can be a pointer to an 1D C-array.
      /// \param[in] m A number of rows of the matrix. Set to 1 if the matrix is a pointer to a 1D C-array.
      /// \param[in] n A number of columns of the matrix. If zero, it is assumed to be equal to m.
      /// \param[in] filename An output filename. If not specified, matrix_name will be used by concatenating it with a suffix '.mat'.
      template<typename T>
      void save_matrix_octave(const std::string& matrix_name, T** matrix, unsigned int m, unsigned int n = 0, const std::string& filename = std::string())
      {
        if(n == 0) n = m;

        //create filename
        std::string fname = filename;
        if(fname.empty())
          fname = matrix_name + ".mat";

        //open file
        std::ofstream fout(fname.c_str());
        if(!fout.is_open())
        {
          throw Hermes::Exceptions::Exception("Unable to save a matrix to a file \"%s\"", fname.c_str());
          return;
        }

        //write header
        fout << std::string("# name: ") << matrix_name << std::endl;
        fout << std::string("# type: matrix") << std::endl;
        fout << std::string("# rows: ") << m << std::endl;
        fout << std::string("# columns: ") << n << std::endl;

        //write contents
        for(unsigned int i = 0; i < m; i++)
        {
          for(unsigned int k = 0; k < n; k++)
            fout << ' ' << matrix[i][k];
          fout << std::endl;
        }

        //finish
        fout.close();
      }

      /// Saves MxM sparse matrix to a octave file format.
      template<typename T>
      void save_sparse_matrix_octave(const std::string& matrix_name, const T* Ax, const int* Ap, const int* Ai,
        unsigned int m, const std::string& filename = std::string())
      {
        // create filename
        std::string fname = filename;
        if(fname.empty())
          fname = matrix_name + ".mat";

        // open file
        std::ofstream fout(fname.c_str());
        if(!fout.is_open())
        {
          throw Hermes::Exceptions::Exception("Unable to save a matrix to a file \"%s\"", fname.c_str());
          return;
        }

        // write header
        fout << std::string("# name: ") << matrix_name << std::endl;
        fout << std::string("# type: sparse matrix") << std::endl;
        fout << std::string("# nnz: ") << Ap[m] << std::endl;
        fout << std::string("# rows: ") << m << std::endl;
        fout << std::string("# columns: ") << m << std::endl;

        // write contents
        for (int j = 0; j < m; j++)
          for (int i = Ap[j]; i < Ap[j + 1]; i++)
            fout << j + 1 << " " << Ai[i] + 1 << " " << Ax[i] << std::endl;

        // finish
        fout.close();
      }

      /// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
      /// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
      template<typename T>
      void transpose(T **matrix, unsigned int m, unsigned int n)
      {
        unsigned int min = std::min(m, n);
        for (unsigned int i = 0; i < min; i++)
          for (unsigned int j = i + 1; j < min; j++)
            std::swap(matrix[i][j], matrix[j][i]);

        if(m < n)
          for (unsigned int i = 0; i < m; i++)
            for (unsigned int j = m; j < n; j++)
              matrix[j][i] = matrix[i][j];
        else if(n < m)
          for (unsigned int i = n; i < m; i++)
            for (unsigned int j = 0; j < n; j++)
              matrix[j][i] = matrix[i][j];
      }

      /// Changes the sign of a matrix
      template<typename T>
      void chsgn(T **matrix, unsigned int m, unsigned int n)
      {
        for (unsigned int i = 0; i < m; i++)
          for (unsigned int j = 0; j < n; j++)
            matrix[i][j] = -matrix[i][j];
      }

      /// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
      /// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
      /// indx[n] is an output vector that records the row permutation effected by the partial
      /// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
      /// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
      /// or invert a matrix.
      HERMES_API void ludcmp(double **a, int n, int *indx, double *d);

      /// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
      /// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
      /// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
      /// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
      /// and can be left in place for successive calls with different right-hand sides b. This routine takes
      /// into account the possibility that b will begin with many zero elements, so it is efficient for use
      /// in matrix inversion.
      template<typename T>
      void lubksb(double **a, int n, int *indx, T *b)
      {
        int i, ip, j;
        T sum;

        for (i = 0; i < n; i++)
        {
          ip = indx[i];
          sum = b[ip];
          b[ip] = b[i];
          for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
          b[i] = sum;
        }
        for (i = n-1; i >= 0; i--)
        {
          sum = b[i];
          for (j = i + 1; j < n; j++) sum -= a[i][j]*b[j];
          b[i] = sum / a[i][i];
        }
      }

      /// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
      /// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
      /// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
      /// elements which are returned in p[n].
      HERMES_API void choldc(double **a, int n, double p[]);

      /// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
      /// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
      /// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
      /// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
      /// for successive calls with different right-hand sides b. b is not modified unless you identify b and
      /// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
      /// the solution x is also complex.
      template<typename T>
      void cholsl(double **a, int n, double p[], T b[], T x[])
      {
        int i, k;
        T sum;

        for (i = 0; i < n; i++)
        {
          sum = b[i];
          k = i;
          while (--k >= 0) sum -= a[i][k] * x[k];
          x[i] = sum / p[i];
        }

        for (i = n-1; i >= 0; i--)
        {
          sum = x[i];
          k = i;
          while (++k < n) sum -= a[k][i] * x[k];
          x[i] = sum / p[i];
        }
      }
    }

    /// Format of file matrix and vector output
    enum EMatrixDumpFormat
    {
      DF_MATLAB_SPARSE, ///< matlab file
      /// \brief plain ascii file
      /// first line is matrix size
      /// second line in number of nonzero values
      /// next lines contains row column and value
      DF_PLAIN_ASCII,
      /// \brief Hermes binary format
      ///
      DF_HERMES_BIN,
      DF_NATIVE,    ///< native format for the linear solver,
      DF_MATRIX_MARKET ///< Matrix Market which can be read by pysparse library
    };

    /// \brief General (abstract) matrix representation in Hermes.
    template<typename Scalar>
    class HERMES_API Matrix : public Hermes::Mixins::Loggable
    {
    public:
      /// get size of matrix
      /// @return size of matrix
      unsigned int get_size() { return this->size;};

      /// constructor of matrix
      /// @param[in] size size of matrix
      Matrix(unsigned int size) { this->size = size;};

      virtual ~Matrix() {};

      Matrix() { this->size = 0;};

      /// allocate the memory for stiffness matrix and right-hand side
      virtual void alloc() = 0;

      /// free the memory associated with stiffness matrix and right-hand side
      virtual void free() = 0;

      /// Get the value from a position
      /// @return the value from the specified position
      /// @param[in] m - the number of row
      /// @param[in] n - the number of column
      virtual Scalar get(unsigned int m, unsigned int n) = 0;

      /// Zero the matrix.
      virtual void zero() = 0;

      /// Add a number to each diagonal entry.
      virtual void add_to_diagonal(Scalar v) = 0;

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
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols) = 0;

      /// dumping matrix and right-hand side
      /// @param[in] file file handle
      /// @param[in] var_name name of variable (will be written to output file)
      /// @param[in] fmt output file format
      /// @return true on succes
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE) = 0;

      /// Get size of matrix
      /// @return size of matrix
      virtual unsigned int get_matrix_size() const = 0;

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
      virtual void finish() { }

      virtual unsigned int get_size() { return this->size; }

      /// Add matrix
      /// @param mat matrix to add
      virtual void add_sparse_matrix(SparseMatrix* mat)
      {
        throw Hermes::Exceptions::Exception("add_sparse_matrix() undefined.");
      };

      /// Add matrix to diagonal
      /// Matrices must be the same type of solver
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
      {
        throw Hermes::Exceptions::Exception("add_sparse_to_diagonal_blocks() undefined.");
      };

      /// Return the number of entries in a specified row
      ///
      /// @param[in] row - index of the row
      /// @return - the number of entries in the row 'row'
      virtual int get_num_row_entries(unsigned int row) { return -1; }

      /// Extract the copy of a row
      ///
      /// @param[in] row - global row to extract
      /// @param[in] len - length of 'vals' and 'idxs' arrays.
      /// @param[out] n_entries - number of nonzero entries extracted.
      /// @param[out] vals - extracted values for this row.
      /// @param[out] idxs - extracted global column indices for the corresponding values.
      virtual void extract_row_copy(unsigned int row, unsigned int len,
        unsigned int &n_entries, double *vals,
        unsigned int *idxs) { }

      /// Return the number of entries in a specified column
      ///
      /// @param[in] col - index of the column
      /// @return - the number of entries in the column 'col'
      virtual int get_num_col_entries(unsigned int col) { return -1; }

      /// Extract the copy of a column
      ///
      /// @param[in] col - global column to extract
      /// @param[in] len - length of 'vals' and 'idxs' arrays.
      /// @param[out] n_entries - number of nonzero entries extracted.
      /// @param[out] vals - extracted values for this column.
      /// @param[out] idxs - extracted global row indices for the corresponding values.
      virtual void extract_col_copy(unsigned int col, unsigned int len,
        unsigned int &n_entries, double *vals,
        unsigned int *idxs) { }

      /// Multiply with a vector.
      virtual void multiply_with_vector(Scalar* vector_in, Scalar* vector_out) {
        throw Hermes::Exceptions::Exception("multiply_with_vector() undefined.");
      };

      /// Multiply with a Scalar.
      virtual void multiply_with_Scalar(Scalar value) {
        throw Hermes::Exceptions::Exception("multiply_with_Scalar() undefined.");
      };

      /// Duplicate sparse matrix (including allocation).
      virtual SparseMatrix* duplicate() { return (SparseMatrix*)NULL;};

      /// Get fill-in.
      virtual double get_fill_in() const = 0;

      unsigned row_storage:1; ///< \todo document
      unsigned col_storage:1; ///< \todo document

      /// get number of nonzero numbers in matrix
      /// @return number of nonzero numbers in matrix
      virtual unsigned int get_nnz() const {
        throw Hermes::Exceptions::Exception("get_nnz() undefined.");
        return 0;
      }

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

    /// \brief General (abstract) vector representation in Hermes.
    template<typename Scalar>
    class HERMES_API Vector : public Hermes::Mixins::Loggable
    {
    public:
      virtual ~Vector() { }

      /// allocate memory for storing ndofs elements
      ///
      /// @param[in] ndofs - number of elements of the vector
      virtual void alloc(unsigned int ndofs) = 0;
      /// free the memory
      virtual void free() = 0;
      /// finish the assembly of the vector
      virtual void finish() { }

      /// Get the value from a position
      /// @return the value form the specified index
      /// @param[in] idx - index which to obtain the value from
      virtual Scalar get(unsigned int idx) = 0;

      /// Extract vector values into user-provided array.
      /// @param[out] v - array which will contain extracted values
      virtual void extract(Scalar *v) const = 0;

      /// Zero the vector
      virtual void zero() = 0;

      /// Multiply by minus one.
      virtual void change_sign() = 0;

      /// set the entry on a specified position
      ///
      /// @param[in] idx - indices where to update
      /// @param[in] y   - value
      virtual void set(unsigned int idx, Scalar y) = 0;

      /// update element on the specified position
      ///
      /// @param[in] idx - indices where to update
      /// @param[in] y   - value
      virtual void add(unsigned int idx, Scalar y) = 0;

      /// Add a vector.
      virtual void add_vector(Vector<Scalar>* vec) = 0;
      /// Add a vector.
      virtual void add_vector(Scalar* vec) = 0;

      /// update subset of the elements
      ///
      /// @param[in] n   - number of positions to update
      /// @param[in] idx - indices where to update
      /// @param[in] y   - values
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y) = 0;

      /// Get vector length.
      unsigned int length() const {return this->size;}

      /// Write to file.
      /// @param[in] file file handle
      /// @param[in] var_name name of variable (will be written to output file)
      /// @param[in] fmt output file format
      /// @return true on succes
      virtual bool dump(FILE *file, const char *var_name,
        EMatrixDumpFormat fmt = DF_MATLAB_SPARSE) = 0;

    protected:
      /// size of vector
      unsigned int size;
    };

    /// \brief Function returning a vector according to the users's choice.
    /// @return created vector
    template<typename Scalar> HERMES_API
      Vector<Scalar>* create_vector();

    /// \brief Function returning a matrix according to the users's choice.
    /// @return created matrix
    template<typename Scalar> HERMES_API
      SparseMatrix<Scalar>*  create_matrix();
  }
}
#endif