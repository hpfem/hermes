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

#ifndef __HERMES_COMMON_MATRIX_H
#define __HERMES_COMMON_MATRIX_H

#include "common.h"
#include "error.h"

/// Creates a new (full) matrix with m rows and n columns with entries of the type T.
/// The entries can be accessed by matrix[i][j]. To delete the matrix, just
/// do "delete matrix".
template<typename T>
T **new_matrix(int m, int n = 0)
{
  if (!n) n = m;
  T **vec = (T **) new char[sizeof(T *) * m + sizeof(T) * m * n];
  MEM_CHECK(vec);
  memset(vec, 0, sizeof(T *) * m + sizeof(T) * m * n);
  T *row = (T *) (vec + m);
  for (int i = 0; i < m; i++, row += n) vec[i] = row;
  return vec;
}

/// Copies a matrix. Both matrices has to be equal to or larger than provideded sizes.
/// Size compatibility check is not done.
template<typename T>
void copy_matrix(T** dest, T** src, int m, int n = 0) 
{
  if (n == 0) n = m;
  for(int i = 0; i < m; i++) {
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
void save_matrix_octave(const std::string& matrix_name, T** matrix, int m, int n = 0, const std::string& filename = std::string()) 
{
  if (n == 0) n = m;

  //create filename
  std::string fname = filename;
  if (fname.empty())
    fname = matrix_name + ".mat";

  //open file
  std::ofstream fout(fname.c_str());
  if (!fout.is_open()) {
    error("Unable to save a matrix to a file \"%s\"", fname.c_str());
    return;
  }

  //write header
  fout << std::string("# name: ") << matrix_name << std::endl;
  fout << std::string("# type: matrix") << std::endl;
  fout << std::string("# rows: ") << m << std::endl;
  fout << std::string("# columns: ") << n << std::endl;

  //write contents
  for(int i = 0; i < m; i++) {
    for(int k = 0; k < n; k++)
      fout << ' ' << matrix[i][k];
    fout << std::endl;
  }

  //finish
  fout.close();
}

/// Saves MxM sparse matrix to a octave file format.
template<typename T>
void save_sparse_matrix_octave(const std::string& matrix_name, const T* Ax, const int* Ap, const int* Ai, 
                               int m, const std::string& filename = std::string()) 
{
  // create filename
  std::string fname = filename;
  if (fname.empty())
    fname = matrix_name + ".mat";

  // open file
  std::ofstream fout(fname.c_str());
  if (!fout.is_open()) {
    error("Unable to save a matrix to a file \"%s\"", fname.c_str());
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
    for (int i = Ap[j]; i < Ap[j+1]; i++)
      fout << j+1 << " " << Ai[i]+1 << " " << Ax[i] << std::endl;

  // finish
  fout.close();
}

/// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
/// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
template<typename T>
void transpose(T **matrix, int m, int n) 
{
  int min = std::min(m, n);
  for (int i = 0; i < min; i++)
    for (int j = i+1; j < min; j++)
       std::swap(matrix[i][j], matrix[j][i]);

  if (m < n)
    for (int i = 0; i < m; i++)
      for (int j = m; j < n; j++)
	matrix[j][i] = matrix[i][j];
  else if (n < m)
    for (int i = n; i < m; i++)
      for (int j = 0; j < n; j++)
	matrix[j][i] = matrix[i][j];
}


/// Changes the sign of a matrix
template<typename T>
void chsgn(T **matrix, int m, int n) 
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      matrix[i][j] = -matrix[i][j];
}

/// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
/// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
/// indx[n] is an output vector that records the row permutation effected by the partial
/// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
/// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
/// or invert a matrix.
void ludcmp(double **a, int n, int *indx, double *d);

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

  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}


/// Simple dot product.
double HERMES_API vec_dot(double *r, double *s, int ndof);

class Vector;
double HERMES_API vec_dot(Vector *r, Vector *s, int ndof);

/// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
/// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
/// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
/// elements which are returned in p[n].
void choldc(double **a, int n, double p[]);

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

  for (i = 0; i < n; i++) {
    sum = b[i];
    k = i;
    while (--k >= 0) sum -= a[i][k] * x[k];
    x[i] = sum / p[i];
  }

  for (i = n-1; i >= 0; i--) {
    sum = x[i];
    k = i;
    while (++k < n) sum -= a[k][i] * x[k];
    x[i] = sum / p[i];
  }
}

//  ///////

enum EMatrixDumpFormat {
	DF_MATLAB_SPARSE,
	DF_PLAIN_ASCII,
	DF_HERMES_BIN,
	DF_NATIVE					// native format for the linear solver
};

class HERMES_API Matrix {
public:
  int get_size() { return this->size;};

  Matrix(int size) { this->size = size;};

  ~Matrix() {};

  Matrix() { this->size = 0;};

  /// allocate the memory for stiffness matrix and right-hand side
  virtual void alloc() = 0;

  /// free the memory associated with stiffness matrix and right-hand side
  virtual void free() = 0;

  /// Get the value from a position
  /// @return the value from the specified position
  /// @param[in] m - the number of row
  /// @param[in] n - the number of column
  virtual scalar get(int m, int n) = 0;

  /// Zero the matrix.
  virtual void zero() = 0;

  /// Add a number to each diagonal entry.
  virtual void add_to_diagonal(scalar v) = 0;

  /// update the stiffness matrix
  ///
  /// @param[in] m    - the row where to update
  /// @param[in] n    - the column where to update
  /// @param[in] v    - value
  virtual void add(int m, int n, scalar v) = 0;

  /// update the stiffness matrix
  ///
  /// @param[in] m         - number of rows of given block
  /// @param[in] n         - number of columns of given block
  /// @param[in] matrix    - block of values
  /// @param[in] rows      - array with row indexes
  /// @param[in] cols      - array with column indexes
  virtual void add(int m, int n, scalar **mat, int *rows, int *cols) = 0;

  /// dumping matrix and right-hand side
  ///
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat = DF_MATLAB_SPARSE) = 0;

  virtual int get_matrix_size() const = 0;

protected:

  int size;  // matrix size
};

class HERMES_API SparseMatrix : public Matrix {
public:
  SparseMatrix();
  SparseMatrix(int size);
  virtual ~SparseMatrix();

  /// prepare memory
  ///
  /// @param[in] ndofs - number of unknowns
  virtual void prealloc(int n);

  /// add indices of nonzero matrix element
  ///
  /// @param[in] row  - row index
  /// @param[in] col  - column index
  virtual void pre_add_ij(int row, int col);

  virtual void finish() { }

  virtual int get_size() { return size; }

  virtual void add_sparse_matrix(SparseMatrix* mat) 
  { 
    error("add_sparse_matrix() undefined.");
  };

  /// Return the number of entries in a specified row
  ///
  /// @param[in] row - index of the row
  /// @return - the number of entries in the row 'row'
  virtual int get_num_row_entries(int row) { return -1; }

  /// Extract the copy of a row
  ///
  /// @param[in] row - global row to extract
  /// @param[in] len - length of 'vals' and 'idxs' arrays.
  /// @param[out] n_entries - number of nonzero entries extracted.
  /// @param[out] vals - extracted values for this row.
  /// @param[out] idxs - extracted global column indices for the corresponding values.
  virtual void extract_row_copy(int row, int len, int &n_entries, double *vals, int *idxs) { }

  /// Return the number of entries in a specified column
  ///
  /// @param[in] row - index of the column
  /// @return - the number of entries in the column 'col'
  virtual int get_num_col_entries(int col) { return -1; }

  /// Extract the copy of a column
  ///
  /// @param[in] row - global column to extract
  /// @param[in] len - length of 'vals' and 'idxs' arrays.
  /// @param[out] n_entries - number of nonzero entries extracted.
  /// @param[out] vals - extracted values for this column.
  /// @param[out] idxs - extracted global row indices for the corresponding values.
  virtual void extract_col_copy(int col, int len, int &n_entries, double *vals, int *idxs) { }

  // virtual function for multiplying a vector by a 
  // derived type of SparseMatrix
  virtual void multiply(scalar* vector_in, scalar* vector_out){ };
	
  //
  virtual double get_fill_in() const = 0;

  unsigned row_storage:1;
  unsigned col_storage:1;

protected:
  static const int PAGE_SIZE = 62;

  struct Page {
    int count;
    int idx[PAGE_SIZE];
    Page *next;
  };

  Page **pages;

  int sort_and_store_indices(Page *page, int *buffer, int *max);
  int get_num_indices();

  // mem stat
  int mem_size;
};

class HERMES_API Vector {
public:
  virtual ~Vector() { }

  /// allocate memory for storing ndofs elements
  ///
  /// @param[in] ndofs - number of elements of the vector
  virtual void alloc(int ndofs) = 0;
  /// free the memory
  virtual void free() = 0;
  // finish the assembly of the vector
  virtual void finish() { }

  /// Get the value from a position
  /// @return the value form the specified index
  /// @param[in] idx - index which to obtain the value from
  virtual scalar get(int idx) = 0;

  /// Extract vector values into user-provided array.
  /// @param[out] v - array which will contain extracted values
  virtual void extract(scalar *v) const = 0;

  /// Zero the vector
  virtual void zero() = 0;

  /// Multiply by minus one.
  virtual void change_sign() = 0;

  /// set the entry on a specified position
  ///
  /// @param[in] idx - indices where to update
  /// @param[in] y   - value
  virtual void set(int idx, scalar y) = 0;

  /// update element on the specified position
  ///
  /// @param[in] idx - indices where to update
  /// @param[in] y   - value
  virtual void add(int idx, scalar y) = 0;

  /// Add a vector.
  virtual void add_vector(Vector* vec) = 0;
  virtual void add_vector(scalar* vec) = 0;

  /// update subset of the elements
  ///
  /// @param[in] n   - number of positions to update
  /// @param[in] idx - indices where to update
  /// @param[in] y   - values
  virtual void add(int n, int *idx, scalar *y) = 0;

  /// Get vector length.
  int length() {return this->size;}

  // Write to file.
  virtual bool dump(FILE *file, const char *var_name, 
                    EMatrixDumpFormat = DF_MATLAB_SPARSE) = 0;
  
protected:
  int size;
};

HERMES_API Vector* create_vector(MatrixSolverType matrix_solver);
HERMES_API SparseMatrix*  create_matrix(MatrixSolverType matrix_solver);

class Solver;
HERMES_API Solver*  create_linear_solver(MatrixSolverType matrix_solver, 
                                         Matrix* matrix, Vector* rhs = NULL);

#endif
