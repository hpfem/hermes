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
/*! \file dense_matrix_operations.h
\brief Dense (small) simply stored matrix operations.
*/
#ifndef __HERMES_COMMON_DENSE_MATRIX_OPERATIONS_H
#define __HERMES_COMMON_DENSE_MATRIX_OPERATIONS_H

#include "common.h"
#include "util/compat.h"
#include "exceptions.h"

namespace Hermes
{
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
      template<typename T>
      HERMES_API void ludcmp(T **a, int n, int *indx, double *d);

      /// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
      /// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
      /// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
      /// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
      /// and can be left in place for successive calls with different right-hand sides b. This routine takes
      /// into account the possibility that b will begin with many zero elements, so it is efficient for use
      /// in matrix inversion.
      template<typename T, typename S>
      HERMES_API void lubksb(T **a, int n, int *indx, S *b);

      /// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
      /// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
      /// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
      /// elements which are returned in p[n].
      template<typename T>
      HERMES_API void choldc(T **a, int n, T p[]);

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
  }
}
#endif
