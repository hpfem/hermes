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
/*! \file umfpack_solver.h
\brief UMFPACK solver interface.
*/
#ifndef __HERMES_COMMON_UMFPACK_SOLVER_H_
#define __HERMES_COMMON_UMFPACK_SOLVER_H_
#include "config.h"
#ifdef WITH_UMFPACK
#include "linear_matrix_solver.h"
#include "matrix.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class HERMES_API UMFPackLinearMatrixSolver;
    template <typename Scalar> class HERMES_API UMFPackIterator;
  }

  namespace Algebra
  {
    using namespace Hermes::Solvers;
    /// \brief General CSC Matrix class.
    /// (can be used in umfpack, in that case use the
    /// UMFPackMatrix subclass, or with EigenSolver, or anything else).
    template <typename Scalar>
    class HERMES_API CSCMatrix : public SparseMatrix<Scalar>
    {
    public:
      /// Creates matrix in CSC format using size, nnz, and the three arrays.
      /// @param[in] size size of matrix (num of rows and columns)
      /// @param[in] nnz number of nonzero values
      /// @param[in] ap index to ap/ax, where each column starts (size is matrix size + 1)
      /// @param[in] ai row indices
      /// @param[in] ax values
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);

      /// \brief Default constructor.
      CSCMatrix();
      /// \brief Constructor with specific size
      /// Calls alloc.
      /// @param[in] size size of matrix (number of rows and columns)
      CSCMatrix(unsigned int size);
      virtual ~CSCMatrix();
      virtual void alloc();
      virtual void free();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(CSCMatrix<Scalar>* mat);
      /// Add matrix to diagonal.
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_to_diagonal_blocks(int num_stages, CSCMatrix<Scalar>* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, CSCMatrix<Scalar>* mat);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;

      // Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);

      // Duplicates a matrix (including allocation).
      CSCMatrix* duplicate();
      // Exposes pointers to the CSC arrays.
      /// @return pointer to #Ap
      int *get_Ap();
      // Exposes pointers to the CSC arrays.
      /// @return pointer to #Ai
      int *get_Ai();
      // Exposes pointers to the CSC arrays.
      /// @return pointer to #Ax
      Scalar *get_Ax();

    protected:
      // UMFPack specific data structures for storing the system matrix (CSC format).
      /// Matrix entries (column-wise).
      Scalar *Ax;
      /// Row indices of values in Ax.
      int *Ai;
      /// Index to Ax/Ai, where each column starts.
      int *Ap;
      /// Number of non-zero entries ( =  Ap[size]).
      unsigned int nnz;
      template <typename T> friend class Hermes::Solvers::UMFPackLinearMatrixSolver;
      template <typename T> friend class Hermes::Solvers::UMFPackIterator;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /// \brief This class is to be used with UMFPack solver only.
    template <typename Scalar>
    class HERMES_API UMFPackMatrix : public CSCMatrix<Scalar>
    {
      template <typename T> friend class Hermes::Solvers::UMFPackLinearMatrixSolver;
      template <typename T> friend class Hermes::Solvers::UMFPackIterator;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /// \brief Class representing the vector for UMFPACK.
    template <typename Scalar>
    class HERMES_API UMFPackVector : public Vector<Scalar>
    {
    public:
      UMFPackVector();
      /// Constructor of vector with specific size.
      /// @param[in] size size of vector
      UMFPackVector(unsigned int size);
      virtual ~UMFPackVector();
      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual Scalar get(unsigned int idx);
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual void change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual void add_vector(Vector<Scalar>* vec);
      virtual void add_vector(Scalar* vec);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");

      /// @return pointer to array with vector data
      /// \sa #v
      Scalar *get_c_array();

    protected:
      /// UMFPack specific data structures for storing the rhs.
      Scalar *v;
      template <typename T> friend class Hermes::Solvers::UMFPackLinearMatrixSolver;
      template <typename T> friend class Hermes::Solvers::UMFPackIterator;
      template<typename T> friend Vector<T>* Hermes::Algebra::create_vector();
    };
  }
  namespace Solvers
  {
    /// \brief Encapsulation of UMFPACK linear solver.
    ///
    /// @ingroup Solvers
    template <typename Scalar>
    class HERMES_API UMFPackLinearMatrixSolver : public DirectSolver<Scalar>
    {
    public:
      /// Constructor of UMFPack solver.
      /// @param[in] m pointer to matrix
      /// @param[in] rhs pointer to right hand side vector
      UMFPackLinearMatrixSolver(UMFPackMatrix<Scalar> *m, UMFPackVector<Scalar> *rhs);
      virtual ~UMFPackLinearMatrixSolver();
      virtual bool solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      UMFPackMatrix<Scalar> *m;
      /// Right hand side vector.
      UMFPackVector<Scalar> *rhs;

      /// \brief Reusable factorization information (A denotes matrix represented by the pointer 'm').
      /// Reordering of matrix A to reduce fill-in during factorization.
      void *symbolic;
      void *numeric;  ///< LU factorization of matrix A.

      /// \todo document
      void free_factorization_data();
      /// \todo document
      bool setup_factorization();
      template <typename T> friend class Hermes::Algebra::CSCMatrix;
      template <typename T> friend class Hermes::Algebra::UMFPackMatrix;
      template <typename T> friend class Hermes::Algebra::UMFPackVector;
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs);
      void check_status(const char *fn_name, int status);
    };

    /// \brief UMFPack matrix iterator. \todo document members
    template <typename Scalar>
    class UMFPackIterator
    {
    protected:
      UMFPackIterator(CSCMatrix<Scalar>* mat);
      bool init();
      void get_current_position(int& i, int& j, Scalar& val);
      bool move_to_position(int i, int j);
      bool move_ptr();
      void add_to_current_position(Scalar val);
    protected:
      int size;
      int nnz;
      int* Ai;
      int* Ap;
      Scalar* Ax;
      int Ai_pos;
      int Ap_pos;
      template <typename T> friend class Hermes::Algebra::CSCMatrix;
      template <typename T> friend class Hermes::Algebra::UMFPackMatrix;
      template <typename T> friend class Hermes::Algebra::UMFPackVector;
    };
  }
}
#endif
#endif