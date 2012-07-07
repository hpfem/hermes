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
/*! \file mumps_solver.h
\brief MUMPS solver interface.
*/
#ifndef __HERMES_COMMON_MUMPS_SOLVER_H_
#define __HERMES_COMMON_MUMPS_SOLVER_H_
#include "config.h"
#ifdef WITH_MUMPS
#include "linear_matrix_solver.h"
#include "matrix.h"

extern "C"
{
#include <mumps_c_types.h>
#include <dmumps_c.h>
#include <zmumps_c.h>
}

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class MumpsSolver;
  }
}
namespace Hermes
{
  namespace Algebra
  {
    /** Type for storing number in Mumps structures */
    template <typename Scalar> struct mumps_type;

    /** Type for storing number in Mumps complex structures */
    template <>
    struct mumps_type<std::complex<double> >
    {
    /** Type for storing mumps struct in Mumps complex structures */
      typedef ZMUMPS_STRUC_C mumps_struct;
    /** Type for storing scalar number in Mumps complex structures */
      typedef ZMUMPS_COMPLEX mumps_Scalar;
    };

    /** Type for storing number in Mumps real structures */
    template <>
    struct mumps_type<double>
    {
    /** Type for storing mumps struct in Mumps real structures */
      typedef DMUMPS_STRUC_C mumps_struct;
    /** Type for storing scalar number in Mumps real structures */
      typedef double mumps_Scalar;
    };

    /** \brief Matrix used with MUMPS solver */
    template <typename Scalar>
    class MumpsMatrix : public SparseMatrix<Scalar>
    {
    public:
      MumpsMatrix();
      virtual ~MumpsMatrix();

      virtual void alloc();
      virtual void free();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;
      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(MumpsMatrix* mat);
      /// Add matrix to diagonal.
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_to_diagonal_blocks(int num_stages, MumpsMatrix* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, MumpsMatrix* mat);

      /// Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      /// Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);
      /// Creates matrix using size, nnz, and the three arrays.
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);
      /// Duplicates a matrix (including allocation).
      MumpsMatrix* duplicate();

    protected:
      /// MUMPS specific data structures for storing the system matrix (CSC format).
      unsigned int nnz;          ///< Number of non-zero elements.
      int *irn;         ///< Row indices.
      int *jcn;         ///< Column indices.
      typename mumps_type<Scalar>::mumps_Scalar *Ax; ///< Matrix entries (column-wise).
      int *Ai;          ///< Row indices of values in Ax.
      unsigned int *Ap;          ///< Index to Ax/Ai, where each column starts.

      friend class Solvers::MumpsSolver<Scalar>;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /** \brief Vector used with MUMPS solver */
    template <typename Scalar>
    class MumpsVector : public Vector<Scalar>
    {
    public:
      MumpsVector();
      virtual ~MumpsVector();

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
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

    protected:
      // MUMPS specific data structures for storing the rhs.
      ///Vector data.
      Scalar *v;

      friend class Solvers::MumpsSolver<Scalar>;
    };
  }
  namespace Solvers
  {
    /// Encapsulation of MUMPS linear solver.
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API MumpsSolver : public DirectSolver<Scalar>
    {
    public:
      /// Constructor of MumpsSolver.
      /// @param[in] m matrix pointer
      /// @param[in] rhs right hand side pointer
      MumpsSolver(MumpsMatrix<Scalar> *m, MumpsVector<Scalar> *rhs);
      virtual ~MumpsSolver();

      virtual bool solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      MumpsMatrix<Scalar> *m;
      /// Right hand side.
      MumpsVector<Scalar> *rhs;

      /// \todo document
      bool setup_factorization();

      /// MUMPS specific structure with solver parameters.
      typename mumps_type<Scalar>::mumps_struct param;

      /// @return false on error
      bool check_status();

      /// (Re)initialize a MUMPS instance.
      /// @return true on succes
      /// \sa #check_status()
      bool reinit();
      /// True if solver is inited.
      bool inited;
    private:
      void mumps_c(typename mumps_type<Scalar>::mumps_struct * param);  //wrapper around dmums_c or zmumps_c
    };
  }
}
#endif
#endif