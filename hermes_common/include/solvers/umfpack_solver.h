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
#include "cs_matrix.h"

extern "C"
{
#include <umfpack.h>
}

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class HERMES_API UMFPackLinearMatrixSolver;
  }

  namespace Algebra
  {
    /// \brief This class is to be used with UMFPack solver only.
    template <typename Scalar>
    class HERMES_API UMFPackMatrix : public CSCMatrix<Scalar>
    {
      // Friends.
      template <typename T> friend class Hermes::Solvers::UMFPackLinearMatrixSolver;
      template <typename T> friend class Hermes::Solvers::CSCIterator;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /// \brief Class representing the vector for UMFPACK.
    template <typename Scalar>
    class HERMES_API UMFPackVector : public Vector<Scalar>
    {
    public:
      /// Default constructor.
      UMFPackVector();
      /// Constructor of vector with specific size.
      /// @param[in] size size of vector
      UMFPackVector(unsigned int size);
      virtual ~UMFPackVector();
      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual Scalar get(unsigned int idx) const;
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual void change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual void set_vector(Vector<Scalar>* vec);
      virtual void set_vector(Scalar* vec);
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
      template <typename T> friend class Hermes::Solvers::CSCIterator;
      template<typename T> friend Vector<T>* Hermes::Algebra::create_vector(bool);
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
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs, bool);
      char* check_status(const char *fn_name, int status);

      /// For output (the array Info is filled regardless of the value of the level).
      void set_output_level(double level);

      /// For reporting.
      double Info [UMFPACK_INFO];
      double Control [UMFPACK_CONTROL];
    };
  }
}
#endif
#endif