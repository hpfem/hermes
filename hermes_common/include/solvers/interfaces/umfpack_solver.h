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
#include "algebra/cs_matrix.h"

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
      UMFPackLinearMatrixSolver(CSCMatrix<Scalar> *m, SimpleVector<Scalar> *rhs);
      virtual ~UMFPackLinearMatrixSolver();
      virtual void solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      CSCMatrix<Scalar> *m;
      /// Right hand side vector.
      SimpleVector<Scalar> *rhs;

      /// \brief Reusable factorization information (A denotes matrix represented by the pointer 'm').
      /// Reordering of matrix A to reduce fill-in during factorization.
      void *symbolic;
      void *numeric;  ///< LU factorization of matrix A.

      /// \todo document
      void free_factorization_data();
      /// \todo document
      bool setup_factorization();
      template <typename T> friend class Hermes::Algebra::CSCMatrix;
      template <typename T> friend class Hermes::Algebra::SimpleVector;
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs, bool use_direct_solver);
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