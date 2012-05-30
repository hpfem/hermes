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
/*! \file amesos_solver.h
\brief AmesosSolver class as an interface to Amesos.
*/
#ifndef __HERMES_COMMON_AMESOS_SOLVER_H_
#define __HERMES_COMMON_AMESOS_SOLVER_H_
#include "../config.h"
#ifdef HAVE_AMESOS
#include "linear_matrix_solver.h"
#include "epetra.h"
#include "compat.h"
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

namespace Hermes
{
  namespace Solvers
  {
    /// \brief Encapsulation of Amesos linear solver.
    ///
    /// @ingroup solvers
    template<typename Scalar>
    class HERMES_API AmesosSolver : public DirectSolver<Scalar>
    {
    public:
      static bool is_available(const char *name);

      /// Returns the current UseTranspose setting.
      bool use_transpose();
      /// If set true, X will be set to the solution of A^T X = B (not A X = B).
      void set_use_transpose(bool use_transpose);

      AmesosSolver(const char *solver_type, EpetraMatrix<Scalar> *m, EpetraVector<Scalar> *rhs);

      virtual ~AmesosSolver();
      virtual bool solve();
      virtual int get_matrix_size();
    protected:
      static Amesos factory;
      Amesos_BaseSolver *solver;
      Epetra_LinearProblem problem;
      EpetraMatrix<Scalar> *m;
      EpetraVector<Scalar> *rhs;

      bool setup_factorization();
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Hermes::MatrixSolverType matrix_solver_type, Matrix<T>* matrix, Vector<T>* rhs);
    };
  }
}
#endif
#endif
