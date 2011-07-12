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
/*! \file solver.h
\brief General linear/nonlinear/iterative solver functionality.
*/
#ifndef __HERMES_COMMON_NONLINEAR_SOLVER_H_
#define __HERMES_COMMON_NONLINEAR_SOLVER_H_

#include "precond.h"
#include "dpinterface.h"

namespace Hermes
{
  namespace Solvers
  {
    /// \brief Base class for defining interface for nonlinear solvers.
    ///
    template <typename Scalar>
    class NonlinearSolver
    {
    public:
      NonlinearSolver(DiscreteProblemInterface<Scalar>* dp);
      NonlinearSolver(DiscreteProblemInterface<Scalar>* dp, Hermes::MatrixSolverType matrix_solver_type);

      ~NonlinearSolver();

      /// Basic solve method.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual bool solve(Scalar* coeff_vec) = 0;

      Scalar *get_sln_vector();

      double get_time();

      /// Sets the attribute verbose_output to the paramater passed.
      void set_verbose_output(bool verbose_output_to_set);
      
      /// Set the name of the iterative method employed by AztecOO (ignored
      /// by the other solvers). 
      /// \param[in] preconditioner_name See the attribute preconditioner.
      void set_iterative_method(const char* iterative_method_name);

      /// Set the name of the preconditioner employed by AztecOO (ignored by
      /// the other solvers).
      /// \param[in] preconditioner_name See the attribute preconditioner.
      void set_preconditioner(const char* preconditioner_name);

    protected:
      DiscreteProblemInterface<Scalar>* dp; ///< FE problem being solved (not NULL in case of using
      ///< NonlinearProblem(DiscreteProblemInterface *) ctor.

      /// The solution vector.
      Scalar* sln_vector;

      /// For use of error measurement.
      int error;

      double time;  ///< time spent on solving (in secs)

      /// Linear solver to use, choices: 
      /// SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
      /// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
      /// Default: SOLVER_UMFPACK.
      Hermes::MatrixSolverType matrix_solver_type;

      /// Verbose output.
      /// Set to 'true' by default.
      bool verbose_output;

      /// Preconditioned solver.
      bool precond_yes;

      /// Name of the iterative method employed by AztecOO (ignored
      /// by the other solvers). 
      /// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
      char* iterative_method;
      
      /// Name of the preconditioner employed by AztecOO (ignored by
      /// the other solvers).
      /// Possibilities: none, jacobi, neumann, least-squares, or a
      ///  preconditioner from IFPACK (see solver/aztecoo.h).
      char* preconditioner;
    };
  }
}
#endif
