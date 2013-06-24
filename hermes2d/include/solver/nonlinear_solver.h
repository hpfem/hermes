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
/*! \file nonlinear_solver.h
\brief General nonlinear solver functionality.
*/
#ifndef __H2D_NONLINEAR_SOLVER_H_
#define __H2D_NONLINEAR_SOLVER_H_

#include "hermes_common.h"
#include "solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// \brief Base class for defining interface for nonlinear solvers.
    ///
    template <typename Scalar>
    class NonlinearSolver : public Solver<Scalar>
    {
    public:
      NonlinearSolver();
      NonlinearSolver(DiscreteProblem<Scalar>* dp);
      NonlinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      NonlinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual ~NonlinearSolver();

      /// Set the name of the iterative method employed by AztecOO (ignored
      /// by the other solvers).
      /// \param[in] preconditioner_name See the attribute preconditioner.
      void set_iterative_method(const char* iterative_method_name);

      /// Set the name of the preconditioner employed by AztecOO (ignored by
      /// the other solvers).
      /// \param[in] preconditioner_name See the attribute preconditioner.
      void set_preconditioner(const char* preconditioner_name);
      
      /// Set the preconditioner object of Epetra_Operator type employed by AztecOO (ignored by
      /// the other solvers).
      /// \param[in] pc Pointer to a wrapper around Epetra_Operator based preconditioners (see precond.h).
      void set_preconditioner(Hermes::Preconditioners::Precond<Scalar> *pc);

      /// Set the maximum number of iterations, thus co-determine when to stop iterations.
      void set_max_allowed_iterations(int max_allowed_iterations);
    
    protected:
      /// Maximum number of iterations allowed.
      int max_allowed_iterations;

      /// There was no initial coefficient vector passed, so this instance had to create one
      /// and this serves as the identificator according to which it will be deleted.
      bool delete_coeff_vec;

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