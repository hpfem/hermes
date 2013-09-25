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
#include "nonlinear_convergence_measurement.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Nonlinear Convergence state.
    enum NonlinearConvergenceState
    {
      Converged,
      NotConverged,
      BelowMinDampingCoeff,
      AboveMaxAllowedResidualNorm,
      AboveMaxIterations,
      Error
    };

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

      /// Clear (reset) all tolerances.
      virtual void clear_tolerances();

      /// Sets the maximum allowed norm of the residual during the calculation.
      /// Default: 1E9
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Interpret the residual as a function.
      /// Translate the residual vector into a residual function (or multiple functions)
      /// in the corresponding finite element space(s) and measure their norm(s) there.
      /// This is more meaningful than just measuring the l2-norm of the residual vector,
      /// since in the FE space not all components in the residual vector have the same weight.
      /// On the other hand, this is slower as it requires global norm calculation, and thus
      /// numerical integration over the entire domain. Therefore this option is off by default.
      void set_residual_as_function();

      /// Find out the convergence state.
      NonlinearConvergenceState get_convergence_state();

      /// Act upon the convergence state.
      /// \return If the main loop in solve() should finalize after this.
      virtual bool handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec) = 0;

      /// Set the residual norm tolerance for ending the Newton's loop.
      /// Default: this->set_tolerance(1e-8, ResidualNormAbsolute);
      /// \param[in] handleMultipleTolerancesAnd If true, multiple tolerances defined will have to be all fulfilled in order to proclaim
      /// solution as a correct one. If false, only one will be enough.
      void set_tolerance(double newton_tol, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd = false);

      /// Get the number of iterations.
      int get_num_iters() const;

    protected:
      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NonlinearSolver"; }

      /// Shared code for constructors.
      void init_nonlinear();

      /// Shortcut method for getting the current iteration.
      virtual int get_current_iteration_number() = 0;

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      double max_allowed_residual_norm;

      bool residual_as_function;

      /// Calculates the residual norm.
      double calculate_residual_norm();

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

      /// Tolerances for all NonlinearConvergenceMeasurementType numbered sequentially as the enum NonlinearConvergenceMeasurementType is.
      double tolerance[NonlinearConvergenceMeasurementTypeCount];

      /// info about set tolerances.
      bool tolerance_set[NonlinearConvergenceMeasurementTypeCount];

      /// If true, multiple tolerances defined will have to be all fulfilled in order to proclaim
      /// solution as a correct one. If false, only one will be enough.
      bool handleMultipleTolerancesAnd;

      int num_iters;
#pragma region OutputAttachable
      // For derived classes - read-only access.
      const OutputParameterDoubleVector& residual_norms() const { return this->p_residual_norms; };
      const OutputParameterDoubleVector& solution_norms() const { return this->p_solution_norms; };
      const OutputParameterDoubleVector& solution_change_norms() const { return this->p_solution_change_norms; };

      // Parameters for OutputAttachable mixin.
      OutputParameterDoubleVector p_residual_norms;
      OutputParameterDoubleVector p_solution_norms;
      OutputParameterDoubleVector p_solution_change_norms;
#pragma endregion
      friend class NonlinearConvergenceMeasurement<Scalar>;
    };
  }

  namespace Exceptions
  {
    class HERMES_API NonlinearException : public Hermes::Exceptions::Exception
    {
    public:
      NonlinearException(Hermes2D::NonlinearConvergenceState convergenceState);

      Hermes2D::NonlinearConvergenceState get_exception_state();

    protected:
      Hermes2D::NonlinearConvergenceState convergenceState;
    };
  }
}
#endif