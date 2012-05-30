// This file is part of Hermes2D
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
/*! \file solver_newton.h
\brief Newton's method.
*/
#ifndef __H2D_SOLVER_NEWTON_H_
#define __H2D_SOLVER_NEWTON_H_

#include "global.h"
#include "discrete_problem.h"
#include "exceptions.h"


namespace Hermes
{
  namespace Hermes2D
  {
    /// Class for Newton's method.
    template<typename Scalar>
    class HERMES_API NewtonSolver : public NonlinearSolver<Scalar>
    {
    public:
      NewtonSolver(DiscreteProblem<Scalar>* dp);
      NewtonSolver(DiscreteProblem<Scalar>* dp, Hermes::MatrixSolverType matrix_solver_type);
      void init_linear_solver();

      ~NewtonSolver();

      /// Solve with user-defined tolerances.
      /// \param[in] residual_as_function Translate the residual vector into a residual function (or multiple functions)
      ///                                 in the corresponding finite element space(s) and measure their norm(s) there.
      ///                                 This is more meaningful than just measuring the l2-norm of the residual vector,
      ///                                 since in the FE space not all components in the residual vector have the same weight.
      ///                                 On the other hand, this is slower as it requires global norm calculation, and thus
      ///                                 numerical integration over the entire domain. Therefore this option is off by default.
      void solve(Scalar* coeff_vec = NULL, double newton_tol = 1e-8, 
                 int newton_max_iter = 100, bool residual_as_function = false);

      /// A solve() method where the jacobian is reused.
      /// Version with user-defined tolerances.
      void solve_keep_jacobian(Scalar* coeff_vec = NULL, double newton_tol = 1e-8, 
                               int newton_max_iter = 100, bool residual_as_function = false);

      /// Sets the maximum allowed norm of the residual during the calculation.
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Call NonlinearSolver::set_iterative_method() and set the method to the linear solver (if applicable).
      virtual void set_iterative_method(const char* iterative_method_name);

      /// Call NonlinearSolver::set_preconditioner() and set the method to the linear solver (if applicable).
      virtual void set_preconditioner(const char* preconditioner_name);

      /// Get times accumulated by this instance of NewtonSolver.
      double get_setup_time() const { return setup_time; }
      double get_assemble_time() const { return assemble_time; }
      double get_solve_time() const { return solve_time; }

      /// Attach an external timer to which this instance of NewtonSolver will accumulate time spent in it.
      void attach_timer(TimePeriod *timer) { this->timer = timer; reset_times(); }

      /// Reset times to zero.
      void reset_times() { setup_time = assemble_time = solve_time = 0.; }

    protected:
      /// Jacobian.
      SparseMatrix<Scalar>* jacobian;

      /// Residual.
      Vector<Scalar>* residual;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* linear_solver;

      /// Used by method solve_keep_jacobian().
      SparseMatrix<Scalar>* kept_jacobian;

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      static double max_allowed_residual_norm;

      /// Times spent in individual phases of the computation.
      double setup_time;
      double assemble_time;
      double solve_time;

      /// Pointer to an external timer to which this instance of NewtonSolver accumulates time spent in it.
      TimePeriod *timer;
    };
  }
}
#endif

