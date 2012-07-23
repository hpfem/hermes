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
    /// @ingroup userSolvingAPI
    /// Class for Newton's method.
    template<typename Scalar>
    class HERMES_API NewtonSolver : public NonlinearSolver<Scalar>, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>
    {
    public:
      NewtonSolver(DiscreteProblem<Scalar>* dp);
      NewtonSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space);
      NewtonSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces);
      void init_linear_solver();

      ~NewtonSolver();

      /// Solve with user-defined tolerances.
      /// \param[in] residual_as_function Translate the residual vector into a residual function (or multiple functions)
      ///                                 in the corresponding finite element space(s) and measure their norm(s) there.
      ///                                 This is more meaningful than just measuring the l2-norm of the residual vector,
      ///                                 since in the FE space not all components in the residual vector have the same weight.
      ///                                 On the other hand, this is slower as it requires global norm calculation, and thus
      ///                                 numerical integration over the entire domain. Therefore this option is off by default.
      void solve(Scalar* coeff_vec = NULL);

      /// A solve() method where the jacobian is reused.
      /// Version with user-defined tolerances.
      void solve_keep_jacobian(Scalar* coeff_vec = NULL);

      /// Sets the maximum allowed norm of the residual during the calculation.
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Sets minimum damping coefficient.
      void set_min_allowed_damping_coeff(double min_allowed_damping_coeff);

      /// Call NonlinearSolver::set_iterative_method() and set the method to the linear solver (if applicable).
      virtual void set_iterative_method(const char* iterative_method_name);

      /// Call NonlinearSolver::set_preconditioner() and set the method to the linear solver (if applicable).
      virtual void set_preconditioner(const char* preconditioner_name);

      /// Get times accumulated by this instance of NewtonSolver.
      double get_setup_time() const { return setup_time; }
      double get_assemble_time() const { return assemble_time; }
      double get_solve_time() const { return solve_time; }

      void set_newton_tol(double newton_tol);
      void set_newton_max_iter(int newton_max_iter);
      void set_residual_as_function();

      /// set time information for time-dependent problems.
      virtual void setTime(double time);
      virtual void setTimeStep(double timeStep);

      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;

    protected:
      /// Jacobian.
      SparseMatrix<Scalar>* jacobian;

      /// Residual.
      Vector<Scalar>* residual;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* linear_solver;

      /// Used by method solve_keep_jacobian().
      SparseMatrix<Scalar>* kept_jacobian;

      double newton_tol;
      int newton_max_iter;
      bool residual_as_function;

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      static double max_allowed_residual_norm;

      static double min_allowed_damping_coeff;

      double currentDampingCofficient;

      /// Times spent in individual phases of the computation.
      double setup_time;
      double assemble_time;
      double solve_time;
      
      /// This instance owns its DP.
      const bool own_dp;
    };
  }
}
#endif