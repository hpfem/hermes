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
    class HERMES_API NewtonSolver : public NonlinearSolver<Scalar>, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, public Hermes::Mixins::OutputAttachable, public Hermes::Hermes2D::Mixins::MatrixRhsOutput<Scalar>, public Hermes::Hermes2D::Mixins::StateQueryable
    {
    public:
      NewtonSolver();
      NewtonSolver(DiscreteProblem<Scalar>* dp);
      NewtonSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space);
      NewtonSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces);
      void init_linear_solver();

      ~NewtonSolver();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NewtonSolver"; }

      /// Solve.
      /// \param[in] coeff_vec Ceofficient vector to start from.
      void solve(Scalar* coeff_vec = NULL);

      /// Solve.
      /// \param[in] initial_guess Solution to start from (which is projected to obtain the initial coefficient vector.
      void solve(Solution<Scalar>* initial_guess);

      /// Solve.
      /// \param[in] initial_guess Solutions to start from (which is projected to obtain the initial coefficient vector.
      void solve(Hermes::vector<Solution<Scalar>*> initial_guess);

      /// Solve which keeps jacobian.
      /// A solve() method where the jacobian is reused.
      void solve_keep_jacobian(Scalar* coeff_vec = NULL);
      
      /// Solve which keeps jacobian.
      /// \param[in] initial_guess Solution to start from (which is projected to obtain the initial coefficient vector.
      void solve_keep_jacobian(Solution<Scalar>* initial_guess);
      
      /// Solve which keeps jacobian.
      /// \param[in] initial_guess Solutions to start from (which is projected to obtain the initial coefficient vector.
      void solve_keep_jacobian(Hermes::vector<Solution<Scalar>*> initial_guess);

      /// Sets the maximum allowed norm of the residual during the calculation.
      /// Default: 1E9
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Sets minimum damping coefficient.
      /// Default: 1E-4
      void set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set);

      /// Call NonlinearSolver::set_iterative_method() and set the method to the linear solver (if applicable).
      virtual void set_iterative_method(const char* iterative_method_name);

      /// Call NonlinearSolver::set_preconditioner() and set the method to the linear solver (if applicable).
      virtual void set_preconditioner(const char* preconditioner_name);

      /// Interpret the residual as a function.
      /// Translate the residual vector into a residual function (or multiple functions)
      /// in the corresponding finite element space(s) and measure their norm(s) there.
      /// This is more meaningful than just measuring the l2-norm of the residual vector,
      /// since in the FE space not all components in the residual vector have the same weight.
      /// On the other hand, this is slower as it requires global norm calculation, and thus
      /// numerical integration over the entire domain. Therefore this option is off by default.
      void set_residual_as_function();

      /// Set the residual norm tolerance for ending the Newton's loop.
      /// Default: 1E-8.
      void set_newton_tol(double newton_tol);

      /// Set the maximum number of Newton's iterations.
      /// Default: 15
      void set_newton_max_iter(int newton_max_iter);

      /// Set time information for time-dependent problems.
      /// See the class Hermes::Mixins::TimeMeasurable.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// See the class Hermes::Hermes2D::Mixins::SettableSpaces.
      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;

      /// Turn on or off manual damping (default is the automatic) and optionally sets manual damping coefficient.
      /// Default: default is the automatic damping, default coefficient if manual damping used is set by this method.
      /// \param[in] onOff on(true)-manual damping, off(false)-automatic damping.
      /// \param[in] coeff The (perpetual) damping coefficient in the case of manual damping. Ignored in the case of automatic damping.
      void set_manual_damping_coeff(bool onOff, double coeff = 1.0);
      
      /// Make the automatic damping start with this coefficient.
      /// This will also be the top bound for the coefficient.
      /// Default: 1.0
      /// \param[in] coeff The initial damping coefficient. Must be > 0 and <= 1.0.
      void set_initial_auto_damping_coeff(double coeff);
      
      /// Set the ratio to the automatic damping.
      /// When the damping coefficient is decided to be descreased or increased, this is the ratio
      /// how it will be changed (this is the bigger ( > 1.0 ) of the two possible values).
      /// I.e. when the damping coefficient is shortened 3 times if deemed too big, make the parameter not 0.333333, but 3.0.
      /// Default: 2.0
      /// \param[in] ratio The ratio (again, it must be > 1.0, and it represents the inverse of the shortening factor).
      void set_auto_damping_ratio(double ratio);

      /// Set the ratio of the current residual norm and the previous residual norm necessary to deem a step 'successful'.
      /// It can be either > 1.0, meaning that even if the norm increased, the step will be 'successful', or < 1.0, meaning
      /// that even though the residual norm goes down, we will further decrease the damping coefficient.
      /// Default: 0.95
      /// param[in] ratio The ratio, must be positive.
      void set_sufficient_improvement_factor(double ratio);

      /// Set how many successful steps are necessary for the damping coefficient to be increased, by multiplication by the parameter
      /// set by set_auto_damping_ratio().
      /// The coefficient is then increased after each 'successful' step, if the sequence of such is not interrupted by an 'unsuccessful' step.
      /// Default: 1
      /// \param[in] steps Number of steps.
      void set_necessary_successful_steps_to_increase(unsigned int steps);

      /// Set the weak forms.
      void set_weak_formulation(const WeakForm<Scalar>* wf);

    protected:
      /// This instance owns its DP.
      const bool own_dp;

      /// Used by method solve_keep_jacobian().
      SparseMatrix<Scalar>* kept_jacobian;

      /// Internal setting of default values (see individual set methods).
      void init_attributes();

      /// Jacobian.
      SparseMatrix<Scalar>* jacobian;

      /// Residual.
      Vector<Scalar>* residual;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* linear_solver;

      double newton_tol;
      int newton_max_iter;
      bool residual_as_function;

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      double max_allowed_residual_norm;
      double min_allowed_damping_coeff;

      double currentDampingCofficient;
      
      /// Manual / auto.
      bool manual_damping;
      /// The ratio between two damping coeffs when changing.
      double auto_damping_ratio;
      /// The initial (and maximum) damping coefficient
      double initial_auto_damping_ratio;
      /// Sufficient improvement for continuing.
      double sufficient_improvement_factor;
      /// necessary number of steps to increase back the damping coeff.
      unsigned int necessary_successful_steps_to_increase;
    };
  }
}
#endif
