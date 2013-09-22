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
#include "nonlinear_solver.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup userSolvingAPI
    /// Class for Newton's method.<br>
    /// Typical usage:<br>
    /// // Initialize Newton's solver.<br>
    /// // Here wf is Hermes2D::WeakForm<double>, space is Hermes2D::Space<double><br>
    /// Hermes::Hermes2D::NewtonSolver<double> newton_solver(&wf, &space);<br>
    /// Set a whole bunch of parameters according to your liking.<br>
    /// See the class documentation for all possible parameters.<br>
    /// newton_solver.set_tolerance(1e-6);<br>
    /// newton_solver.set_max_allowed_iterations(15);<br>
    /// newton_solver.set_max_allowed_residual_norm(1e6);<br>
    /// newton_solver.set_min_allowed_damping_coeff(1e-3);<br>
    /// <br>
    /// // Solve the linear problem.<br>
    /// try<br>
    /// {<br>
    ///&nbsp;// Just call solve().<br>
    ///&nbsp;newton_solver.solve();<br>
    /// <br>
    ///&nbsp;// Get the solution vector from the solver.<br>
    ///&nbsp;double* sln_vector = newton_solver.get_sln_vector();<br>
    /// <br>
    ///&nbsp;// Translate the solution vector into the previously initialized Solution<double> using the static method vector_to_solution.<br>
    ///&nbsp;Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, &space, &sln);<br>
    /// }<br>
    /// // All kinds of Exceptions may happen (Linear algebraic solver, some bad parameters, some data not initialized...)<br>
    /// catch(Hermes::Exceptions::Exception& e)<br>
    /// {<br>
    ///&nbsp;e.print_msg();<br>
    ///&nbsp;return -1;<br>
    /// }<br>
    /// // For illustrative purposes, otherwise one can just catch std::exception directly, as Hermes::Exceptions::Exception derive from it.<br>
    /// catch(std::exception& e)<br>
    /// {<br>
    ///&nbsp;std::cout << e.what(); <br>
    ///&nbsp;return -1;<br>
    /// }<br>
    template<typename Scalar>
    class HERMES_API NewtonSolver : public Hermes::Hermes2D::NonlinearSolver<Scalar>
    {
    public:
      NewtonSolver();
      NewtonSolver(DiscreteProblem<Scalar>* dp);
      NewtonSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      NewtonSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual ~NewtonSolver();

      // See the base class for details, the following serves only for avoiding C++ name-hiding.
      using Solver<Scalar>::solve;
      /// Solve.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual void solve(Scalar* coeff_vec);

#pragma region damping-public
      /// Sets minimum damping coefficient.
      /// Default: 1E-4
      void set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set);

      /// Turn on or off manual damping (default is the automatic) and optionally sets manual damping coefficient.
      /// Default: default is the automatic damping, default coefficient if manual damping used is set by this method.
      /// \param[in] onOff on(true)-manual damping, off(false)-automatic damping.
      /// \param[in] coeff The (perpetual) damping coefficient in the case of manual damping. Ignored in the case of automatic damping.
      void set_manual_damping_coeff(bool onOff, double coeff);
      
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
#pragma endregion

#pragma region jacobian_recalculation-public
      /// Set the ratio of the current residual norm and the previous residual norm necessary to deem a step 'successful'.
      /// IMPORTANT: it is truly a FACTOR, i.e. the two successive residual norms are put in a fraction and this number is
      /// then compared to the ratio set by this method.
      void set_sufficient_improvement_factor_jacobian(double ratio);

      /// Set maximum number of steps (Newton iterations) that a jacobian can be reused if it is deemed a 'successful' reusal
      /// with respect to the improvement factor.
      void set_max_steps_with_reused_jacobian(unsigned int steps);
#pragma endregion
      
    protected:
      /// Common constructors code.
      /// Internal setting of default values (see individual set methods).
      void init_newton();

      /// \return Whether or not should the processing continue.
      virtual void on_damping_factor_updated();
      /// \return Whether or not should the processing continue.
      virtual void on_reused_jacobian_step_begin();
      /// \return Whether or not should the processing continue.
      virtual void on_reused_jacobian_step_end();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NewtonSolver"; }

      void init_solving(Scalar*& coeff_vec);

      bool do_initial_step_return_finished(Scalar* coeff_vec);
      void assemble_residual(Scalar* coeff_vec);
      void solve_linear_system(Scalar* coeff_vec);

      void finalize_solving(Scalar* coeff_vec);
      void deinit_solving(Scalar* coeff_vec);
      
      /// Act upon the convergence state.
      /// \return If the main loop in solve() should finalize after this.
      bool handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec);

      /// Calculates the new damping coefficient.
      bool calculate_damping_factor(unsigned int& successful_steps);

      /// Shortcut method for getting the current iteration.
      int get_current_iteration_number();

      /// Output info about the step.
      void step_info();

#pragma region damping-private
      /// Manual / auto.
      bool manual_damping;

      /// Manual.
      double manual_damping_factor;

      /// Auto.
      /// The ratio between two damping coeffs when changing.
      double auto_damping_ratio;
      /// The initial (and maximum) damping coefficient
      double initial_auto_damping_factor;
      /// Sufficient improvement for continuing.
      double sufficient_improvement_factor;
      /// necessary number of steps to increase back the damping coeff.
      unsigned int necessary_successful_steps_to_increase;
      /// Minimum allowed damping coeff.
      double min_allowed_damping_coeff;      
#pragma endregion

#pragma region jacobian_recalculation-private
      /// For deciding if the jacobian is reused at this point.
      bool force_reuse_jacobian_values(unsigned int& successful_steps_with_reused_jacobian);
      /// For deciding if the reused jacobian did not bring residual increase at this point.
      bool jacobian_reused_okay(unsigned int& successful_steps_with_reused_jacobian);

      double sufficient_improvement_factor_jacobian;
      unsigned int max_steps_with_reused_jacobian;
      
      /// Backup vector for unsuccessful reuse of Jacobian.
      Vector<Scalar>* residual_back;
#pragma endregion

#pragma region OutputAttachable
      // For derived classes - read-only access.
      const OutputParameterUnsignedInt& successful_steps_damping() const { return this->p_successful_steps_damping; };
      const OutputParameterUnsignedInt& successful_steps_jacobian() const { return this->p_successful_steps_jacobian; };
      
      const OutputParameterDoubleVector& damping_factors() const { return this->p_damping_factors; };
      const OutputParameterBool& residual_norm_drop() const { return this->p_residual_norm_drop; };
      const OutputParameterBoolVector& iterations_with_recalculated_jacobian() const { return this->p_iterations_with_recalculated_jacobian; };

    private:
      // Parameters for OutputAttachable mixin.
      OutputParameterBoolVector p_iterations_with_recalculated_jacobian;
      OutputParameterUnsignedInt p_successful_steps_damping;
      OutputParameterUnsignedInt p_successful_steps_jacobian;
      OutputParameterDoubleVector p_damping_factors;
      OutputParameterBool p_residual_norm_drop;
#pragma endregion

    private:
			Scalar* coeff_vec_back;

      friend class NonlinearConvergenceMeasurement<Scalar>;
    };
  }
}
#endif
