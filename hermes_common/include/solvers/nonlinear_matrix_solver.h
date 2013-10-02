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
/*! \file nonlinear_matrix_solver.h
\brief General nonlinear solver functionality.
*/
#ifndef __HERMES_COMMON_NONLINEAR_MATRIX_SOLVER_H_
#define __HERMES_COMMON_NONLINEAR_MATRIX_SOLVER_H_

#include "precond.h"
#include "exceptions.h"
#include "algebra/cs_matrix.h"
#include "algebra/vector.h"
#include "mixins.h"
#include "algebra/algebra_mixins.h"
#include "solvers/nonlinear_convergence_measurement.h"
#include "solvers/matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    /// Convergence measurement strategies.
    /// Count of the types - for solver to hold arrays of such a length.
    const int NonlinearConvergenceMeasurementTypeCount = 7;

    /// This specifies the quantity that is compared to newton_tolerance (settable by set_tolerance()).
    enum NonlinearConvergenceMeasurementType
    {
      ResidualNormRelativeToInitial = 0x0001,
      ResidualNormRelativeToPrevious = 0x0002,
      ResidualNormRatioToInitial = 0x0004,
      ResidualNormRatioToPrevious = 0x0008,
      ResidualNormAbsolute = 0x0010,
      SolutionChangeAbsolute = 0x0020,
      SolutionChangeRelative = 0x0040
    };

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

    template<typename Scalar> class HERMES_API NonlinearConvergenceMeasurement;
    
    /// \brief Base class for defining interface for nonlinear solvers.
    ///
    template <typename Scalar>
    class NonlinearMatrixSolver : 
      public virtual Hermes::Solvers::MatrixSolver<Scalar>,
      public virtual Hermes::Mixins::OutputAttachable,
      public virtual Hermes::Mixins::StateQueryable, 
      public virtual Hermes::Mixins::TimeMeasurable
    {
    public:
      NonlinearMatrixSolver();
      virtual ~NonlinearMatrixSolver();

      /// Solve.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual void solve(Scalar* coeff_vec);

      /// Set the maximum number of iterations, thus co-determine when to stop iterations.
      void set_max_allowed_iterations(int max_allowed_iterations);

      /// Clear (reset) all tolerances.
      virtual void clear_tolerances();

      /// Sets the maximum allowed norm of the residual during the calculation.
      /// Default: 1E9
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Set the residual norm tolerance for ending the Newton's loop.
      /// Default: this->set_tolerance(1e-8, ResidualNormAbsolute);
      /// \param[in] handleMultipleTolerancesAnd If true, multiple tolerances defined will have to be all fulfilled in order to proclaim
      /// solution as a correct one. If false, only one will be enough.
      void set_tolerance(double newton_tol, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd = false);
      
      /// Get the number of iterations.
      int get_num_iters() const;
      
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

      virtual void assemble_residual() = 0;
      virtual void assemble_jacobian() = 0;
      virtual void assemble() = 0;
      
      /// \return Whether or not should the processing continue.
      virtual void on_damping_factor_updated();
      /// \return Whether or not should the processing continue.
      virtual void on_reused_jacobian_step_begin();
      /// \return Whether or not should the processing continue.
      virtual void on_reused_jacobian_step_end();

      /// Act upon the convergence state.
      /// \return If the main loop in solve() should finalize after this.
      virtual bool handle_convergence_state_return_finished(NonlinearConvergenceState state);

      /// Initialization - called at the beginning of solving.
      /// Very important e.g. for assigning DOFs before assembling.
      virtual void init_solving(Scalar* coeff_vec);

      /// Find out the convergence state.
      virtual NonlinearConvergenceState get_convergence_state();

      /// Initial step.
      bool do_initial_step_return_finished();
      
      /// Solve the step's linear system.
      virtual void solve_linear_system();

      /// Update the solution.
      /// This is a method that serves the purpose of distinguishing methods that solve for increment (Newton), or for solution (Picard).
      virtual double update_solution_return_change_norm(Scalar* linear_system_solution) = 0;

      /// Internal.
      void finalize_solving();

      /// Internal.
      virtual void deinit_solving();
      
      /// Calculates the new damping coefficient.
      bool calculate_damping_factor(unsigned int& successful_steps);

      /// Shortcut method for getting the current iteration.
      int get_current_iteration_number();

      /// Output info about the step.
      void step_info();
      
      /// Norm for convergence.
      double calculate_residual_norm();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NonlinearMatrixSolver"; }

      /// Shared code for constructors.
      void init_nonlinear();

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      double max_allowed_residual_norm;

      /// Maximum number of iterations allowed.
      int max_allowed_iterations;

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
      const OutputParameterUnsignedInt& iteration() const { return this->p_iteration; };
      const OutputParameterDoubleVector& residual_norms() const { return this->p_residual_norms; };
      const OutputParameterDoubleVector& solution_norms() const { return this->p_solution_norms; };
      const OutputParameterDoubleVector& solution_change_norms() const { return this->p_solution_change_norms; };
      const OutputParameterUnsignedInt& successful_steps_damping() const { return this->p_successful_steps_damping; };
      const OutputParameterUnsignedInt& successful_steps_jacobian() const { return this->p_successful_steps_jacobian; };
      const OutputParameterDoubleVector& damping_factors() const { return this->p_damping_factors; };
      const OutputParameterBool& residual_norm_drop() const { return this->p_residual_norm_drop; };
      const OutputParameterBoolVector& iterations_with_recalculated_jacobian() const { return this->p_iterations_with_recalculated_jacobian; };
    
      /// Parameters for OutputAttachable mixin.
      /// Should be private, but then it does not work.
      OutputParameterDoubleVector p_residual_norms;
      OutputParameterDoubleVector p_solution_norms;
      OutputParameterDoubleVector p_solution_change_norms;
      OutputParameterBoolVector p_iterations_with_recalculated_jacobian;
      OutputParameterUnsignedInt p_successful_steps_damping;
      OutputParameterUnsignedInt p_successful_steps_jacobian;
      OutputParameterDoubleVector p_damping_factors;
      OutputParameterBool p_residual_norm_drop;
      OutputParameterUnsignedInt p_iteration;
#pragma endregion

			Scalar* previous_sln_vector;
      bool use_initial_guess_for_iterative_solvers;
      friend class NonlinearConvergenceMeasurement<Scalar>;
    };
  }

  namespace Exceptions
  {
    class HERMES_API NonlinearException : public Hermes::Exceptions::Exception
    {
    public:
      NonlinearException(Hermes::Solvers::NonlinearConvergenceState convergenceState);

      Hermes::Solvers::NonlinearConvergenceState get_exception_state();

    protected:
      Hermes::Solvers::NonlinearConvergenceState convergenceState;
    };
  }
}
#endif