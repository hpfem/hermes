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
/*! \file nonlinear_solver.cpp
\brief General nonlinear solver functionality.
*/
#include "solvers/nonlinear_matrix_solver.h"
#include "common.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    NonlinearMatrixSolver<Scalar>::NonlinearMatrixSolver() : MatrixSolver<Scalar>()
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::init_nonlinear()
    {
      this->handleMultipleTolerancesAnd = false;
      this->max_allowed_iterations = 20;
      this->max_allowed_residual_norm = 1E9;
      this->num_iters = 0;
      this->previous_sln_vector = nullptr;
      this->use_initial_guess_for_iterative_solvers = false;
      this->clear_tolerances();
    }

    template<typename Scalar>
    NonlinearMatrixSolver<Scalar>::~NonlinearMatrixSolver()
    {
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::isOkay() const
    {
      bool toleranceSet = false;
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        if(this->tolerance_set[i])
          toleranceSet = true;
      if(!toleranceSet)
      {
        throw Exceptions::Exception("No tolerance set in NonlinearMatrixSolver.");
        return false;
      }
      return true;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_max_allowed_iterations(int max_allowed_iterations_)
    {
      if(max_allowed_iterations_ < 1)
        throw Exceptions::ValueException("max_allowed_iterations", max_allowed_iterations_, 1);
      this->max_allowed_iterations = max_allowed_iterations_;
    }

    template<typename Scalar>
    double NonlinearMatrixSolver<Scalar>::calculate_residual_norm()
    {
      return get_l2_norm(this->get_residual());
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::init_solving(Scalar* coeff_vec)
    {
      this->check();
      this->tick();

      // Number of DOFs.
      assert(this->problem_size > 0);

      if(this->sln_vector != nullptr)
      {
        delete [] this->sln_vector;
        this->sln_vector = nullptr;
      }

      this->sln_vector = new Scalar[this->problem_size];

      if(coeff_vec == nullptr)
        memset(this->sln_vector, 0, this->problem_size*sizeof(Scalar));
      else
        memcpy(this->sln_vector, coeff_vec, this->problem_size*sizeof(Scalar));

      // previous_sln_vector
      this->previous_sln_vector = (Scalar*)calloc(this->problem_size, sizeof(Scalar));

      // Backup vector for unsuccessful reuse of Jacobian.
      residual_back = create_vector<Scalar>();
      residual_back->alloc(this->problem_size);

      this->previous_jacobian = nullptr;
      this->previous_residual = nullptr;

      this->on_initialization();
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_tolerance(double tolerance_, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
    {
      this->handleMultipleTolerancesAnd = handleMultipleTolerancesAnd;

      if(tolerance_ < 0.0)
        throw Exceptions::ValueException("tolerance", tolerance_, 0.0);

      switch(toleranceType)
      {
      case ResidualNormRelativeToInitial:
        {
          this->tolerance[0] = tolerance_;
          this->tolerance_set[0] = true;
        }
        break;
      case ResidualNormRelativeToPrevious:
        {
          this->tolerance[1] = tolerance_;
          this->tolerance_set[1] = true;
        }
        break;
      case ResidualNormRatioToInitial:
        {
          this->tolerance[2] = tolerance_;
          this->tolerance_set[2] = true;
        }
        break;
      case ResidualNormRatioToPrevious:
        {
          this->tolerance[3] = tolerance_;
          this->tolerance_set[3] = true;
        }
        break;
      case ResidualNormAbsolute:
        {
          this->tolerance[4] = tolerance_;
          this->tolerance_set[4] = true;
        }
        break;
      case SolutionChangeAbsolute:
        {
          this->tolerance[5] = tolerance_;
          this->tolerance_set[5] = true;
        }
        break;
      case SolutionChangeRelative:
        {
          this->tolerance[6] = tolerance_;
          this->tolerance_set[6] = true;
        }
        break;
      default:
        throw Exceptions::Exception("Unknown NonlinearConvergenceMeasurementType in NonlinearMatrixSolver::set_tolerance.");
      }
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set)
    {
      if(max_allowed_residual_norm_to_set < 0.0)
        throw Exceptions::ValueException("max_allowed_residual_norm_to_set", max_allowed_residual_norm_to_set, 0.0);
      this->max_allowed_residual_norm = max_allowed_residual_norm_to_set;
    }

    template<typename Scalar>
    int NonlinearMatrixSolver<Scalar>::get_num_iters() const
    {
      return this->num_iters;
    }

    template<typename Scalar>
    NonlinearConvergenceState NonlinearMatrixSolver<Scalar>::get_convergence_state()
    {
      if(this->get_current_iteration_number() >= this->max_allowed_iterations)
        return AboveMaxIterations;

      if(NonlinearConvergenceMeasurement<Scalar>::converged(this))
        return Converged;
      else
        return NotConverged;

      return Error;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::clear_tolerances()
    {
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        this->tolerance[i] = std::numeric_limits<double>::max();
      memset(this->tolerance_set, 0, sizeof(bool)*NonlinearConvergenceMeasurementTypeCount);
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_sufficient_improvement_factor_jacobian(double ratio)
    {
      if(ratio < 0.0)
        throw Exceptions::ValueException("sufficient_improvement_factor_jacobian", sufficient_improvement_factor_jacobian, 0.0);
      this->sufficient_improvement_factor_jacobian = ratio;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_max_steps_with_reused_jacobian(unsigned int steps)
    {
      this->max_steps_with_reused_jacobian = steps;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set)
    {
      if(min_allowed_damping_coeff_to_set < 0.0)
        throw Exceptions::ValueException("min_allowed_damping_coeff_to_set", min_allowed_damping_coeff_to_set, 0.0);
      this->min_allowed_damping_coeff = min_allowed_damping_coeff_to_set;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_manual_damping_coeff(bool onOff, double coeff)
    {
      if(coeff <= 0.0 || coeff > 1.0)
        throw Exceptions::ValueException("coeff", coeff, 0.0, 1.0);
      if(onOff)
      {
        this->manual_damping = true;
        this->manual_damping_factor = coeff;
      }
      else
        this->manual_damping = false;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_initial_auto_damping_coeff(double coeff)
    {
      if(coeff <= 0.0 || coeff > 1.0)
        throw Exceptions::ValueException("coeff", coeff, 0.0, 1.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      else
        this->initial_auto_damping_factor = coeff;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_auto_damping_ratio(double ratio)
    {
      if(ratio <= 1.0)
        throw Exceptions::ValueException("ratio", ratio, 1.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->auto_damping_ratio = ratio;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_sufficient_improvement_factor(double ratio)
    {
      if(ratio <= 0.0)
        throw Exceptions::ValueException("ratio", ratio, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->sufficient_improvement_factor = ratio;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_necessary_successful_steps_to_increase(unsigned int steps)
    {
      if(steps < 1)
        throw Exceptions::ValueException("necessary_successful_steps_to_increase", steps, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->necessary_successful_steps_to_increase = steps;
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::handle_convergence_state_return_finished(NonlinearConvergenceState state)
    {
      // If we have not converged and everything else is ok, we finish.
      if(state == NotConverged)
        return false;

      // And now the finishing states (both good and bad).
      this->finalize_solving();

      // Act upon the state.
      switch(state)
      {
      case Converged:
        this->info("\tNonlinearSolver: done.\n");
        break;
      case AboveMaxIterations:
        throw Exceptions::NonlinearException(AboveMaxIterations);
        break;
      case BelowMinDampingCoeff:
        throw Exceptions::NonlinearException(BelowMinDampingCoeff);
        break;
      case AboveMaxAllowedResidualNorm:
        throw Exceptions::NonlinearException(AboveMaxAllowedResidualNorm);
        break;
      case Error:
        throw Exceptions::Exception("Unknown exception in NonlinearMatrixSolver.");
        break;
      default:
        throw Exceptions::Exception("Unknown ConvergenceState in NonlinearMatrixSolver.");
        break;
      }

      // Return that we should finish.
      return true;
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::calculate_damping_factor(unsigned int& successful_steps)
    {
      Hermes::vector<double>& damping_factors_vector = this->get_parameter_value(p_damping_factors);

      if(this->manual_damping)
      {
        damping_factors_vector.push_back(this->manual_damping_factor);
        return true;
      }
      
      if(this->damping_factor_condition())
      {
        if(++successful_steps >= this->necessary_successful_steps_to_increase)
        {
          double new_damping_factor = std::min(this->initial_auto_damping_factor, this->auto_damping_ratio * damping_factors_vector.back());
          this->info("\t\tstep successful, new damping factor: %g.", new_damping_factor);
          damping_factors_vector.push_back(new_damping_factor);
        }
        else
        {
          this->info("\t\tstep successful, keep damping factor: %g.", damping_factors_vector.back());
          damping_factors_vector.push_back(damping_factors_vector.back());
        }

        return true;
      }
      else
      {
        double current_damping_factor = damping_factors_vector.back();
        damping_factors_vector.pop_back();
        successful_steps = 0;
        if(current_damping_factor <= this->min_allowed_damping_coeff)
        {
          this->warn("\t\tNOT successful, damping factor at minimum level: %g.", min_allowed_damping_coeff);
          this->info("\t\tto decrease the minimum level, use set_min_allowed_damping_coeff()");
          throw Exceptions::NonlinearException(BelowMinDampingCoeff);
        }
        else
        {
          double new_damping_factor = (1. / this->auto_damping_ratio) * current_damping_factor;
          this->warn("\t\tNOT successful, step restarted with factor: %g.", new_damping_factor);
          damping_factors_vector.push_back(new_damping_factor);
        }

        return false;
      }
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::deinit_solving()
    {
      ::free(this->previous_sln_vector);
      delete residual_back;
      this->problem_size = -1;
      if (this->previous_jacobian)
        delete this->previous_jacobian;
      if (this->previous_residual)
        delete this->previous_residual;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::finalize_solving()
    {
      this->tick();
      this->num_iters = this->get_current_iteration_number();
      this->info("\tNonlinearSolver: solution duration: %f s.\n", this->last());
      this->on_finish();
      this->deinit_solving();
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::force_reuse_jacobian_values(unsigned int& successful_steps_with_reused_jacobian)
    {
      if(successful_steps_with_reused_jacobian >= this->max_steps_with_reused_jacobian)
      {
        successful_steps_with_reused_jacobian = 0;
        return false;
      }
      successful_steps_with_reused_jacobian++;
      return true;
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::do_initial_step_return_finished()
    {
      // Store the initial norm.
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(this->sln_vector, this->problem_size));

      // Assemble the system.
      if(this->jacobian_reusable && this->constant_jacobian)
      {
        this->linear_matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
        this->assemble_residual(false);
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(false);
      }
      else
      {
        this->linear_matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        this->assemble(false, false);
        this->jacobian_reusable = true;
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(true);
      }

      // Get the residual norm and act upon it.
      double residual_norm = this->calculate_residual_norm();
      this->info("\n\tNonlinearSolver: initial residual norm: %g", residual_norm);
      this->get_parameter_value(this->p_residual_norms).push_back(residual_norm);

      this->solve_linear_system();

      if(this->handle_convergence_state_return_finished(this->get_convergence_state()))
        return true;

      if(this->on_initial_step_end() == false)
      {
        this->info("\tNonlinearSolver: aborted.");
        this->finalize_solving();
        return true;
      }

      return false;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::solve_linear_system()
    {
      // store the previous solution to previous_sln_vector.
      memcpy(this->previous_sln_vector, this->sln_vector, sizeof(Scalar)*this->problem_size);

      // Solve, if the solver is iterative, give him the initial guess.
      this->linear_matrix_solver->solve(this->use_initial_guess_for_iterative_solvers ? this->sln_vector : nullptr);

      // 1. store the solution.
      double solution_change_norm = this->update_solution_return_change_norm(this->linear_matrix_solver->get_sln_vector());

      // 2. store the solution change.
      this->get_parameter_value(this->p_solution_change_norms).push_back(solution_change_norm);

      // 3. store the solution norm.
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(this->sln_vector, this->problem_size));
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::step_info()
    {
      // Output.
      this->info("\n\tNonlinearSolver: iteration %d,", this->get_current_iteration_number());
      this->info("\t\tresidual norm: %g,", this->get_parameter_value(this->p_residual_norms).back());

      double solution_norm = this->get_parameter_value(this->p_solution_norms).back();
      double previous_solution_norm = this->get_parameter_value(this->p_solution_norms)[this->get_parameter_value(this->p_solution_norms).size() - 2];
      double solution_change_norm = this->get_parameter_value(this->p_solution_change_norms).back();
      this->info("\t\tsolution norm: %g,", solution_norm);
      this->info("\t\tsolution change norm: %g.", solution_change_norm);
      this->info("\t\trelative solution change: %g.", solution_change_norm / previous_solution_norm);
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      // Initialization.
      this->init_solving(coeff_vec);

#pragma region parameter_setup
      // Initialize parameters.
      unsigned int successful_steps_damping = 0;
      unsigned int successful_steps_jacobian = 0;
      Hermes::vector<bool> iterations_with_recalculated_jacobian;
      Hermes::vector<double> residual_norms;
      Hermes::vector<double> solution_norms;
      Hermes::vector<double> solution_change_norms;
      Hermes::vector<double> damping_factors;
      bool residual_norm_drop;
      unsigned int iteration = 1;

      // Initial damping factor.
      damping_factors.push_back(this->manual_damping ? manual_damping_factor : initial_auto_damping_factor);

      // Link parameters.
      this->set_parameter_value(this->p_residual_norms, &residual_norms);
      this->set_parameter_value(this->p_solution_norms, &solution_norms);
      this->set_parameter_value(this->p_solution_change_norms, &solution_change_norms);
      this->set_parameter_value(this->p_successful_steps_damping, &successful_steps_damping);
      this->set_parameter_value(this->p_successful_steps_jacobian, &successful_steps_jacobian);
      this->set_parameter_value(this->p_residual_norm_drop, &residual_norm_drop);
      this->set_parameter_value(this->p_damping_factors, &damping_factors);
      this->set_parameter_value(this->p_iterations_with_recalculated_jacobian, &iterations_with_recalculated_jacobian);
      this->set_parameter_value(this->p_iteration, &iteration);
#pragma endregion

      // Does the initial step & checks the convergence & deallocates.
      if(this->do_initial_step_return_finished())
        return;

      // Main Nonlinear loop 
      while (true)
      {
        // Handle the event of step beginning.
        if(!this->on_step_begin())
        {
          this->info("\tNonlinearSolver: aborted.");
          this->finalize_solving();
          return;
        }

#pragma region damping_factor_loop
        this->info("\n\tNonlinearSolver: Damping factor handling:");
        // Loop searching for the damping factor.
        do
        {
          // Assemble just the residual.
          this->assemble_residual(false);
          // Current residual norm.
          this->get_parameter_value(this->p_residual_norms).push_back(this->calculate_residual_norm());

          // Test convergence - if in this loop we found a solution.
          this->info("\t\t\tconvergence test");
          if(this->handle_convergence_state_return_finished(this->get_convergence_state()))
            return;

          // Inspect the damping factor.
          this->info("\t\tprobing the damping factor...");
          try
          {
            // Calculate damping factor, and return whether or not was this a successful step.
            residual_norm_drop = this->calculate_damping_factor(successful_steps_damping);
          }
          catch (Exceptions::NonlinearException& e)
          {
            this->finalize_solving();
            throw;
            return;
          }

          if(!residual_norm_drop)
          {
            // Delete the previous residual and solution norm.
            residual_norms.pop_back();
            solution_norms.pop_back();

            // Adjust the previous solution change norm.
            solution_change_norms.back() /= this->auto_damping_ratio;

            // Try with the different damping factor.
            // Important thing here is the factor used that must be calculated from the current one and the previous one.
            // This results in the following relation (since the damping factor is only updated one way).
            for (int i = 0; i < this->problem_size; i++)
              this->sln_vector[i] = this->previous_sln_vector[i] + (this->sln_vector[i] - this->previous_sln_vector[i]) / this->auto_damping_ratio;

            // Add new solution norm.
            solution_norms.push_back(get_l2_norm(this->sln_vector, this->problem_size));
          }
        }
        while (!residual_norm_drop);
#pragma endregion

        // Damping factor was updated, handle the event.
        this->on_damping_factor_updated();

#pragma region jacobian_reusage_loop
        this->info("\n\tNonlinearSolver: Jacobian handling:");
        // Loop until jacobian is not reusable anymore.
        // The whole loop is skipped if the jacobian is not suitable for being reused at all.
        while(this->jacobian_reusable && (this->constant_jacobian || force_reuse_jacobian_values(successful_steps_jacobian)))
        {
          this->residual_back->set_vector(this->get_residual());

          // Info & handle the situation as necessary.
          this->info("\t\treusing Jacobian.");
          this->on_reused_jacobian_step_begin();

          // Solve the system.
          this->linear_matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
          this->solve_linear_system();

          // Assemble next residual for both reusage and convergence test.
          this->assemble_residual(false);
          // Current residual norm.
          this->get_parameter_value(this->p_residual_norms).push_back(this->calculate_residual_norm());

          // Test whether it was okay to reuse the jacobian.
          if(!this->jacobian_reused_okay(successful_steps_jacobian))
          {
            this->warn("\t\treused Jacobian disapproved.");
            this->get_parameter_value(this->p_residual_norms).pop_back();
            this->get_parameter_value(this->p_solution_norms).pop_back();
            this->get_parameter_value(this->p_solution_change_norms).pop_back();
            memcpy(this->sln_vector, this->previous_sln_vector, sizeof(Scalar)*this->problem_size);
            this->get_residual()->set_vector(residual_back);
            break;
          }

          // Increase the iteration count.
          this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(false);
          this->get_parameter_value(this->p_iteration)++;

          // Output successful reuse info.
          this->step_info();

          // Handle the event of end of a step.
          this->on_reused_jacobian_step_end();
          if(!this->on_step_end())
          {
            this->info("\tNonlinearSolver: aborted.");
            this->finalize_solving();
            return;
          }

          // Test convergence - if in this iteration we found a solution.
          if(this->handle_convergence_state_return_finished(this->get_convergence_state()))
            return;
        }
#pragma endregion

        // Reassemble the jacobian once not reusable anymore.
        this->info("\t\tre-calculating Jacobian.");
        this->assemble_jacobian(true);

        // Set factorization schemes.
        if(this->jacobian_reusable)
          this->linear_matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_REORDERING);
        else
          this->linear_matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

        // Solve the system, state that the jacobian is reusable should it be desirable.
        this->solve_linear_system();
        this->jacobian_reusable = true;

        // Increase the iteration count.
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(true);
        this->get_parameter_value(this->p_iteration)++;

        // Output info.
        this->step_info();

        // Handle the event of end of a step.
        if(!this->on_step_end())
        {
          this->info("\tNonlinearSolver: aborted.");
          this->finalize_solving();
          return;
        }

        // Test convergence - if in this iteration we found a solution.
        if(this->handle_convergence_state_return_finished(this->get_convergence_state()))
          return;
      }
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::damping_factor_condition()
    {
      double residual_norm = *(this->get_parameter_value(this->residual_norms()).end() - 1);
      double previous_residual_norm = *(this->get_parameter_value(this->residual_norms()).end() - 2);

      return (residual_norm < previous_residual_norm * this->sufficient_improvement_factor);
    }

    template<typename Scalar>
    int NonlinearMatrixSolver<Scalar>::get_current_iteration_number()
    {
      return this->get_parameter_value(this->p_iteration);
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::jacobian_reused_okay(unsigned int& successful_steps_with_reused_jacobian)
    {
      double residual_norm = *(this->get_parameter_value(this->p_residual_norms).end() - 1);
      double previous_residual_norm = *(this->get_parameter_value(this->p_residual_norms).end() - 2);

      if((residual_norm / previous_residual_norm) > this->sufficient_improvement_factor_jacobian)
      {
        successful_steps_with_reused_jacobian = 0;
        return false;
      }
      else
        return true;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::on_damping_factor_updated()
    {
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::on_reused_jacobian_step_begin()
    {
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::on_reused_jacobian_step_end()
    {
    }

    template class HERMES_API NonlinearMatrixSolver<double>;
    template class HERMES_API NonlinearMatrixSolver<std::complex<double> >;
  }

  namespace Exceptions
  {
    NonlinearException::NonlinearException(Solvers::NonlinearConvergenceState convergenceState) : Exception("NonlinearException"), convergenceState(convergenceState)
    {
    }

    Solvers::NonlinearConvergenceState NonlinearException::get_exception_state()
    {
      return this->convergenceState;
    }
  }
}