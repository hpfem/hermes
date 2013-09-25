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
#include "solver/newton_solver.h"
#include "projections/ogprojection.h"
#include "hermes_common.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver() : NonlinearSolver<Scalar>()
    {
      init_newton();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(DiscreteProblem<Scalar>* dp) : NonlinearSolver<Scalar>(dp)
    {
      init_newton();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : NonlinearSolver<Scalar>(wf, space)
    {
      init_newton();
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : NonlinearSolver<Scalar>(wf, spaces)
    {
      init_newton();
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::isOkay() const
    {
      return NonlinearSolver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_newton()
    {
      this->min_allowed_damping_coeff = 1E-4;
      this->manual_damping = false;
      this->auto_damping_ratio = 2.0;
      this->initial_auto_damping_factor = 1.0;
      this->sufficient_improvement_factor = 0.95;
      this->necessary_successful_steps_to_increase = 3;

      this->sufficient_improvement_factor_jacobian = 1e-1;
      this->max_steps_with_reused_jacobian = 3;

      this->coeff_vec_back = NULL;

      this->set_tolerance(1e-8, ResidualNormAbsolute);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_sufficient_improvement_factor_jacobian(double ratio)
    {
      if(ratio < 0.0)
        throw Exceptions::ValueException("sufficient_improvement_factor_jacobian", sufficient_improvement_factor_jacobian, 0.0);
      this->sufficient_improvement_factor_jacobian = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_max_steps_with_reused_jacobian(unsigned int steps)
    {
      this->max_steps_with_reused_jacobian = steps;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set)
    {
      if(min_allowed_damping_coeff_to_set < 0.0)
        throw Exceptions::ValueException("min_allowed_damping_coeff_to_set", min_allowed_damping_coeff_to_set, 0.0);
      this->min_allowed_damping_coeff = min_allowed_damping_coeff_to_set;
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::~NewtonSolver()
    {
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_manual_damping_coeff(bool onOff, double coeff)
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
    void NewtonSolver<Scalar>::set_initial_auto_damping_coeff(double coeff)
    {
      if(coeff <= 0.0 || coeff > 1.0)
        throw Exceptions::ValueException("coeff", coeff, 0.0, 1.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      else
        this->initial_auto_damping_factor = coeff;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_auto_damping_ratio(double ratio)
    {
      if(ratio <= 1.0)
        throw Exceptions::ValueException("ratio", ratio, 1.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->auto_damping_ratio = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_sufficient_improvement_factor(double ratio)
    {
      if(ratio <= 0.0)
        throw Exceptions::ValueException("ratio", ratio, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->sufficient_improvement_factor = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_necessary_successful_steps_to_increase(unsigned int steps)
    {
      if(steps < 1)
        throw Exceptions::ValueException("necessary_successful_steps_to_increase", steps, 0.0);
      if(this->manual_damping)
        this->warn("Manual damping is turned on and you called set_initial_auto_damping_coeff(), turn off manual damping first by set_manual_damping_coeff(false);");
      this->necessary_successful_steps_to_increase = steps;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec)
    {
      // If we have not converged and everything else is ok, we finish.
      if(state == NotConverged)
        return false;

      // And now the finishing states (both good and bad).
      this->finalize_solving(coeff_vec);

      // Act upon the state.
      switch(state)
      {
      case Converged:
        this->info("\tNewton: done.\n");
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
        throw Exceptions::Exception("Unknown exception in NewtonSolver.");
        break;
      default:
        throw Exceptions::Exception("Unknown ConvergenceState in NewtonSolver.");
        break;
      }

      // Return that we should finish.
      return true;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::calculate_damping_factor(unsigned int& successful_steps)
    {
      Hermes::vector<double>& damping_factors = this->get_parameter_value(p_damping_factors);

      if(this->manual_damping)
      {
        damping_factors.push_back(this->manual_damping_factor);
        return true;
      }

      double residual_norm = *(this->get_parameter_value(this->p_residual_norms).end() - 1);
      double previous_residual_norm = *(this->get_parameter_value(this->p_residual_norms).end() - 2);

      if(residual_norm < previous_residual_norm * this->sufficient_improvement_factor)
      {
        if(++successful_steps >= this->necessary_successful_steps_to_increase)
        {
          double new_damping_factor = std::min(this->initial_auto_damping_factor, this->auto_damping_ratio * damping_factors.back());
          this->info("\t\tstep successful, new damping factor: %g.", new_damping_factor);
          damping_factors.push_back(new_damping_factor);
        }
        else
        {
          this->info("\t\tstep successful, keep damping factor: %g.", damping_factors.back());
          damping_factors.push_back(damping_factors.back());
        }

        return true;
      }
      else
      {
        double current_damping_factor = damping_factors.back();
        damping_factors.pop_back();
        successful_steps = 0;
        if(current_damping_factor <= this->min_allowed_damping_coeff)
        {
          this->warn("\t\tNOT improved, damping factor at minimum level: %g.", min_allowed_damping_coeff);
          this->info("\t To decrease the minimum level, use set_min_allowed_damping_coeff()");
          throw Exceptions::NonlinearException(BelowMinDampingCoeff);
        }
        else
        {
          double new_damping_factor = (1. / this->auto_damping_ratio) * current_damping_factor;
          this->warn("\t\tNOT improved, step restarted with factor: %g.", new_damping_factor);
          damping_factors.push_back(new_damping_factor);
        }

        return false;
      }
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_solving(Scalar*& coeff_vec)
    {
      this->check();
      this->tick();

      // Number of DOFs.
      this->ndof = Space<Scalar>::assign_dofs(this->get_spaces());

      // coeff_vec
      this->delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = (Scalar*)calloc(this->ndof, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      // coeff_vec_back
      this->coeff_vec_back = (Scalar*)calloc(this->ndof, sizeof(Scalar));

      // Backup vector for unsuccessful reuse of Jacobian.
      residual_back = create_vector<Scalar>();
      residual_back->alloc(this->ndof);

      // sln_vector
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }
      this->sln_vector = new Scalar[this->ndof];

      this->on_initialization();

      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      // UMFPACK reporting.
      if(this->do_UMFPACK_reporting)
        memset(this->UMFPACK_reporting_data, 0, 3 * sizeof(double));
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::deinit_solving(Scalar* coeff_vec)
    {
      if(this->delete_coeff_vec)
      {
        ::free(coeff_vec);
        this->delete_coeff_vec = false;
      }

      ::free(coeff_vec_back);
      delete residual_back;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::finalize_solving(Scalar* coeff_vec)
    {
      memcpy(this->sln_vector, coeff_vec, this->ndof * sizeof(Scalar));
      this->tick();
      this->num_iters = this->get_current_iteration_number();
      this->info("\tNewton: solution duration: %f s.\n", this->last());
      this->on_finish();
      this->deinit_solving(coeff_vec);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::force_reuse_jacobian_values(unsigned int& successful_steps_with_reused_jacobian)
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
    bool NewtonSolver<Scalar>::do_initial_step_return_finished(Scalar* coeff_vec)
    {
      // Store the initial norm.
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(coeff_vec, this->ndof));

      // Assemble the system.
      if(this->jacobian_reusable && this->reuse_jacobian_values())
      {
        this->matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
        this->dp->assemble(coeff_vec, this->residual);
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(false);
      }
      else
      {
        this->matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        this->dp->assemble(coeff_vec, this->jacobian, this->residual);
        this->jacobian_reusable = true;
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(true);
      }
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);

      this->residual->change_sign();

      // Output.
      this->process_vector_output(this->residual, 1);
      this->process_matrix_output(this->jacobian, 1);

      // Get the residual norm and act upon it.
      double residual_norm = this->calculate_residual_norm();
      this->info("\n\tNewton: initial residual norm: %g", residual_norm);
      if(residual_norm > this->max_allowed_residual_norm)
      {
        this->finalize_solving(coeff_vec);
        throw Exceptions::NonlinearException(AboveMaxAllowedResidualNorm);
        return true;
      }
      else
        this->get_parameter_value(this->p_residual_norms).push_back(this->calculate_residual_norm());

      this->solve_linear_system(coeff_vec);

      return (this->on_initial_step_end() == false);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::step_info()
    {
      // Output.
      this->info("\n\tNewton: iteration %d,", this->get_current_iteration_number());
      this->info("\t\tresidual norm: %g,", this->get_parameter_value(this->p_residual_norms).back());
      this->info("\t\tsolution norm: %g,", this->get_parameter_value(this->p_solution_norms).back());
      this->info("\t\tsolution change norm: %g.", this->get_parameter_value(this->p_solution_change_norms).back());
      this->info("\t\trelative solution change: %g.", this->get_parameter_value(this->p_solution_change_norms).back() / this->get_parameter_value(this->p_solution_norms)[this->get_parameter_value(this->p_solution_norms).size() - 2]);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::assemble_residual(Scalar* coeff_vec)
    {
      // Assemble just the residual vector & change its sign.
      this->dp->assemble(coeff_vec, this->residual);
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);
      this->residual->change_sign();

      // Output to disk.
      this->process_vector_output(this->residual, this->get_current_iteration_number());

      // Current residual norm.
      this->get_parameter_value(this->p_residual_norms).push_back(this->calculate_residual_norm());
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve_linear_system(Scalar* coeff_vec)
    {
      // store the previous coeff_vec to coeff_vec_back.
      memcpy(coeff_vec_back, coeff_vec, sizeof(Scalar)*this->ndof);

      // Solve, if the solver is iterative, give him the initial guess.
      this->matrix_solver->solve(coeff_vec);
      this->handle_UMFPACK_reports();

      // Get current damping factor.
      double current_damping_factor = this->get_parameter_value(this->p_damping_factors).back();

      // store the solution norm change.
      // obtain the solution increment.
      Scalar* sln_vector_local = this->matrix_solver->get_sln_vector();

      // 1. store the solution.
      for (int i = 0; i < this->ndof; i++)
        coeff_vec[i] += current_damping_factor * sln_vector_local[i];

      // 2. store the solution change.
      this->get_parameter_value(this->p_solution_change_norms).push_back(current_damping_factor * get_l2_norm(sln_vector_local, this->ndof));

      // 3. store the solution norm.
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(coeff_vec, this->ndof));
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
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
#pragma endregion

      if(this->do_initial_step_return_finished(coeff_vec))
      {
        this->info("\tNewton: aborted.");
        this->finalize_solving(coeff_vec);
        return;
      }

      // Main Newton loop 
      while (true)
      {
        // Handle the event of step beginning.
        if(!this->on_step_begin())
        {
          this->info("\tNewton: aborted.");
          this->finalize_solving(coeff_vec);
          return;
        }

#pragma region damping_factor_loop
        this->info("\n\tNewton: Damping factor handling:");
        // Loop searching for the damping factor.
        do
        {
          // Assemble just the residual.
          this->assemble_residual(coeff_vec);

          // Test convergence - if in this loop we found a solution.
          this->info("\t\ttest convergence...");
          this->step_info();
          if(this->handle_convergence_state_return_finished(this->get_convergence_state(), coeff_vec))
            return;
          else
            this->info("\t\thas not converged.");

          // Inspect the damping factor.
          try
          {
            // Calculate damping factor, and return whether or not was this a successful step.
            residual_norm_drop = this->calculate_damping_factor(successful_steps_damping);
          }
          catch (Exceptions::NonlinearException& e)
          {
            this->finalize_solving(coeff_vec);
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
            for (int i = 0; i < this->ndof; i++)
              coeff_vec[i] = coeff_vec_back[i] + (coeff_vec[i] - coeff_vec_back[i]) / this->auto_damping_ratio;

            // Add new solution norm.
            solution_norms.push_back(get_l2_norm(coeff_vec, this->ndof));
          }
        }
        while (!residual_norm_drop);
#pragma endregion

        // Damping factor was updated, handle the event.
        this->on_damping_factor_updated();

#pragma region jacobian_reusage_loop
        this->info("\n\tNewton: Jacobian handling:");
        // Loop until jacobian is not reusable anymore.
        // The whole loop is skipped if the jacobian is not suitable for being reused at all.
        while(this->jacobian_reusable && (this->reuse_jacobian_values() || force_reuse_jacobian_values(successful_steps_jacobian)))
        {
          this->residual_back->set_vector(this->residual);

          // Info & handle the situation as necessary.
          this->info("\t\treusing Jacobian.");
          this->on_reused_jacobian_step_begin();

          // Solve the system.
          this->matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
          this->solve_linear_system(coeff_vec);
          // Assemble next residual for both reusage and convergence test.
          this->assemble_residual(coeff_vec);
          // Test whether it was okay to reuse the jacobian.
          if(!this->jacobian_reused_okay(successful_steps_jacobian))
          {
            this->warn("\t\treused Jacobian disapproved.");
            this->get_parameter_value(this->p_residual_norms).pop_back();
            this->get_parameter_value(this->p_solution_norms).pop_back();
            this->get_parameter_value(this->p_solution_change_norms).pop_back();
            memcpy(coeff_vec, coeff_vec_back, sizeof(Scalar)*this->ndof);
            this->residual->set_vector(residual_back);
            break;
          }

          // Increase the iteration count.
          this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(false);

          // Output successful reuse info.
          this->step_info();

          // Handle the event of end of a step.
          this->on_reused_jacobian_step_end();
          if(!this->on_step_end())
          {
            this->info("\tNewton: aborted.");
            this->finalize_solving(coeff_vec);
            return;
          }

          // Test convergence - if in this iteration we found a solution.
          if(this->handle_convergence_state_return_finished(this->get_convergence_state(), coeff_vec))
            return;
        }
#pragma endregion

        // Reassemble the jacobian once not reusable anymore.
        this->info("\t\tre-calculating Jacobian.");
        this->dp->assemble(coeff_vec, this->jacobian);
        if(this->report_cache_hits_and_misses)
          this->add_cache_hits_and_misses(this->dp);

        // Set factorization schemes.
        if(this->jacobian_reusable)
          this->matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_REORDERING);
        else
          this->matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

        // Solve the system, state that the jacobian is reusable should it be desirable.
        this->solve_linear_system(coeff_vec);
        this->jacobian_reusable = true;

        // Handle the event of end of a step.
        if(!this->on_step_end())
        {
          this->info("\tNewton: aborted.");
          this->finalize_solving(coeff_vec);
          return;
        }

        // Increase the iteration count.
        this->get_parameter_value(this->p_iterations_with_recalculated_jacobian).push_back(true);
        // Output info.
        this->step_info();
      }
    }

    template<typename Scalar>
    int NewtonSolver<Scalar>::get_current_iteration_number()
    {
      return this->get_parameter_value(p_iterations_with_recalculated_jacobian).size() + 1;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::jacobian_reused_okay(unsigned int& successful_steps_with_reused_jacobian)
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
    void NewtonSolver<Scalar>::on_damping_factor_updated()
    {
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::on_reused_jacobian_step_begin()
    {
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::on_reused_jacobian_step_end()
    {
    }

    template class HERMES_API NewtonSolver<double>;
    template class HERMES_API NewtonSolver<std::complex<double> >;
  }
}
