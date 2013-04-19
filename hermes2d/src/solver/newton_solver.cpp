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
      this->current_convergence_measurement = ResidualNormAbsolute;
      this->newton_tolerance = 1e-8;
      this->max_allowed_iterations = 15;
      this->residual_as_function = false;
      this->max_allowed_residual_norm = 1E9;
      this->min_allowed_damping_coeff = 1E-4;
      this->manual_damping = false;
      this->auto_damping_ratio = 2.0;
      this->initial_auto_damping_coefficient = 1.0;
      this->sufficient_improvement_factor = 0.95;
      this->necessary_successful_steps_to_increase = 3;

      this->sufficient_improvement_factor_jacobian = 1e-1;
      this->max_steps_with_reused_jacobian = 3;

      this->coeff_vec_back = NULL;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_sufficient_improvement_factor_jacobian(double ratio)
    {
      this->sufficient_improvement_factor_jacobian = ratio;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_max_steps_with_reused_jacobian(unsigned int steps)
    {
      this->max_steps_with_reused_jacobian = steps;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_tolerance(double tolerance_)
    {
      this->newton_tolerance = tolerance_;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set)
    {
      if(max_allowed_residual_norm_to_set <= 0.0)
        throw Exceptions::ValueException("max_allowed_residual_norm_to_set", max_allowed_residual_norm_to_set, 0.0);
      this->max_allowed_residual_norm = max_allowed_residual_norm_to_set;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_min_allowed_damping_coeff(double min_allowed_damping_coeff_to_set)
    {
      if(min_allowed_damping_coeff_to_set <= 0.0)
        throw Exceptions::ValueException("min_allowed_damping_coeff_to_set", min_allowed_damping_coeff_to_set, 0.0);
      this->min_allowed_damping_coeff = min_allowed_damping_coeff_to_set;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_residual_as_function()
    {
      this->residual_as_function = true;
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
        this->manual_damping_coefficient = coeff;
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
        this->initial_auto_damping_coefficient = coeff;
    }
      
    template<typename Scalar>
    void NewtonSolver<Scalar>::set_auto_damping_ratio(double ratio)
    {
      if(ratio <= 0.0)
        throw Exceptions::ValueException("ratio", ratio, 0.0);
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
    void NewtonSolver<Scalar>::set_convergence_measurement(NewtonSolverConvergenceMeasurement measurement)
    {
      this->current_convergence_measurement = measurement;
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::test_convergence(Scalar* coeff_vec)
    {
      // Get the state.
      NewtonSolver<Scalar>::ConvergenceState state;

      unsigned int iteration = this->get_parameter_value(this->p_iteration);
      double residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 1];
      double current_damping_coefficient = this->get_parameter_value(p_current_damping_coefficient);

      if(residual_norm > this->max_allowed_residual_norm)
        state = AboveMaxAllowedResidualNorm;

      if(iteration >= this->max_allowed_iterations)
        state = AboveMaxIterations;

      if(newtonConverged<Scalar>(this))
        state = Converged;
      else
        state = NotConverged;

      // Act upon the state.
      switch(state)
      {
      case Converged:
        // We want to return the solution in a different structure.
        memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
        this->finalize_solving(coeff_vec);
        return true;
        break;

      case AboveMaxIterations:
        memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
        throw NewtonException(AboveMaxIterations);
        this->finalize_solving(coeff_vec);
        return false;
        break;

      case BelowMinDampingCoeff:
        memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
        throw NewtonException(BelowMinDampingCoeff);
        this->finalize_solving(coeff_vec);
        return false;
        break;

      case AboveMaxAllowedResidualNorm:
        memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
        throw NewtonException(AboveMaxAllowedResidualNorm);
        this->finalize_solving(coeff_vec);
        return false;
        break;

      case Error:
        throw Exceptions::Exception("Unknown exception in NewtonSolver.");
        this->finalize_solving(coeff_vec);
        return false;
        break;

      default:
        // The only state here is NotConverged which yields staying in the loop.
        return false;
        break;
      }

      return false;
    }

    template<typename Scalar>
    double NewtonSolver<Scalar>::calculate_residual_norm()
    {
      // Measure the residual norm.
      if(residual_as_function)
      {
        // Prepare solutions for measuring residual norm.
        Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions;
        Hermes::vector<MeshFunction<Scalar>*> solutionsPtrs;
        Hermes::vector<bool> dir_lift_false;
        for (unsigned int i = 0; i < this->dp->get_spaces().size(); i++) 
        {
          MeshFunctionSharedPtr<Scalar> sharedPtr(new Solution<Scalar>());
          solutions.push_back(sharedPtr);
          solutionsPtrs.push_back(sharedPtr.get());
          dir_lift_false.push_back(false);
        }

        Solution<Scalar>::vector_to_solutions(this->residual, this->dp->get_spaces(), solutions, dir_lift_false);

        // Calculate the norm.
        return Global<Scalar>::calc_norms(solutionsPtrs);
      }
      else
      {
        // Calculate the l2-norm of residual vector, this is the traditional way.
        return Global<Scalar>::get_l2_norm(this->residual);
      }
    }

    template<typename Scalar>
    double NewtonSolver<Scalar>::calculate_damping_coefficient(bool& residual_norm_drop, unsigned int& successful_steps)
    {
      if(this->manual_damping)
        return this->manual_damping_coefficient;

      int iteration = this->get_parameter_value(p_iteration);
      double residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 1];
      double previous_residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 2];
      double current_damping_coefficient = this->get_parameter_value(p_current_damping_coefficient);

      if(residual_norm < previous_residual_norm * this->sufficient_improvement_factor)
      {
        residual_norm_drop = true;
        if(++successful_steps >= this->necessary_successful_steps_to_increase)
        {
          current_damping_coefficient = std::min(this->initial_auto_damping_coefficient, this->auto_damping_ratio * current_damping_coefficient);
          this->info("\tNewton: step successful, damping coefficient: %g.", current_damping_coefficient);
          return current_damping_coefficient;
        }
        if(residual_norm < previous_residual_norm)
        {
          this->info("\tNewton: step successful, damping coefficient: %g.", current_damping_coefficient);
          return current_damping_coefficient;
        }
      }
      else
      {
        successful_steps = 0;
        residual_norm_drop = false;
        if(current_damping_coefficient <= this->min_allowed_damping_coeff)
        {
          this->warn("\tNewton: results NOT improved, current damping coefficient is at the minimum possible level: %g.", min_allowed_damping_coeff);
          this->info("\t  If you want to decrease the minimum level, use the method set_min_allowed_damping_coeff()");
          throw NewtonException(BelowMinDampingCoeff);
          return this->min_allowed_damping_coeff;
        }
        else
        {
          current_damping_coefficient = (1. / this->auto_damping_ratio) * current_damping_coefficient;
          this->warn("\tNewton: results NOT improved, step restarted with damping coefficient: %g.", current_damping_coefficient);
          return current_damping_coefficient;
        }
      }
      return current_damping_coefficient;
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
        coeff_vec = (Scalar*)calloc(ndof, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      // coeff_vec_back
      this->coeff_vec_back = (Scalar*)calloc(ndof, sizeof(Scalar));

      // sln_vector
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }
      this->sln_vector = new Scalar[ndof];

      this->on_initialization();
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
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::finalize_solving(Scalar* coeff_vec)
    {
      this->tick();
      this->info("\tNewton: solution duration: %f s.\n", this->last());
      this->on_finish();
      this->deinit_solving(coeff_vec);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::force_reuse_jacobian_values(unsigned int& successful_steps_with_reused_jacobian)
    {
      int iteration = this->get_parameter_value(p_iteration);
      if(iteration == 1)
        return false;

      double residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 1];
      double previous_residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 2];
      
      if(successful_steps_with_reused_jacobian >= this->max_steps_with_reused_jacobian)
      {
        successful_steps_with_reused_jacobian = 0;
        return false;
      }
      if((residual_norm / previous_residual_norm) > this->sufficient_improvement_factor_jacobian)
      {
        successful_steps_with_reused_jacobian = 0;
        return false;
      }
      successful_steps_with_reused_jacobian++;
      return true;
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonException::NewtonException(typename NewtonSolver<Scalar>::ConvergenceState convergenceState) : convergenceState(convergenceState)
    {
    }

    template<typename Scalar>
    typename NewtonSolver<Scalar>::ConvergenceState NewtonSolver<Scalar>::NewtonException::get_exception_state()
    {
      return this->convergenceState;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::do_initial_step(Scalar* coeff_vec)
    {
      // Store the initial norm.
      this->get_parameter_value(p_solution_norms).push_back(Global<Scalar>::get_l2_norm(coeff_vec, this->ndof));

      // Assemble the system.
      if(this->jacobian_reusable && this->reuse_jacobian_values())
      {
        this->dp->assemble(coeff_vec, this->residual);
      }
      else
      {
        this->dp->assemble(coeff_vec, this->jacobian, this->residual);
        this->jacobian_reusable = true;
      }
      this->residual->change_sign();

      // Output.
      this->process_vector_output(this->residual, 1);
      this->process_matrix_output(this->jacobian, 1);

      // Get the residual norm and act upon it.
      double residual_norm = this->calculate_residual_norm();
      this->info("\tNewton: initial residual norm: %g", residual_norm);
      if(residual_norm > this->max_allowed_residual_norm)
      {
        this->finalize_solving(coeff_vec);
        throw NewtonException(AboveMaxAllowedResidualNorm);
      }
      else
        this->get_parameter_value(p_residual_norms).push_back(this->calculate_residual_norm());

      this->solve_linear_system(coeff_vec);

      if(!this->on_initial_step_end())
      {
        memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
        this->info("Aborted");
        this->finalize_solving(coeff_vec);
        return;
      }

      this->get_parameter_value(p_iteration)++;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::assemble_residual(Scalar* coeff_vec)
    {
      // Assemble just the residual vector & change its sign.
      this->dp->assemble(coeff_vec, this->residual);
      this->residual->change_sign();

      double residual_norm = this->calculate_residual_norm();
      // Output.
      this->info("\tNewton: iteration %d, residual norm: %g", this->get_parameter_value(p_iteration), residual_norm);

      if(residual_norm > this->max_allowed_residual_norm)
      {
        this->finalize_solving(coeff_vec);
        throw NewtonException(AboveMaxAllowedResidualNorm);
      }
      else
      {
        // Output to disk.
        this->process_vector_output(this->residual, this->get_parameter_value(p_iteration));
        // Current residual norm && current_damping_coefficient.
        this->get_parameter_value(p_residual_norms).push_back(this->calculate_residual_norm());
      }
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve_linear_system(Scalar* coeff_vec)
    {
      if(this->matrix_solver->solve())
      {
        // store the previous coeff_vec to coeff_vec_back.
        memcpy(coeff_vec_back, coeff_vec, sizeof(Scalar)*ndof);

        // obtain the solution increment.
        Scalar* sln_vector_local = this->matrix_solver->get_sln_vector();

        double damping_coeff = this->get_parameter_value(p_current_damping_coefficient);

        // store the solution norm change.
        this->get_parameter_value(p_solution_change_norm) = damping_coeff * Global<Scalar>::get_l2_norm(sln_vector_local, ndof);

        // add the increment to the solution.
        for (int i = 0; i < ndof; i++)
          coeff_vec[i] += damping_coeff * sln_vector_local[i];

        // store the solution norm.
        this->get_parameter_value(p_solution_norms).push_back(Global<Scalar>::get_l2_norm(coeff_vec, this->ndof));
      }
      else
      {
        this->deinit_solving(coeff_vec);
        throw Exceptions::LinearMatrixSolverException();
      }
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      // Initialization.
      this->init_solving(coeff_vec);

#pragma region parameter_setup
      // The Newton's loop.
      unsigned int it = 1;
      unsigned int successful_steps_damping = 0;
      unsigned int successful_steps_jacobian = 0;
      Hermes::vector<double> residual_norms;
      Hermes::vector<double> solution_norms;
      double solution_change_norm;
      double current_damping_coefficient = this->manual_damping ? manual_damping_coefficient : initial_auto_damping_coefficient;
      bool residual_norm_drop = true;

      this->set_parameter_value(this->p_residual_norms, &residual_norms);
      this->set_parameter_value(this->p_solution_norms, &solution_norms);
      this->set_parameter_value(this->p_solution_change_norm, &solution_change_norm);
      this->set_parameter_value(this->p_successful_steps_damping, &successful_steps_damping);
      this->set_parameter_value(this->p_successful_steps_jacobian, &successful_steps_jacobian);
      this->set_parameter_value(this->p_iteration, &it);
      this->set_parameter_value(this->p_residual_norm_drop, &residual_norm_drop);
      this->set_parameter_value(this->p_current_damping_coefficient, &current_damping_coefficient);
#pragma endregion

      this->do_initial_step(coeff_vec);

      // Main Newton loop 
      while (true)
      {
        // User method.
        this->on_step_begin();
        
        // Loop searching for the damping coefficient.
        do
        {
          // Assemble just the residual.
          this->assemble_residual(coeff_vec);

          // Test convergence - if in this loop we found a solution.
          if(this->test_convergence(coeff_vec))
            return;

          // Inspect the damping coefficient.
          try
          {
            current_damping_coefficient = this->calculate_damping_coefficient(residual_norm_drop, successful_steps_damping);
          }
          catch (NewtonException& e)
          {
            this->finalize_solving(coeff_vec);
            throw;
          }

          if(!residual_norm_drop)
          {
            // Delete the previous residual norm.
            residual_norms.pop_back();

            // Try with the different damping coefficient.
            for (int i = 0; i < ndof; i++)
              coeff_vec[i] = coeff_vec_back[i] + current_damping_coefficient * (coeff_vec[i] - coeff_vec_back[i]);
          }
        }
        while (!residual_norm_drop);

        // Loop until jacobian is reusable.
        while(this->jacobian_reusable && (this->reuse_jacobian_values() || force_reuse_jacobian_values(successful_steps_jacobian)))
        {
          // Solve the system.
          this->matrix_solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
          this->solve_linear_system(coeff_vec);

          // User method call.
          if(!this->on_step_end())
          {
            memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
            this->info("Aborted");
            this->finalize_solving(coeff_vec);
            return;
          }

          // Increase the iteration count.
          it++;

          // Assemble next residual for convergence test.
          this->assemble_residual(coeff_vec);

          // Test convergence - if in this loop we found a solution.
          if(this->test_convergence(coeff_vec))
            return;
        }

        // Reassemble the jacobian once not reusable anymore.
        this->dp->assemble(coeff_vec, this->jacobian);

        // Set factorization schemes.
        if(this->jacobian_reusable)
          this->matrix_solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING_AND_SCALING);
        else
        this->matrix_solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);

        this->solve_linear_system(coeff_vec);
        this->jacobian_reusable = true;

        if(!this->on_step_end())
        {
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          this->info("Aborted");
          this->finalize_solving(coeff_vec);
          return;
        }

        // Increase the iteration count.
        it++;
      }
    }

    template class HERMES_API NewtonSolver<double>;
    template class HERMES_API NewtonSolver<std::complex<double> >;
  }
}