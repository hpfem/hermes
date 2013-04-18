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
    typename NewtonSolver<Scalar>::ConvergenceState NewtonSolver<Scalar>::get_convergence_state()
    {
      // If maximum allowed residual norm is exceeded, fail.
      unsigned int iteration = this->get_parameter_value(this->p_iteration);
      double residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 1];
      
      if(residual_norm > this->max_allowed_residual_norm)
        return AboveMaxAllowedResidualNorm;

      if(iteration >= this->max_allowed_iterations)
        return AboveMaxIterations;

      if(newtonConverged<Scalar>(this))
        return Converged;
      else
        return NotConverged;
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
        if(current_damping_coefficient < this->min_allowed_damping_coeff)
        {
          this->warn("\tNewton: results NOT improved, current damping coefficient is at the minimum possible level: %g.", min_allowed_damping_coeff);
          this->info("\t  If you want to decrease the minimum level, use the method set_min_allowed_damping_coeff()");
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
    void NewtonSolver<Scalar>::init_solving(int ndof, Scalar*& coeff_vec, Scalar*& coeff_vec_back)
    {
      this->check();
      this->tick();

      this->delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = (Scalar*)calloc(ndof, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      coeff_vec_back = (Scalar*)calloc(ndof, sizeof(Scalar));

      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      this->sln_vector = new Scalar[ndof];

      this->on_initialization();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::deinit_solving(Scalar* coeff_vec, Scalar*& coeff_vec_back)
    {
      if(this->delete_coeff_vec)
      {
        ::free(coeff_vec);
        this->delete_coeff_vec = false;
      }

      ::free(coeff_vec_back);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::finalize_solving(Scalar* coeff_vec, Scalar*& coeff_vec_back)
    {
      this->tick();
      this->info("\tNewton: solution duration: %f s.\n", this->last());
      this->on_finish();
      this->deinit_solving(coeff_vec, coeff_vec_back);
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::force_reuse_jacobian_values(unsigned int& successful_steps_with_reused_jacobian)
    {
      int iteration = this->get_parameter_value(p_iteration);
      if(iteration == 1)
        return false;

      double residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 1];
      double previous_residual_norm = this->get_parameter_value(p_residual_norms)[iteration - 2];

#ifdef _DEBUG
      // The following should hold, as if the previous norm was smaller, the step
      // should be restarted with updated damping coefficient.
      assert(previous_residual_norm > residual_norm);
#endif
      
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
    void NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      int ndof = Space<Scalar>::assign_dofs(this->get_spaces());

      Scalar* coeff_vec_back;
      this->init_solving(ndof, coeff_vec, coeff_vec_back);

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

      while (true)
      {
        // User method.
        this->on_step_begin();

        // Assemble just the residual vector.
        this->dp->assemble(coeff_vec, this->residual);

        // Output.
        this->process_vector_output(this->residual, it);
      
        // Current residual norm && current_damping_coefficient.
        residual_norms.push_back(this->calculate_residual_norm());
        solution_norms.push_back(Global<Scalar>::get_l2_norm(coeff_vec, ndof));

        // Initial residual norm.
        if(it == 1)
        {
          this->info("\tNewton: initial residual norm: %g", residual_norms.back());
          solution_change_norm = Global<Scalar>::get_l2_norm(coeff_vec, ndof);
        }
        else
        {
          this->info("\tNewton: iteration %d, residual norm: %g", it - 1, residual_norms.back());
          current_damping_coefficient = this->calculate_damping_coefficient(residual_norm_drop, successful_steps_damping);
          
          // The bad case handling - step has to be restarted.
          if(!residual_norm_drop)
          {
            for (int i = 0; i < ndof; i++)
              coeff_vec[i] = coeff_vec_back[i] + current_damping_coefficient * (coeff_vec[i] - coeff_vec_back[i]);
            if(!this->on_step_end())
            {
              memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
              this->info("Aborted");
              this->finalize_solving(coeff_vec, coeff_vec_back);
              return;
            }
            continue;
          }
        }

        // Find out the state with respect to all residual norms.
        NewtonSolver<Scalar>::ConvergenceState state = get_convergence_state();

        switch(state)
        {
        case Converged:
          // We want to return the solution in a different structure.
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          this->finalize_solving(coeff_vec, coeff_vec_back);
          return;
          break;

        case AboveMaxIterations:
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          throw NewtonException(AboveMaxIterations);
          this->finalize_solving(coeff_vec, coeff_vec_back);
          return;
          break;

        case BelowMinDampingCoeff:
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          throw NewtonException(BelowMinDampingCoeff);
          this->finalize_solving(coeff_vec, coeff_vec_back);
          return;
          break;

        case AboveMaxAllowedResidualNorm:
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          throw NewtonException(AboveMaxAllowedResidualNorm);
          this->finalize_solving(coeff_vec, coeff_vec_back);
          return;
          break;

        case Error:
          throw Exceptions::Exception("Unknown exception in NewtonSolver.");
          this->finalize_solving(coeff_vec, coeff_vec_back);
          return;
          break;

        default:
          // The only state here is NotConverged which yields staying in the loop.
          break;
        }

        // Assemble the jacobian when necessary (nonconstant jacobian, not reusable, ...).
        bool force_reuse_jacobian = this->force_reuse_jacobian_values(successful_steps_jacobian);
        this->conditionally_assemble(coeff_vec, force_reuse_jacobian, false);

        // Output.
        this->process_matrix_output(this->jacobian, it);

        // Multiply the residual vector with -1 since the matrix
        // equation reads J(Y^n) \deltaY^{n + 1} = -F(Y^n).
        this->residual->change_sign();

        // Solve the linear system.
        if(!this->matrix_solver->solve())
        {
          this->deinit_solving(coeff_vec, coeff_vec_back);
          throw Exceptions::LinearMatrixSolverException();
        }

        // Add \deltaY^{n + 1} to Y^n.
        // The good case - step does not have to be restarted.
#ifdef _DEBUG
        assert(residual_norm_drop);
#endif
        memcpy(coeff_vec_back, coeff_vec, sizeof(Scalar)*ndof);
        Scalar* sln_vector = this->matrix_solver->get_sln_vector();
        solution_change_norm = current_damping_coefficient * Global<Scalar>::get_l2_norm(sln_vector, ndof);
        for (int i = 0; i < ndof; i++)
          coeff_vec[i] += current_damping_coefficient * sln_vector[i];

        // User method call.
        if(!this->on_step_end())
        {
          memcpy(this->sln_vector, coeff_vec, ndof * sizeof(Scalar));
          this->info("Aborted");
          this->finalize_solving(coeff_vec, coeff_vec_back);
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