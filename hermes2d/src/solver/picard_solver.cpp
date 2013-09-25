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
#include "solver/picard_solver.h"
#include "projections/ogprojection.h"
#include "exact_solution.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver() : NonlinearSolver<Scalar>()
    {
      init_picard();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp) : NonlinearSolver<Scalar>(dp)
    {
      init_picard();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : NonlinearSolver<Scalar>(wf, space)
    {
      init_picard();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : NonlinearSolver<Scalar>(wf, spaces)
    {
      init_picard();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::~PicardSolver()
    {
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init_picard()
    {
      num_last_vectors_used = 3;
      anderson_beta = 1.0;
      anderson_is_on = false;
      this->dp->set_linear(false, false);
      this->set_tolerance(1e-3, SolutionChangeRelative);
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::isOkay() const
    {
      if(num_last_vectors_used <= 1)
      {
        throw Hermes::Exceptions::Exception("Picard: Bad number of last iterations to be used (must be at least one).");
        return false;
      }

      return NonlinearSolver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::calculate_anderson_coeffs()
    {
      // If num_last_vectors_used is 2, then there is only one residual, and thus only one alpha coeff which is 1.0.
      if(num_last_vectors_used == 2)
      {
        anderson_coeffs[0] = 1.0;
        return;
      }

      // In the following, num_last_vectors_used is at least three.
      // Thematrix problem will have dimension num_last_vectors_used - 2.
      int n = num_last_vectors_used - 2;

      // Allocate the matrix system for the Anderson coefficients.
      Scalar** mat = new_matrix<Scalar>(n, n);
      Scalar* rhs = new Scalar[n];

      // Set up the matrix and rhs vector.
      for (int i = 0; i < n; i++)
      {
        // Calculate i-th entry of the rhs vector.
        rhs[i] = 0;
        for (int k = 0; k < this->ndof; k++)
        {
          Scalar residual_n_k = previous_vectors[n + 1][k] - previous_vectors[n][k];
          Scalar residual_i_k = previous_vectors[i + 1][k] - previous_vectors[i][k];
          rhs[i] += residual_n_k * (residual_n_k - residual_i_k);
        }
        for (int j = 0; j < n; j++)
        {
          Scalar val = 0;
          for (int k = 0; k < this->ndof; k++)
          {
            Scalar residual_n_k = previous_vectors[n + 1][k] - previous_vectors[n][k];
            Scalar residual_i_k = previous_vectors[i + 1][k] - previous_vectors[i][k];
            Scalar residual_j_k = previous_vectors[j + 1][k] - previous_vectors[j][k];
            val += (residual_n_k - residual_i_k) * (residual_n_k - residual_j_k);
          }

          mat[i][j] = val;
        }
      }

      // Solve the matrix system.
      double d;
      int* perm = new int[n];
      ludcmp(mat, n, perm, &d);
      lubksb<Scalar>(mat, n, perm, rhs);

      // Use the result to define the Anderson coefficients. Remember that
      // n were computed and the last one is 1.0 minus the sum of the 'n' numbers.
      Scalar sum = 0;
      for (int i = 0; i < n; i++)
      {
        anderson_coeffs[i] = rhs[i];
        sum += rhs[i];
      }
      anderson_coeffs[n] = 1.0 - sum;

      // Clean up.
      delete [] mat;
      delete [] rhs;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_num_last_vector_used(int num)
    {
      this->num_last_vectors_used = num;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_anderson_beta(double beta)
    {
      this->anderson_beta = beta;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::use_Anderson_acceleration(bool to_set)
    {
      anderson_is_on = to_set;
    }

    template<typename Scalar>
    int PicardSolver<Scalar>::get_current_iteration_number()
    {
      return this->get_parameter_value(this->p_iteration);
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec)
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
        this->info("\tPicard: done.\n");
        break;
      case AboveMaxIterations:
        throw Exceptions::NonlinearException(AboveMaxIterations);
        break;
      case Error:
        throw Exceptions::Exception("Unknown exception in PicardSolver.");
        break;
      default:
        throw Exceptions::Exception("Unknown ConvergenceState in PicardSolver.");
        break;
      }

      // Return that we should finish.
      return true;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init_anderson()
    {
      if (anderson_is_on) 
      {
        previous_vectors = new Scalar*[num_last_vectors_used];
        for (int i = 0; i < num_last_vectors_used; i++)
          previous_vectors[i] = new Scalar[this->ndof];
        anderson_coeffs = new Scalar[num_last_vectors_used-1];
        memcpy(previous_vectors[0], this->sln_vector, this->ndof*sizeof(Scalar));
      }
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::deinit_anderson()
    {
      if (anderson_is_on)
      {
        for (int i = 0; i < num_last_vectors_used; i++)
          delete [] previous_vectors[i];
        delete [] previous_vectors;
        delete [] anderson_coeffs;
      }
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::handle_previous_vectors(unsigned int& vec_in_memory)
    {
      // If Anderson is used, store the new vector in the memory.
      if (anderson_is_on)
      {
        // If memory not full, just add the vector.
        if (vec_in_memory < num_last_vectors_used)
          memcpy(previous_vectors[vec_in_memory++], this->sln_vector, this->ndof * sizeof(Scalar));
        else
        {
          // If memory full, shift all vectors back, forgetting the oldest one.
          // Save this->sln_vector[] as the newest one.
          Scalar* oldest_vec = previous_vectors[0];

          for (int i = 0; i < num_last_vectors_used-1; i++)
            previous_vectors[i] = previous_vectors[i + 1];

          previous_vectors[num_last_vectors_used-1] = oldest_vec;

          memcpy(oldest_vec, this->sln_vector, this->ndof*sizeof(Scalar));
        }

        if (vec_in_memory >= num_last_vectors_used)
        {
          // Calculate Anderson coefficients.
          this->calculate_anderson_coeffs();

          // Calculate new vector and store it in this->sln_vector[].
          for (int i = 0; i < this->ndof; i++)
          {
            this->sln_vector[i] = 0.;
            for (int j = 1; j < num_last_vectors_used; j++)
              this->sln_vector[i] += anderson_coeffs[j-1] * previous_vectors[j][i] - (1.0 - anderson_beta) * anderson_coeffs[j - 1] * (previous_vectors[j][i] - previous_vectors[j - 1][i]);
          }
        }
      }
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::step_info()
    {
      // Output.
      this->info("\n\tPicard: iteration %d,", this->get_current_iteration_number());
      double sln_change_norm = this->get_parameter_value(this->p_solution_change_norms).back();
      double sln_norm = this->get_parameter_value(this->p_solution_norms).back();
        
      this->info("\n\tPicard: solution change (L2 norm): %g (%g%%).", sln_change_norm, 100. * (sln_change_norm / sln_norm));
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::calculate_error(Scalar* coeff_vec)
    {
      // This is the new sln_vector.
      Scalar* new_sln_vector = this->matrix_solver->get_sln_vector();
     
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(new_sln_vector, this->ndof));

      // sln_vector still stores the old solution.
      // !!!! coeff_vec stores the Anderson-generated previous solution.
      double abs_error = 0.;
      for (int i = 0; i < this->ndof; i++)
        abs_error += std::abs((this->sln_vector[i] - new_sln_vector[i]) * (this->sln_vector[i] - new_sln_vector[i]));
      abs_error = std::sqrt(abs_error);

      this->get_parameter_value(this->p_solution_change_norms).push_back(abs_error);

      // only now we can update the sln_vector.
      memcpy(this->sln_vector, new_sln_vector, sizeof(Scalar)*this->ndof);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init_solving(Scalar*& coeff_vec)
    {
      this->check();
      this->tick();

      // Number of DOFs.
      this->ndof = Space<Scalar>::assign_dofs(this->get_spaces());

      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      this->sln_vector = new Scalar[this->ndof];

      if(coeff_vec == NULL)
        memset(this->sln_vector, 0, this->ndof*sizeof(Scalar));
      else
        memcpy(this->sln_vector, coeff_vec, this->ndof*sizeof(Scalar));

      this->delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = (Scalar*)calloc(this->ndof, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      this->on_initialization();

      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::do_initial_step_return_finished(Scalar* coeff_vec)
    {
      // Store the initial norm.
      this->get_parameter_value(this->p_solution_norms).push_back(get_l2_norm(coeff_vec, this->ndof));
      
      // Solve the linear system.
      this->solve_linear_system(coeff_vec);

      // Calculate errors.
      this->calculate_error(coeff_vec);

      // Use the solution vector for Anderson.
      this->handle_previous_vectors(this->get_parameter_value(this->p_vec_in_memory));

      // Info.
      this->step_info();

      // coeff_vec stores the previous iteration - after this, for the first ordinary step, it will hold the initial step solution.
      memcpy(coeff_vec, this->sln_vector, sizeof(Scalar)*this->ndof);

      // Test convergence - if the first iteration is already a solution.
      if(this->handle_convergence_state_return_finished(this->get_convergence_state(), coeff_vec))
        return true;

      return (this->on_initial_step_end() == false);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve_linear_system(Scalar* coeff_vec)
    {
      // Assemble the residual and also jacobian when necessary (nonconstant jacobian, not reusable, ...).
      this->conditionally_assemble(coeff_vec);
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);

      // Solve, if the solver is iterative, give him the initial guess.
      this->matrix_solver->solve(coeff_vec);
      this->handle_UMFPACK_reports();
      this->process_matrix_output(this->jacobian, this->get_current_iteration_number()); 
      this->process_vector_output(this->residual, this->get_current_iteration_number());

      this->get_parameter_value(this->p_residual_norms).push_back(this->calculate_residual_norm());
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->init_solving(coeff_vec);

      this->init_anderson();

#pragma region parameter_setup
      // Initialize parameters.
      unsigned int vec_in_memory = 1; ///< There is already one vector in the memory.
      unsigned int it = 1;
      Hermes::vector<double> solution_norms;
      Hermes::vector<double> solution_change_norms;
      Hermes::vector<double> residual_norms;

      // Link parameters.
      this->set_parameter_value(this->p_iteration, &it);
      this->set_parameter_value(this->p_vec_in_memory, &vec_in_memory);
      this->set_parameter_value(this->p_residual_norms, &residual_norms);
      this->set_parameter_value(this->p_solution_norms, &solution_norms);
      this->set_parameter_value(this->p_solution_change_norms, &solution_change_norms);
#pragma endregion

      /// Initial iteratios is handled separately (though it is completely identical - this is just to reflect Newton solver).
      if(this->do_initial_step_return_finished(coeff_vec))
      {
        this->info("\tPicard: aborted.");
        // No vector passed, sln_vector in this case contains the solution.
        this->finalize_solving(coeff_vec);
        return;
      }
      else
        it++; 

      while (true)
      {
        // Handle the event of step beginning.
        if(!this->on_step_begin())
        {
          this->info("\tPicard: aborted.");
          this->finalize_solving(coeff_vec);
          return;
        }

        // Solve.
        this->solve_linear_system(coeff_vec);

        // Calculate error.
        this->calculate_error(coeff_vec);

        // Use the solution vector for Anderson.
        this->handle_previous_vectors(this->get_parameter_value(this->p_vec_in_memory));

        // Output for the user.
        this->step_info();

        // Test convergence - if in this iteration we found a solution.
        if(this->handle_convergence_state_return_finished(this->get_convergence_state(), coeff_vec))
          return;

        // Handle the event of end of a step.
        if(!this->on_step_end())
        {
          this->info("Aborted");
          this->finalize_solving(coeff_vec);
          return;
        }

        // Increase counter of iterations.
        it++;

        // Renew the last iteration vector.
        memcpy(coeff_vec, this->sln_vector, this->ndof*sizeof(Scalar));
      }
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::finalize_solving(Scalar* coeff_vec)
    {
      this->tick();
      this->num_iters = this->get_current_iteration_number();
      this->info("\tPicard: solution duration: %f s.\n", this->last());
      this->on_finish();
      this->deinit_solving(coeff_vec);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::deinit_solving(Scalar* coeff_vec)
    {
      if(this->delete_coeff_vec)
      {
        ::free(coeff_vec);
        this->delete_coeff_vec = false;
      }
    }

    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}