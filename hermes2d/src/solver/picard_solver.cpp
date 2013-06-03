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
      this->picard_tolerance = 1e-4;
      this->max_allowed_iterations = 50;
      num_last_vectors_used = 3;
      anderson_beta = 1.0;
      anderson_is_on = false;
      this->dp->set_linear(false, false);
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::isOkay() const
    {
      if(!NonlinearSolver<Scalar>::isOkay())
        return false;

      if(num_last_vectors_used <= 1)
      {
        throw Hermes::Exceptions::Exception("Picard: Bad number of last iterations to be used (must be at least one).");
        return false;
      }

      return true;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_verbose_output_linear_solver(bool to_set)
    {
      this->verbose_output_linear_solver = to_set;
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
      double** mat = new_matrix<double>(n, n);
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

          // FIXME: This is not a nice way to cast Scalar to double. Not mentioning
          // that this will not work for Scalar = complex.
          double* ptr = (double*)(&val);
          mat[i][j] = *ptr;
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
    void PicardSolver<Scalar>::set_tolerance(double tol)
    {
      this->picard_tolerance = tol;
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
    typename PicardSolver<Scalar>::ConvergenceState PicardSolver<Scalar>::get_convergence_state(double relative_error, int iteration)
    {
      if(iteration >= this->max_allowed_iterations)
        return AboveMaxIterations;

      if(relative_error < this->picard_tolerance && iteration > 1)
        return Converged;
      else
        return NotConverged;

      return Error;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init_anderson()
    {
      if (anderson_is_on) 
      {
        previous_vectors = new Scalar*[num_last_vectors_used];
        for (int i = 0; i < num_last_vectors_used; i++)
          previous_vectors[i] = new Scalar[ndof];
        anderson_coeffs = new Scalar[num_last_vectors_used-1];
        memcpy(previous_vectors[0], this->sln_vector, ndof*sizeof(Scalar));
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
        {
          memcpy(previous_vectors[vec_in_memory], this->sln_vector, ndof*sizeof(Scalar));
          vec_in_memory++;
        }
        else
        {
          // If memory full, shift all vectors back, forgetting the oldest one.
          // Save this->sln_vector[] as the newest one.
          Scalar* oldest_vec = previous_vectors[0];

          for (int i = 0; i < num_last_vectors_used-1; i++)
            previous_vectors[i] = previous_vectors[i + 1];

          previous_vectors[num_last_vectors_used-1] = oldest_vec;

          memcpy(oldest_vec, this->sln_vector, ndof*sizeof(Scalar));

          // Calculate Anderson coefficients.
          calculate_anderson_coeffs();

          // Calculate new vector and store it in this->sln_vector[].
          for (int i = 0; i < ndof; i++)
          {
            this->sln_vector[i] = 0;
            for (int j = 1; j < num_last_vectors_used; j++)
            {
              this->sln_vector[i] += anderson_coeffs[j-1] * previous_vectors[j][i] - (1.0 - anderson_beta) * anderson_coeffs[j-1] * (previous_vectors[j][i] - previous_vectors[j-1][i]);
            }
          }
        }
      }
    }

    template<typename Scalar>
    double PicardSolver<Scalar>::calculate_relative_error(Scalar* coeff_vec)
    {
      double last_iter_vec_norm = get_l2_norm(coeff_vec, ndof);
      if(last_iter_vec_norm < Hermes::epsilon)
      {
        this->warn("\tPicard: a very small error threshold met, the loop should end.");
        return last_iter_vec_norm;
      }

      double abs_error = 0;
      for (int i = 0; i < ndof; i++)
        abs_error += std::abs((this->sln_vector[i] - coeff_vec[i]) * (this->sln_vector[i] - coeff_vec[i]));
      abs_error = sqrt(abs_error);

      double rel_error = abs_error / last_iter_vec_norm;

      return rel_error;
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

      this->sln_vector = new Scalar[ndof];

      if(coeff_vec == NULL)
        memset(this->sln_vector, 0, ndof*sizeof(Scalar));
      else
        memcpy(this->sln_vector, coeff_vec, ndof*sizeof(Scalar));

      this->delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = (Scalar*)calloc(ndof, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      this->on_initialization();

      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->init_solving(coeff_vec);

      this->init_anderson();

      unsigned int it = 1;
      unsigned int vec_in_memory = 1;   // There is already one vector in the memory.
      this->set_parameter_value(this->p_iteration, &it);
      this->set_parameter_value(this->p_vec_in_memory, &vec_in_memory);

      while (true)
      {
        this->on_step_begin();

        // Assemble the residual and also jacobian when necessary (nonconstant jacobian, not reusable, ...).
        this->conditionally_assemble(coeff_vec);
        if(this->report_cache_hits_and_misses)
          this->add_cache_hits_and_misses(this->dp);

        this->process_matrix_output(this->jacobian, it); 
        this->process_vector_output(this->residual, it);

        this->on_step_end();

        // Solve the linear system.
        if(!this->matrix_solver->solve())
          throw Exceptions::LinearMatrixSolverException();

        memcpy(this->sln_vector, this->matrix_solver->get_sln_vector(), sizeof(Scalar)*ndof);

        this->handle_previous_vectors(vec_in_memory);

        if(it > 1)
        {
          double rel_error = this->calculate_relative_error(coeff_vec);

          // Output for the user.
          this->info("\tPicard: iteration %d, nDOFs %d, relative error %g%%", it, ndof, rel_error * 100);

          // Find out the state with respect to all residual norms.
          PicardSolver<Scalar>::ConvergenceState state = get_convergence_state(rel_error, it);

          switch(state)
          {
          case Converged:
            this->deinit_solving(coeff_vec);
            return;
            break;

          case AboveMaxIterations:
            throw Exceptions::ValueException("iterations", it, this->max_allowed_iterations);
            this->deinit_solving(coeff_vec);
            return;
            break;

          case Error:
            throw Exceptions::Exception("Unknown exception in PicardSolver.");
            this->deinit_solving(coeff_vec);
            return;
            break;

          default:
            // The only state here is NotConverged which yields staying in the loop.
            break;
          }
        }

        if(it == 1)
        {
          if(!this->on_initial_step_end())
          {
            this->info("Aborted");
            this->deinit_solving(coeff_vec);
            return;
          }
        }
        else
        {
          if(!this->on_step_end())
          {
            this->info("Aborted");
            this->deinit_solving(coeff_vec);
            return;
          }
        }

        // Increase counter of iterations.
        it++;

        // Renew the last iteration vector.
        memcpy(coeff_vec, this->sln_vector, ndof*sizeof(Scalar));
      }
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
