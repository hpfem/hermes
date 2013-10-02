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

#include "picard_matrix_solver.h"
#include "dense_matrix_operations.h"

using namespace Hermes::Algebra;
using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Solvers
  {

    template<typename Scalar>
    PicardMatrixSolver<Scalar>::PicardMatrixSolver() : NonlinearMatrixSolver<Scalar>()
    {
      init_picard();
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::init_picard()
    {
      this->min_allowed_damping_coeff = 1E-4;
      this->manual_damping = true;
      this->auto_damping_ratio = 2.0;
      this->manual_damping_factor = 1.0;
      this->initial_auto_damping_factor = 1.0;
      this->sufficient_improvement_factor = 1.05;
      this->necessary_successful_steps_to_increase = 3;

      this->sufficient_improvement_factor_jacobian = 1e-1;
      this->max_steps_with_reused_jacobian = 0;

      this->set_tolerance(1e-3, SolutionChangeRelative);

      this->num_last_vectors_used = 3;
      this->anderson_beta = 1.0;
      this->anderson_is_on = false;
      this->vec_in_memory = 0;

      this->use_initial_guess_for_iterative_solvers = true;
    }

    template<typename Scalar>
    double PicardMatrixSolver<Scalar>::update_solution_return_change_norm(Scalar* linear_system_solution)
    {
      double current_damping_factor = this->get_parameter_value(this->p_damping_factors).back();

      double solution_change_norm = 0.;
      for (int i = 0; i < this->problem_size; i++)
      {
        solution_change_norm += std::pow(std::abs(linear_system_solution[i] - this->sln_vector[i]), 2.);
        this->sln_vector[i] += current_damping_factor * (linear_system_solution[i] - this->sln_vector[i]);
      }

      return std::sqrt(solution_change_norm) * current_damping_factor;
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::solve_linear_system()
    {
      NonlinearMatrixSolver<Scalar>::solve_linear_system();
      this->handle_previous_vectors();
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::init_solving(Scalar*& coeff_vec)
    {
      NonlinearMatrixSolver<Scalar>::init_solving(coeff_vec);
      this->init_anderson();
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::deinit_solving()
    {
      this->deinit_anderson();
      NonlinearMatrixSolver<Scalar>::deinit_solving();
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::set_num_last_vector_used(int num)
    {
      this->num_last_vectors_used = num;
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::set_anderson_beta(double beta)
    {
      this->anderson_beta = beta;
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::use_Anderson_acceleration(bool to_set)
    {
      anderson_is_on = to_set;
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::init_anderson()
    {
      if (anderson_is_on) 
      {
        previous_Anderson_sln_vector = new Scalar[this->problem_size];
        previous_vectors = new Scalar*[num_last_vectors_used];
        for (int i = 0; i < num_last_vectors_used; i++)
          previous_vectors[i] = new Scalar[this->problem_size];
        anderson_coeffs = new Scalar[num_last_vectors_used-1];
        memcpy(previous_vectors[0], this->sln_vector, this->problem_size*sizeof(Scalar));
        this->vec_in_memory = 1;
      }
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::deinit_anderson()
    {
      if (anderson_is_on)
      {
        delete [] previous_Anderson_sln_vector;
        for (int i = 0; i < num_last_vectors_used; i++)
          delete [] previous_vectors[i];
        delete [] previous_vectors;
        delete [] anderson_coeffs;
      }
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::handle_previous_vectors()
    {
      // If Anderson is used, store the new vector in the memory.
      if (anderson_is_on)
      {
        // If memory not full, just add the vector.
        if (this->vec_in_memory < num_last_vectors_used)
          memcpy(previous_vectors[this->vec_in_memory++], this->sln_vector, this->problem_size * sizeof(Scalar));
        else
        {
          // If memory full, shift all vectors back, forgetting the oldest one.
          // Save this->sln_vector[] as the newest one.
          Scalar* oldest_vec = previous_vectors[0];

          for (int i = 0; i < num_last_vectors_used-1; i++)
            previous_vectors[i] = previous_vectors[i + 1];

          previous_vectors[num_last_vectors_used-1] = oldest_vec;

          memcpy(oldest_vec, this->sln_vector, this->problem_size*sizeof(Scalar));
        }

        if (this->vec_in_memory >= num_last_vectors_used)
        {
          // Calculate Anderson coefficients.
          this->calculate_anderson_coeffs();

          // Calculate new vector and store it in this->Picard.
          for (int i = 0; i < this->problem_size; i++)
          {
            this->previous_Anderson_sln_vector[i] = 0.;
            for (int j = 1; j < num_last_vectors_used; j++)
              this->previous_Anderson_sln_vector[i] += anderson_coeffs[j-1] * previous_vectors[j][i] - (1.0 - anderson_beta) * anderson_coeffs[j - 1] * (previous_vectors[j][i] - previous_vectors[j - 1][i]);
          }
        }
      }
    }

    template<typename Scalar>
    void PicardMatrixSolver<Scalar>::calculate_anderson_coeffs()
    {
      // If num_last_vectors_used is 2, then there is only one residual, and thus only one alpha coeff which is 1.0.
      if(num_last_vectors_used == 2)
      {
        anderson_coeffs[0] = 1.0;
        return;
      }

      // In the following, num_last_vectors_used is at least three.
      // Thematrix problem will have problem_size num_last_vectors_used - 2.
      int n = num_last_vectors_used - 2;

      // Allocate the matrix system for the Anderson coefficients.
      Scalar** mat = new_matrix<Scalar>(n, n);
      Scalar* rhs = new Scalar[n];
      // Set up the matrix and rhs vector.
      for (int i = 0; i < n; i++)
      {
        // Calculate i-th entry of the rhs vector.
        rhs[i] = 0;
        for (int k = 0; k < this->problem_size; k++)
        {
          Scalar residual_n_k = previous_vectors[n + 1][k] - previous_vectors[n][k];
          Scalar residual_i_k = previous_vectors[i + 1][k] - previous_vectors[i][k];
          rhs[i] += residual_n_k * (residual_n_k - residual_i_k);
        }
        for (int j = 0; j < n; j++)
        {
          Scalar val = 0;
          for (int k = 0; k < this->problem_size; k++)
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

    template class HERMES_API PicardMatrixSolver<double>;
    template class HERMES_API PicardMatrixSolver<std::complex<double> >;
  }
}