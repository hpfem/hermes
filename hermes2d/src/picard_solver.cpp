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
#include "picard_solver.h"
#include "projections/ogprojection.h"
#include "exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver()
      : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>()), verbose_output_linear_solver(false), own_dp(true)
    {
      init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp)
      : NonlinearSolver<Scalar>(dp), verbose_output_linear_solver(false), own_dp(false)
    {
      init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space)
      : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>(wf, space)), verbose_output_linear_solver(false), own_dp(true)
    {
      init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces)
      : NonlinearSolver<Scalar>(new DiscreteProblem<Scalar>(wf, spaces)), verbose_output_linear_solver(false), own_dp(true)
    {
      init();
    }
    
    template<typename Scalar>
    void PicardSolver<Scalar>::set_weak_formulation(const WeakForm<Scalar>* wf)
    {
      (static_cast<DiscreteProblem<Scalar>*>(this->dp))->set_weak_formulation(wf);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init()
    {
      tol = 1e-4;
      max_iter = 50;
      num_last_vectors_used = 3;
      anderson_beta = 1.0;
      anderson_is_on = false;

      matrix = create_matrix<Scalar>();
      rhs = create_vector<Scalar>();
      linear_solver = create_linear_solver<Scalar>(matrix, rhs);
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::isOkay() const
    {
      if(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_weak_formulation() == NULL)
        return false;
      if(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size() == 0)
        return false;
      return true;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_time(double time)
    {
      Hermes::vector<Space<Scalar>*> spaces;
      for(unsigned int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size(); i++)
        spaces.push_back(const_cast<Space<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(i)));

      Space<Scalar>::update_essential_bc_values(spaces, time);
      const_cast<WeakForm<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->wf)->set_current_time(time);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_time_step(double time_step)
    {
      const_cast<WeakForm<Scalar>*>(static_cast<DiscreteProblem<Scalar>*>(this->dp)->wf)->set_current_time_step(time_step);
    }

    template<typename Scalar>
    PicardSolver<Scalar>::~PicardSolver()
    {
      delete matrix;
      delete rhs;
      delete linear_solver;
      if(own_dp)
        delete this->dp;
      else
        static_cast<DiscreteProblem<Scalar>*>(this->dp)->have_matrix = false;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_spaces(Hermes::vector<const Space<Scalar>*> spaces)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_spaces(spaces);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_space(const Space<Scalar>* space)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_space(space);
    }

    template<typename Scalar>
    Hermes::vector<const Space<Scalar>*> PicardSolver<Scalar>::get_spaces() const
    {
      return static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_verbose_output_linear_solver(bool to_set)
    {
      this->verbose_output_linear_solver = to_set;
    }

    template<typename Scalar>
    void calculate_anderson_coeffs(Scalar** previous_vectors, Scalar* anderson_coeffs, int num_last_vectors_used, int ndof)
    {
      if(num_last_vectors_used <= 1) throw Hermes::Exceptions::Exception("Picard: Anderson acceleration makes sense only if at least two last iterations are used.");

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
        for (int k = 0; k < ndof; k++)
        {
          Scalar residual_n_k = previous_vectors[n + 1][k] - previous_vectors[n][k];
          Scalar residual_i_k = previous_vectors[i + 1][k] - previous_vectors[i][k];
          rhs[i] += residual_n_k * (residual_n_k - residual_i_k);
        }
        for (int j = 0; j < n; j++)
        {
          Scalar val = 0;
          for (int k = 0; k < ndof; k++)
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

      return;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_picard_tol(double tol)
    {
      this->tol = tol;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_picard_max_iter(int max_iter)
    {
      this->max_iter = max_iter;
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
    void PicardSolver<Scalar>::solve(Solution<Scalar>* initial_guess)
    {
      Hermes::vector<Solution<Scalar>*> vectorToPass;
      vectorToPass.push_back(initial_guess);
      this->solve(vectorToPass);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve(Hermes::vector<Solution<Scalar>*> initial_guess)
    {
      int ndof = this->dp->get_num_dofs();
      Scalar* coeff_vec = new Scalar[ndof];
      OGProjection<Scalar> ogProjection;
      ogProjection.project_global(static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces(), initial_guess, coeff_vec);
      this->solve(coeff_vec);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->check();
      this->tick();

      // Sanity check.
      if(num_last_vectors_used < 1)
        throw Hermes::Exceptions::Exception("Picard: Bad number of last iterations to be used (must be at least one).");

      // Preliminaries.
      int num_spaces = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces().size();
      int ndof = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_num_dofs();
      Hermes::vector<const Space<Scalar>* > spaces = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
      Hermes::vector<bool> add_dir_lift;
      for(unsigned int i = 0; i < spaces.size(); i++)
        add_dir_lift.push_back(false);

      // Delete solution vector if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      this->sln_vector = new Scalar[ndof];

      bool delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = new Scalar[ndof];
        memset(coeff_vec, 0, ndof*sizeof(Scalar));
        delete_coeff_vec = true;
      }

      memcpy(this->sln_vector, coeff_vec, ndof*sizeof(Scalar));

      // Save the coefficient vector, it will be used to calculate increment error
      // after a new coefficient vector is calculated.
      Scalar* last_iter_vector = new Scalar[ndof];
      for (int i = 0; i < ndof; i++)
        last_iter_vector[i] = this->sln_vector[i];

      // If Anderson is used, allocate memory for vectors and coefficients.
      Scalar** previous_vectors = NULL;      // To store num_last_vectors_used last coefficient vectors.
      Scalar* anderson_coeffs = NULL;        // To store num_last_vectors_used - 1 Anderson coefficients.
      if (anderson_is_on)
      {
        previous_vectors = new Scalar*[num_last_vectors_used];
        for (int i = 0; i < num_last_vectors_used; i++) previous_vectors[i] = new Scalar[ndof];
        anderson_coeffs = new Scalar[num_last_vectors_used-1];
      }

      // If Anderson is used, save the initial coefficient vector in the memory.
      if (anderson_is_on)
        for (int i = 0; i < ndof; i++) previous_vectors[0][i] = this->sln_vector[i];

      int it = 1;
      int vec_in_memory = 1;   // There is already one vector in the memory.

      this->on_initialization();

      while (true)
      {
        this->on_step_begin();

        this->dp->assemble(last_iter_vector, matrix, rhs);
        if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= it))
        {
          char* fileName = new char[this->matrixFilename.length() + 5];
          if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->matrixFilename.c_str(), it);
          FILE* matrix_file = fopen(fileName, "w+");

          matrix->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
          fclose(matrix_file);
          delete [] fileName;
        }
        if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= it))
        {
          char* fileName = new char[this->RhsFilename.length() + 5];
          if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), it);
          else
            sprintf(fileName, "%s%i", this->RhsFilename.c_str(), it);
          FILE* rhs_file = fopen(fileName, "w+");
          rhs->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          fclose(rhs_file);
          delete [] fileName;
        }

        this->on_step_end();

        //rhs->change_sign();

        // Solve the linear system.
        if(!linear_solver->solve())
          throw Exceptions::LinearMatrixSolverException();

        memcpy(this->sln_vector, linear_solver->get_sln_vector(), sizeof(Scalar)*ndof);

        // If Anderson is used, store the new vector in the memory.
        if (anderson_is_on)
        {
          // If memory not full, just add the vector.
          if (vec_in_memory < num_last_vectors_used)
          {
            for (int i = 0; i < ndof; i++) previous_vectors[vec_in_memory][i] = this->sln_vector[i];
            vec_in_memory++;
          }
          else
          {
            // If memory full, shift all vectors back, forgetting the oldest one.
            // Save this->sln_vector[] as the newest one.
            Scalar* oldest_vec = previous_vectors[0];
            for (int i = 0; i < num_last_vectors_used-1; i++) previous_vectors[i] = previous_vectors[i + 1];
            previous_vectors[num_last_vectors_used-1] = oldest_vec;
            for (int j = 0; j < ndof; j++) previous_vectors[num_last_vectors_used-1][j] = this->sln_vector[j];
          }
        }

        // If there is enough vectors in the memory, calculate Anderson coeffs.
        if (anderson_is_on && vec_in_memory >= num_last_vectors_used)
        {
          // Calculate Anderson coefficients.
          calculate_anderson_coeffs(previous_vectors, anderson_coeffs, num_last_vectors_used, ndof);

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

        // Calculate relative error between last_iter_vector[] and this->sln_vector[].
        // FIXME: this is wrong in the complex case (complex conjugation must be used).
        // FIXME: This will crash if norm of last_iter_vector[] is zero.
        double last_iter_vec_norm = 0;
        for (int i = 0; i < ndof; i++)
          last_iter_vec_norm += std::abs(last_iter_vector[i] * last_iter_vector[i]);

        last_iter_vec_norm = sqrt(last_iter_vec_norm);

        double abs_error = 0;
        for (int i = 0; i < ndof; i++) abs_error += std::abs((this->sln_vector[i] - last_iter_vector[i]) * (this->sln_vector[i] - last_iter_vector[i]));
        abs_error = sqrt(abs_error);

        double rel_error = abs_error / last_iter_vec_norm;

        // Output for the user.
        if(std::abs(last_iter_vec_norm) < 1e-12)
          this->info("\tPicard: iteration %d, nDOFs %d, starting from zero vector.", it, ndof);
        else
          this->info("\tPicard: iteration %d, nDOFs %d, relative error %g%%", it, ndof, rel_error);

        // Stopping because error is sufficiently low.
        if(rel_error < tol)
        {
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (anderson_is_on)
          {
            for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
            delete [] previous_vectors;
            delete [] anderson_coeffs;
          }
          
          static_cast<DiscreteProblem<Scalar>*>(this->dp)->have_matrix = false;

          this->tick();
          this->info("\tPicard: solution duration: %f s.\n", this->last());
          this->on_finish();
          return;
        }

        // Stopping because maximum number of iterations reached.
        if(it >= max_iter)
        {
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (anderson_is_on)
          {
            for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
            delete [] previous_vectors;
            delete [] anderson_coeffs;
          }
          static_cast<DiscreteProblem<Scalar>*>(this->dp)->have_matrix = false;

          this->tick();
          this->info("\tPicard: solution duration: %f s.\n", this->last());

          this->on_finish();
          throw Hermes::Exceptions::Exception("\tPicard: maximum allowed number of Picard iterations exceeded.");
          return;
        }
        this->on_step_end();

        // Increase counter of iterations.
        it++;

        // Renew the last iteration vector.
        for (int i = 0; i < ndof; i++)
          last_iter_vector[i] = this->sln_vector[i];
      }
    }
    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}
