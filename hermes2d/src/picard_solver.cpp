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
#include "newton_solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter, 
          Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type)
    {
      if(dp->get_spaces().size() != 1)
        error("Mismatched number of spaces and solutions in PicardSolver.");
      this->slns_prev_iter.push_back(sln_prev_iter);
    }
    
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Hermes::vector<Solution<Scalar>* > slns_prev_iter, 
          Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type)
    {
      int n = slns_prev_iter.size();
      if(dp->get_spaces().size() != n)
        error("Mismatched number of spaces and solutions in PicardSolver.");
      for (int i=0; i<n; i++) 
      {
        this->slns_prev_iter.push_back(slns_prev_iter[i]);
      }
    }
    
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter) 
           : NonlinearSolver<Scalar>(dp, SOLVER_UMFPACK)
    {
      int n = slns_prev_iter.size();
      if(dp->get_spaces().size() != n)
        error("Mismatched number of spaces and solutions in PicardSolver.");
      this->slns_prev_iter.push_back(sln_prev_iter);
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Hermes::vector<Solution<Scalar>* > slns_prev_iter) 
           : NonlinearSolver<Scalar>(dp, SOLVER_UMFPACK)
    {
      int n = slns_prev_iter.size();
      if(dp->get_spaces().size() != n)
        error("Mismatched number of spaces and solutions in PicardSolver.");
      for (int i=0; i<n; i++) 
      {
        this->slns_prev_iter.push_back(slns_prev_iter[i]);
      }
    }
    
    template<typename Scalar>
    bool PicardSolver<Scalar>::solve()
    {
      return solve(1e-8, 100);
    }

    template<typename Scalar>
    void calculate_anderson_coeffs(Scalar** previous_vectors, Scalar* anderson_coeffs, int num_last_vectors_used, int ndof)
    {
      if (num_last_vectors_used <= 1) error("Anderson acceleration makes sense only if at least two last iterations are used.");

      // If num_last_vectors_used is 2, then there is only one residual, and thus only one alpha coeff which is 1.0.
      if (num_last_vectors_used == 2) 
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
          Scalar residual_n_k = previous_vectors[n+1][k] - previous_vectors[n][k];
          Scalar residual_i_k = previous_vectors[i+1][k] - previous_vectors[i][k];
          rhs[i] += residual_n_k * (residual_n_k - residual_i_k);
	}
        for (int j = 0; j < n; j++)
	{
          Scalar val = 0;
          for (int k = 0; k < ndof; k++) 
  	  { 
            Scalar residual_n_k = previous_vectors[n+1][k] - previous_vectors[n][k];
            Scalar residual_i_k = previous_vectors[i+1][k] - previous_vectors[i][k];
            Scalar residual_j_k = previous_vectors[j+1][k] - previous_vectors[j][k];
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
    bool PicardSolver<Scalar>::solve(double tol, int max_iter, int num_last_vectors_used, 
                                     double anderson_beta)
    {
      // Sanity check.
      if (num_last_vectors_used < 1) 
        error("PicardSolver: Bad number of last iterations to be used (must be at least one).");

      // Preliminaries.
      bool anderson_is_on = (num_last_vectors_used > 1);
      int num_spaces = this->slns_prev_iter.size();
      int ndof = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_num_dofs();
      Hermes::vector<Space<Scalar>* > spaces = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
      NewtonSolver<Scalar> newton(static_cast<DiscreteProblem<Scalar>*>(this->dp), this->matrix_solver_type);
      newton.set_verbose_output(false);

      // Delete solution vector if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      // Project slns_prev_iter on the FE space(s) to obtain initial 
      // coefficient vector for the Picard's method.
      info("Projecting to obtain initial vector for the Picard's method.");
      this->sln_vector = new Scalar[ndof];
      OGProjection<Scalar>::project_global(spaces, this->slns_prev_iter, this->sln_vector, this->matrix_solver_type);

      // Save the coefficient vector, it will be used to calculate increment error
      // after a new coefficient vector is calculated.
      Scalar* last_iter_vector = new Scalar[ndof];
      for (int i = 0; i < ndof; i++) last_iter_vector[i] = this->sln_vector[i];

      // Important: This makes the Solution(s) slns_prev_iter compatible with this->sln_vector.
      Solution<Scalar>::vector_to_solutions(this->sln_vector, spaces, this->slns_prev_iter);

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

      // Picard's loop.
      int it = 1;
      int vec_in_memory = 1;   // There is already one vector in the memory.
      while (true)
      {
        // Perform Newton's iteration to solve the Picard's linear problem.
        if (!newton.solve(this->sln_vector, tol, max_iter))
	{
          warn("Newton's iteration in the Picard's method failed.");
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (anderson_is_on) 
	  {
            for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
            delete [] previous_vectors;
            delete [] anderson_coeffs;
          }
          return false;
        }
        for (int i = 0; i < ndof; i++) this->sln_vector[i] = newton.get_sln_vector()[i];

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
            for (int i = 0; i < num_last_vectors_used-1; i++) previous_vectors[i] = previous_vectors[i+1];
            previous_vectors[num_last_vectors_used-1] = oldest_vec;
            for (int j = 0; j < ndof; j++) previous_vectors[num_last_vectors_used-1][j] = this->sln_vector[j];
          }
        }

        // If there is enough vectors in the memory, calculate Anderson coeffs.
	if (anderson_is_on && vec_in_memory >= num_last_vectors_used)
	{
          // Calculate Anderson coefficients. 
          calculate_anderson_coeffs(previous_vectors, anderson_coeffs, num_last_vectors_used, ndof);

          /*
          // Debug - do not delete.
          printf("alpha: ");
          for (int i = 0; i < num_last_vectors_used-1; i++) printf("%g ", anderson_coeffs[i]);
          printf("\n");
          */

          // Calculate new vector and store it in this->sln_vector[].
          for (int i = 0; i < ndof; i++)
	  {
            this->sln_vector[i] = 0;
            for (int j = 1; j < num_last_vectors_used; j++)
	    {
              this->sln_vector[i] += anderson_coeffs[j-1] * previous_vectors[j][i] 
		- (1.0 - anderson_beta) * anderson_coeffs[j-1] * (previous_vectors[j][i] - previous_vectors[j-1][i]); 
	    }
          }
        }

        // Calculate relative error between last_iter_vector[] and this->sln_vector[].
        // FIXME: this is wrong in the complex case (complex conjugation must be used).
        // FIXME: This will crash is norm of last_iter_vector[] is zero.
        double last_iter_vec_norm = 0;
        for (int i = 0; i < ndof; i++) last_iter_vec_norm += std::abs(last_iter_vector[i] * last_iter_vector[i]);
        last_iter_vec_norm = sqrt(last_iter_vec_norm);
        double abs_error = 0;
        for (int i = 0; i < ndof; i++) abs_error += std::abs((this->sln_vector[i] - last_iter_vector[i]) * 
							     (this->sln_vector[i] - last_iter_vector[i]));
        abs_error = sqrt(abs_error);
        double rel_error = abs_error / last_iter_vec_norm;

        // Output for the user.
        if (this->verbose_output) 
          info("---- Picard iter %d, ndof %d, rel. error %g%%", it, ndof, rel_error);

        // Stopping because error is sufficiently low.
        if (rel_error < tol)
        {
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (anderson_is_on) 
	  {
            for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
            delete [] previous_vectors;
            delete [] anderson_coeffs;
          }
          return true;
        }

        // Stopping because maximum number of iterations reached.
        if (it >= max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Picard iterations exceeded, returning false.");
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (anderson_is_on) 
	  {
            for (int i = 0; i < num_last_vectors_used; i++) delete [] previous_vectors[i];
            delete [] previous_vectors;
            delete [] anderson_coeffs;
          }          
          return false;
        }

        // Increase counter of iterations.
        it++;

        // Renew the last iteration vector.
        for (int i = 0; i < ndof; i++) last_iter_vector[i] = this->sln_vector[i];

        // Translate the last coefficient vector into previous Solution(s).
        Solution<Scalar>::vector_to_solutions(this->sln_vector, spaces,  slns_prev_iter);
      }
    }
    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}
