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
          Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type),
      sln_prev_iter(sln_prev_iter)
    {
      if(dp->get_spaces().size() > 1)
        error("PicardSolver so far works only for scalar equations.");
    }
    
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter) 
           : NonlinearSolver<Scalar>(dp, SOLVER_UMFPACK),
      sln_prev_iter(sln_prev_iter)
    {
      if(dp->get_spaces().size() > 1)
        error("PicardSolver so far works only for scalar equations.");
    }
    
    template<typename Scalar>
    bool PicardSolver<Scalar>::solve()
    {
      return solve(1e-8, 100);
    }

    template<typename Scalar>
    void calculate_anderson_coeffs(Scalar** anderson_vec, Scalar* anderson_coeffs, int num_last_iter_used, int ndof)
    {
      if (num_last_iter_used <= 1) error("Anderson acceleration makes sense only if at least two last iterations are used.");

      // If num_last_iter_used is 2, then there is only one residual, and thus only one alpha coeff which is 1.0.
      if (num_last_iter_used == 2) 
      {
        anderson_coeffs[0] = 1.0;
        return;
      }

      // In the following, num_last_iter_used is at least three.
      int n = num_last_iter_used - 2;

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
          Scalar residual_n_k = anderson_vec[n+1][k] - anderson_vec[n][k];
          Scalar residual_i_k = anderson_vec[i+1][k] - anderson_vec[i][k];
          rhs[i] += residual_n_k * (residual_n_k - residual_i_k);
	}
        for (int j = 0; j < n; j++)
	{ 
          mat[i][j] = 0;
          for (int k = 0; k < ndof; k++) 
  	  { 
            Scalar residual_n_k = anderson_vec[n+1][k] - anderson_vec[n][k];
            Scalar residual_i_k = anderson_vec[i+1][k] - anderson_vec[i][k];
            Scalar residual_j_k = anderson_vec[j+1][k] - anderson_vec[j][k];
            // FIXME: the following line is wrong in complex case since the last part should be 
            // complex-conjugate. Something like (conj(residual_n_k) - conj(residual_j_k)).
            Scalar val = (residual_n_k - residual_i_k) * (residual_n_k - residual_j_k);
            mat[i][j] += std::abs(val); // In fact 'val' is a real number also in complex case.
	  }
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
    bool PicardSolver<Scalar>::solve(double tol, int max_iter, int num_last_iter_used, 
                                     double anderson_beta)
    {
      if (num_last_iter_used < 1) 
        error("PicardSolver: Bad number of last iterations to be used (must be at least one).");

      int it = 0;
      int ndof = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs();
      Space<Scalar>* space = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0);
      NewtonSolver<Scalar> newton(static_cast<DiscreteProblem<Scalar>*>(this->dp), this->matrix_solver_type);
      newton.set_verbose_output(false);

      // Delete the old solution vector, if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      // Project sln_prev_iter on the FE space to obtain initial 
      // coefficient vector for the Picard's method.
      info("Projecting to obtain initial vector for the Picard's method.");
      this->sln_vector = new Scalar[ndof];
      OGProjection<Scalar>::project_global(space, this->sln_prev_iter, this->sln_vector, this->matrix_solver_type);

      // This makes the Solution sln_prev_iter compatible with this->sln_vector.
      Solution<Scalar>::vector_to_solution(this->sln_vector, space, this->sln_prev_iter);

      /*
      Views::ScalarView<Scalar> view("Projected initial condition", new Views::WinGeom(0, 0, 600, 500));
      view.show(this->sln_prev_iter);
      Views::View::wait();
      */

      // If the number of last iterations used is greater than one, Anderson acceleration will be employed.
      // In this case, allocate memory for Anderson vectors and coeffs.
      Scalar** anderson_vec = NULL;      // To store num_last_iter_used last coefficient vectors.
      Scalar* anderson_coeffs = NULL;    // To store num_last_iter_used - 1 Anderson coefficients.
      if (num_last_iter_used > 1) 
      {
        anderson_vec = new Scalar*[num_last_iter_used];
        for (int i = 0; i < num_last_iter_used; i++) anderson_vec[i] = new Scalar[ndof];
        anderson_coeffs = new Scalar[num_last_iter_used-1];
      }

      // Saving this->sln_vector as the first Anderson vector.
      int anderson_counter = 0;
      if (num_last_iter_used > 1) 
      {
        for (int i = 0; i < ndof; i++) anderson_vec[anderson_counter][i] = this->sln_vector[i];
      }
      anderson_counter++;

      // Picard's loop.
      while (true)
      {
        // Perform Newton's iteration to solve the Picard's linear problem.
        if (!newton.solve(this->sln_vector, tol, max_iter))
	{
          warn("Newton's iteration in the Picard's method failed.");
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (num_last_iter_used > 1) 
	  {
            for (int i = 0; i < num_last_iter_used; i++) delete [] anderson_vec[i];
            delete [] anderson_vec;
            delete [] anderson_coeffs;
          }
          return false;
        }

        // Copy the Newton solution vector into this->sln_vector;
        Scalar* ptr = newton.get_sln_vector();
        for (int i = 0; i < ndof; i++) this->sln_vector[i] = ptr[i];

        if (num_last_iter_used > 1) 
	{
          // Store solution vectors until there is num_last_iter_used of them.
          if (anderson_counter < num_last_iter_used)
  	  {
            for (int i = 0; i < ndof; i++) anderson_vec[anderson_counter][i] = this->sln_vector[i];
            anderson_counter++;
	  }
          else
	  {
            // Save pointer to the oldest vector, move the second to the first, third to the 
            // second etc., and copy the saved oldest vector at the newest position. Then 
            // overwrite it with the new vector.
            Scalar* oldest_vec = anderson_vec[0];
            for (int i = 0; i < num_last_iter_used-1; i++) anderson_vec[i] = anderson_vec[i+1];
            anderson_vec[num_last_iter_used-1] = oldest_vec;
            for (int j = 0; j < ndof; j++) anderson_vec[num_last_iter_used-1][j] = this->sln_vector[j];

            // Calculate Anderson coefficients. 
            calculate_anderson_coeffs(anderson_vec, anderson_coeffs, num_last_iter_used, ndof);
            /*
            printf("Anderson coeffs: ");
            for (int m = 0; m < num_last_iter_used - 1;  m++) printf("%g ", std::abs(anderson_coeffs[m]));
            printf("\n");
            */

            // Calculate new approximation.
            for (int i = 0; i < ndof; i++)
	    {
              this->sln_vector[i] = 0;
              for (int j = 0; j < num_last_iter_used-1; j++)
	      {
                this->sln_vector[i] += anderson_coeffs[j] * anderson_vec[j+1][i]; 
	      }
            }
  	  }
        }

        // Translate the new coefficient vector into a Solution.
        Solution<Scalar> sln_new;
        Solution<Scalar>::vector_to_solution(this->sln_vector, space, &sln_new);

        // Calculate relative error.
        double rel_error = Global<Scalar>::calc_rel_error(sln_prev_iter, &sln_new, HERMES_L2_NORM);

        if (this->verbose_output) 
          info("---- Picard iter %d, ndof %d, rel. error %g%%", it + 1, ndof, rel_error);

        // Stopping criterion.
        if (rel_error < tol)
        {
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (num_last_iter_used > 1) 
	  {
            for (int i = 0; i < num_last_iter_used; i++) delete [] anderson_vec[i];
            delete [] anderson_vec;
            delete [] anderson_coeffs;
          }
          return true;
        }
        if (it++ >= max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Picard iterations exceeded, returning false.");
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          if (num_last_iter_used > 1) 
	  {
            for (int i = 0; i < num_last_iter_used; i++) delete [] anderson_vec[i];
            delete [] anderson_vec;
            delete [] anderson_coeffs;
          }          
          return false;
        }

        // Copy the solution into the previous iteration.
        sln_prev_iter->copy(&sln_new);
      }
    }
    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}
