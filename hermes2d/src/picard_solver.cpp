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
#include "linear_solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblemLinear<Scalar>* dp, Solution<Scalar>* sln_prev_iter)
        : NonlinearSolver<Scalar>(dp), verbose_output_linear_solver(false)
    {
      int n = slns_prev_iter.size();
      if(dp->get_spaces().size() != n)
        throw Hermes::Exceptions::Exception("Mismatched number of spaces and solutions in PicardSolver.");
      this->slns_prev_iter.push_back(sln_prev_iter);
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblemLinear<Scalar>* dp, Hermes::vector<Solution<Scalar>* > slns_prev_iter)
        : NonlinearSolver<Scalar>(dp), verbose_output_linear_solver(false)
    {
      int n = slns_prev_iter.size();
      if(dp->get_spaces().size() != n)
        throw Hermes::Exceptions::Exception("Mismatched number of spaces and solutions in PicardSolver.");
      for (int i = 0; i<n; i++)
      {
        this->slns_prev_iter.push_back(slns_prev_iter[i]);
      }
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_verbose_output_linear_solver(bool to_set)
    {
      this->verbose_output_linear_solver = to_set;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::solve()
    {
      return solve(1e-8, 100);
    }

    template<typename Scalar>
    void calculate_anderson_coeffs(Scalar** previous_vectors, Scalar* anderson_coeffs, int num_last_vectors_used, int ndof)
    {
      if(num_last_vectors_used <= 1) throw Hermes::Exceptions::Exception("Anderson acceleration makes sense only if at least two last iterations are used.");

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
    bool PicardSolver<Scalar>::solve(double tol, int max_iter, int num_last_vectors_used,
        double anderson_beta)
    {
      // Sanity check.
      if(num_last_vectors_used < 1)
        throw Hermes::Exceptions::Exception("PicardSolver: Bad number of last iterations to be used (must be at least one).");

      // Preliminaries.
      int num_spaces = this->slns_prev_iter.size();
      int ndof = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_num_dofs();
      Hermes::vector<const Space<Scalar>* > spaces = static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
      Hermes::vector<bool> add_dir_lift;
      for(unsigned int i = 0; i < spaces.size(); i++)
        add_dir_lift.push_back(false);
      LinearSolver<Scalar> linear_solver(static_cast<DiscreteProblemLinear<Scalar>*>(this->dp));

      linear_solver.set_verbose_output(this->verbose_output_linear_solver);
      linear_solver.set_verbose_callback(this->get_verbose_callback());

      // Delete solution vector if there is any.
      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      // Project slns_prev_iter on the FE space(s) to obtain initial
      // coefficient vector for the Picard's method.
      this->info("Projecting to obtain initial vector for the Picard's method.");
      this->sln_vector = new Scalar[ndof];

      OGProjection<Scalar> ogProjection;
      ogProjection.project_global(spaces, this->slns_prev_iter, this->sln_vector);

      // Save the coefficient vector, it will be used to calculate increment error
      // after a new coefficient vector is calculated.
      Scalar* last_iter_vector = new Scalar[ndof];
      for (int i = 0; i < ndof; i++)
        last_iter_vector[i] = this->sln_vector[i];

      // Important: This makes the Solution(s) slns_prev_iter compatible with this->sln_vector.
      Solution<Scalar>::vector_to_solutions(this->sln_vector, spaces, this->slns_prev_iter);

      int it = 1;

      while (true)
      {
        linear_solver.solve();

        this->sln_vector = linear_solver.get_sln_vector();

        // Calculate relative error between last_iter_vector[] and this->sln_vector[].
        // FIXME: this is wrong in the complex case (complex conjugation must be used).
        // FIXME: This will crash is norm of last_iter_vector[] is zero.
        double last_iter_vec_norm = 0;
        for (int i = 0; i < ndof; i++)
          last_iter_vec_norm += std::abs(last_iter_vector[i] * last_iter_vector[i]);

        last_iter_vec_norm = sqrt(last_iter_vec_norm);

        double abs_error = 0;
        for (int i = 0; i < ndof; i++) abs_error += std::abs((this->sln_vector[i] - last_iter_vector[i]) * (this->sln_vector[i] - last_iter_vector[i]));
          abs_error = sqrt(abs_error);

        double rel_error = abs_error / last_iter_vec_norm;

        // Output for the user.
        this->info("---- Picard iter %d, ndof %d, rel. error %g%%", it, ndof, rel_error);

        // Stopping because error is sufficiently low.
        if(rel_error < tol)
        {
          delete [] last_iter_vector;
          Solution<Scalar>::vector_to_solutions(this->sln_vector, spaces,  slns_prev_iter);
          return true;
        }

        // Stopping because maximum number of iterations reached.
        if(it >= max_iter)
        {
          this->info("Maximum allowed number of Picard iterations exceeded, returning false.");
          delete [] last_iter_vector;
          // If Anderson acceleration was employed, release memory for the Anderson vectors and coeffs.
          return false;
        }

        // Increase counter of iterations.
        it++;

        // Renew the last iteration vector.
        for (int i = 0; i < ndof; i++)
          last_iter_vector[i] = this->sln_vector[i];

        // Translate the last coefficient vector into previous Solution(s).
        Solution<Scalar>::vector_to_solutions(this->sln_vector, spaces,  slns_prev_iter);
      }
    }
    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}