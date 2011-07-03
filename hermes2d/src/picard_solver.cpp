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
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type)
    {
      if(dp->get_spaces().size() > 1)
        error("PicardSolver so far works only for scalar equations.");
    }
    
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp) : NonlinearSolver<Scalar>(dp, SOLVER_UMFPACK)
    {
      if(dp->get_spaces().size() > 1)
        error("PicardSolver so far works only for scalar equations.");
    }
    
    template<typename Scalar>
    bool PicardSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      return solve(coeff_vec, 1e-8, 100);
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::solve(Scalar* coeff_vec, double tol, int max_iter)
    {
      int it = 0;
      while (true)
      {
        // Translate the previous coefficient vector into a Solution sln_prev_iter.
        Solution<Scalar> sln_prev_iter;
        Solution<Scalar>::vector_to_solution(coeff_vec, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0), &sln_prev_iter);

        // Perform Newton's iteration to solve the linear problem.
        // Translate the resulting coefficient vector into the Solution sln_new.
        NewtonSolver<Scalar> newton(static_cast<DiscreteProblem<Scalar>*>(this->dp), this->matrix_solver_type);
        Solution<Scalar> sln_new;
        if (!newton.solve(coeff_vec, tol, max_iter))
          error("Newton's iteration failed.");
        else
          Solution<Scalar>::vector_to_solution(newton.get_sln_vector(), static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0), &sln_new);

        double rel_error = Global<Scalar>::calc_rel_error(&sln_prev_iter, &sln_new, HERMES_L2_NORM);

        if (this->verbose_output) 
          info("---- Picard iter %d, ndof %d, rel. error %g%%", it+1, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs(), rel_error);

        // Stopping criterion.
        if (rel_error < tol)
          return true;

        if (it++ >= max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Picard iterations exceeded, returning false.");
          return false;
        }
      }
    }
    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}
