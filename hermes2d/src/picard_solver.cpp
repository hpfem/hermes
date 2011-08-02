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
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter, Hermes::MatrixSolverType matrix_solver_type) : NonlinearSolver<Scalar>(dp, matrix_solver_type),
      sln_prev_iter(sln_prev_iter)
    {
      if(dp->get_spaces().size() > 1)
        error("PicardSolver so far works only for scalar equations.");
    }
    
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter) : NonlinearSolver<Scalar>(dp, SOLVER_UMFPACK),
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
    bool PicardSolver<Scalar>::solve(double tol, int max_iter)
    {
      int it = 0;
      NewtonSolver<Scalar> newton(static_cast<DiscreteProblem<Scalar>*>(this->dp), 
          this->matrix_solver_type);
      while (true)
      {
        // Delete the old solution vector, if there is any.
        if(this->sln_vector != NULL)
        {
          delete [] this->sln_vector;
          this->sln_vector = NULL;
        }

        // Project the previous solution in order to get a vector of coefficient.
        this->sln_vector = new Scalar[static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs()];
        memset(this->sln_vector, 0, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs() * sizeof(Scalar));

        // Perform Newton's iteration to solve the linear problem.
        // Translate the resulting coefficient vector into the Solution sln_new.
        
        Solution<Scalar> sln_new;
        if (!newton.solve(this->sln_vector, tol, max_iter))
	      {
          warn("Newton's iteration in the Picard's method failed.");
          return false;
        }
        else
          Solution<Scalar>::vector_to_solution(newton.get_sln_vector(),
            static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0), &sln_new);

        double rel_error = Global<Scalar>::calc_rel_error(sln_prev_iter, &sln_new, 
                           HERMES_L2_NORM);

        if (this->verbose_output) 
          info("---- Picard iter %d, ndof %d, rel. error %g%%", 
               it+1, static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs(), 
               rel_error);

        // Stopping criterion.
        if (rel_error < tol)
        {
          for (int i = 0; i < static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_space(0)->get_num_dofs(); i++)
            this->sln_vector[i] = newton.get_sln_vector()[i];
          return true;
        }
        if (it++ >= max_iter)
        {
          if (this->verbose_output) 
            info("Maximum allowed number of Picard iterations exceeded, returning false.");
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
