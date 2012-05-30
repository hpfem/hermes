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
/*! \file solver_picard.h
\brief Picard's method.
*/
#ifndef __H2D_SOLVER_PICARD_H_
#define __H2D_SOLVER_PICARD_H_

#include "global.h"
#include "projections/ogprojection.h"
#include "discrete_problem.h"
#include "views/scalar_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Class for the Picard's method.
    template<typename Scalar>
    class HERMES_API PicardSolver : public NonlinearSolver<Scalar>
    {
    public:
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Solution<Scalar>* sln_prev_iter);
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Hermes::vector<Solution<Scalar>* > slns_prev_iter);
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Solution<Scalar>* sln_prev_iter,
                   Hermes::MatrixSolverType matrix_solver_type);
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Hermes::vector<Solution<Scalar>* > sln_prev_iter,
                   Hermes::MatrixSolverType matrix_solver_type);
      
      /// Sets the attribute verbose_output for the inner Newton's loop to the paramater passed.
      void set_verbose_output_inner_newton(bool verbose_output_to_set);

      /// Solve with default tolerances.
      virtual bool solve();

      /// Solve with user-defined tolerances.
      /// num_last_vectors_used ... number of last vectors used for Anderson acceleration.
      bool solve(double tol, int max_iter, int num_last_vectors_used = 3, double anderson_beta = 1.0);
    private:
      Hermes::vector<Solution<Scalar>* > slns_prev_iter;
      bool verbose_output_inner_newton;
    };
  }
}
#endif

