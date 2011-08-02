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

#include "hermes2d_common_defs.h"
#include "discrete_problem.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Class for Newton's method functions.
    /// \todo Enable more norms than just L2 norm.
    /// \todo Enable more equations.
    template<typename Scalar>
    class HERMES_API PicardSolver : public NonlinearSolver<Scalar>
    {
    public:
      PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter);
      PicardSolver(DiscreteProblem<Scalar>* dp, Solution<Scalar>* sln_prev_iter, Hermes::MatrixSolverType matrix_solver_type);

      /// Solve with default tolerances.
      virtual bool solve();

      /// Solve with user-defined tolerances.
      bool solve(double tol, int max_iter);
    private:
      Solution<Scalar>* sln_prev_iter;
    };
  }
}
#endif

