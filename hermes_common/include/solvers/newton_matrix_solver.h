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
/*! \file newton_matrix_solver.h
\brief Newton's method for algebraic equations.
*/
#ifndef __HERMES_COMMON_NEWTON_MATRIX_SOLVER_H_
#define __HERMES_COMMON_NEWTON_MATRIX_SOLVER_H_

#include "solvers/nonlinear_matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    class HERMES_API NewtonMatrixSolver : public NonlinearMatrixSolver<Scalar>
    {
    public:
      NewtonMatrixSolver();
      virtual ~NewtonMatrixSolver() {};

    protected:
      virtual double update_solution_return_change_norm(Scalar* linear_system_solution);
      
      /// Find out the convergence state.
      virtual NonlinearConvergenceState get_convergence_state();

      virtual bool damping_factor_condition();

      /// Common constructors code.
      /// Internal setting of default values (see individual set methods).
      void init_newton();

      /// State querying helpers.
      inline std::string getClassName() const { return "NewtonMatrixSolver"; }
    };
  }
}
#endif
