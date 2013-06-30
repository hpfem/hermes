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
/*! \file solver_newton.h
\brief Newton's method.
*/

#include "solver/newton_solver_convergence_measurement.h"
#include "solver/newton_solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    bool NewtonSolverConvergenceMeasurement<Scalar>::converged(NewtonSolver<Scalar>* newton)
    {
      // get iteration.
      unsigned int iteration = newton->get_parameter_value(newton->iteration());
      if(iteration < 2)
        return false;

      const Hermes::vector<double>& residual_norms = newton->get_parameter_value(newton->residual_norms());
      const Hermes::vector<double>& solution_norms = newton->get_parameter_value(newton->solution_norms());
      const Hermes::vector<double>& solution_change_norms = newton->get_parameter_value(newton->solution_change_norms());

#ifdef _DEBUG
      assert(residual_norms.size() > 1);
      assert(solution_norms.size() > 1);
#endif

      double initial_residual_norm = residual_norms[0];
      double previous_residual_norm = residual_norms[iteration - 2];
      double current_residual_norm = residual_norms[iteration - 1];

      double initial_solution_norm = solution_norms[0];
      double previous_solution_norm = solution_norms[iteration - 2];
      double current_solution_norm = solution_norms[iteration - 1];
      double current_solution_change_norm = solution_change_norms[iteration - 2];

      bool converged;
      if(newton->handleMultipleTolerancesAnd)
        converged = true;
      else
        converged = false;

      double convergence_decision_value[NewtonSolverConvergenceMeasurementTypeCount];
      convergence_decision_value[0] = ((initial_residual_norm - current_residual_norm) / initial_residual_norm);
      convergence_decision_value[1] = ((previous_residual_norm - current_residual_norm) / previous_residual_norm);
      convergence_decision_value[2] = (current_residual_norm / initial_residual_norm);
      convergence_decision_value[3] = (current_residual_norm / previous_residual_norm);
      convergence_decision_value[4] = current_residual_norm;
      convergence_decision_value[5] = current_solution_change_norm;
      convergence_decision_value[6] = (current_solution_change_norm / previous_solution_norm);

      for(int i = 0; i < NewtonSolverConvergenceMeasurementTypeCount; i++)
      {
        bool converged_this_tolerance = (convergence_decision_value[i] < newton->newton_tolerance[i]);
        if(newton->handleMultipleTolerancesAnd)
          converged = converged && converged_this_tolerance;
        else
          if(converged_this_tolerance)
            return true;
      }

      return converged;
    }

    template class HERMES_API NewtonSolverConvergenceMeasurement<double>;
    template class HERMES_API NewtonSolverConvergenceMeasurement<std::complex<double> >;
  }
}