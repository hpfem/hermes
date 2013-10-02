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
/*! \file nonlinear_convergence_measurement.h
\brief nonlinear_convergence_measurement.
*/

#include "solvers/nonlinear_convergence_measurement.h"
#include "solvers/newton_matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    bool NonlinearConvergenceMeasurement<Scalar>::converged(NonlinearMatrixSolver<Scalar>* nonlinear_solver)
    {
      // get iteration.
      unsigned int iteration = nonlinear_solver->get_current_iteration_number();

      const Hermes::vector<double>& residual_norms = nonlinear_solver->get_parameter_value(nonlinear_solver->residual_norms());
      const Hermes::vector<double>& solution_norms = nonlinear_solver->get_parameter_value(nonlinear_solver->solution_norms());
      const Hermes::vector<double>& solution_change_norms = nonlinear_solver->get_parameter_value(nonlinear_solver->solution_change_norms());
      int residual_norms_count = residual_norms.size();
      int solution_norms_count = solution_norms.size();
      int solution_change_norms_count = solution_change_norms.size();

      double initial_residual_norm = residual_norms[0];
      double current_residual_norm = residual_norms.back();
      double previous_residual_norm = iteration == 1 ? current_residual_norm : residual_norms[residual_norms_count - 2];

      double initial_solution_norm = solution_norms[0];
      double previous_solution_norm = solution_norms[solution_norms_count - 2];
      double current_solution_norm = solution_norms.back();
      double current_solution_change_norm = solution_change_norms.back();

      bool converged;
      if(nonlinear_solver->handleMultipleTolerancesAnd)
        converged = true;
      else
        converged = false;

      double convergence_decision_value[NonlinearConvergenceMeasurementTypeCount];
      convergence_decision_value[0] = ((initial_residual_norm - current_residual_norm) / initial_residual_norm);
      convergence_decision_value[1] = ((previous_residual_norm - current_residual_norm) / previous_residual_norm);
      convergence_decision_value[2] = (current_residual_norm / initial_residual_norm);
      convergence_decision_value[3] = (current_residual_norm / previous_residual_norm);
      convergence_decision_value[4] = current_residual_norm;
      convergence_decision_value[5] = current_solution_change_norm;
      convergence_decision_value[6] = (current_solution_change_norm / previous_solution_norm);

      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
      {
        if(!nonlinear_solver->tolerance_set[i])
          continue;

        if(i == 1 && iteration == 1)
        {
          if(nonlinear_solver->handleMultipleTolerancesAnd)
            return false;
          else
            continue;
        }

        bool converged_this_tolerance = (convergence_decision_value[i] < nonlinear_solver->tolerance[i]);
        if(nonlinear_solver->handleMultipleTolerancesAnd)
          converged = converged && converged_this_tolerance;
        else
          if(converged_this_tolerance)
            return true;
      }

      return converged;
    }

    template class HERMES_API NonlinearConvergenceMeasurement<double>;
    template class HERMES_API NonlinearConvergenceMeasurement<std::complex<double> >;
  }
}