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

#include "newton_matrix_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    NewtonMatrixSolver<Scalar>::NewtonMatrixSolver() : NonlinearMatrixSolver<Scalar>()
    {
      init_newton();
    }

    template<typename Scalar>
    void NewtonMatrixSolver<Scalar>::init_newton()
    {
      this->min_allowed_damping_coeff = 1E-4;
      this->manual_damping = false;
      this->auto_damping_ratio = 2.0;
      this->initial_auto_damping_factor = 1.0;
      this->sufficient_improvement_factor = 0.95;
      this->necessary_successful_steps_to_increase = 3;

      this->sufficient_improvement_factor_jacobian = 1e-1;
      this->max_steps_with_reused_jacobian = 3;

      this->set_tolerance(1e-8, ResidualNormAbsolute);
    }

    template<typename Scalar>
    NonlinearConvergenceState NewtonMatrixSolver<Scalar>::get_convergence_state()
    {
      double residual_norm = this->get_parameter_value(this->p_residual_norms).back();

      if(residual_norm > this->max_allowed_residual_norm)
        return AboveMaxAllowedResidualNorm;
      else
        return NonlinearMatrixSolver<Scalar>::get_convergence_state();
    }

    template<typename Scalar>
    bool NewtonMatrixSolver<Scalar>::damping_factor_condition()
    {
      double residual_norm = *(this->get_parameter_value(this->residual_norms()).end() - 1);
      double previous_residual_norm = *(this->get_parameter_value(this->residual_norms()).end() - 2);

      return (residual_norm < previous_residual_norm * this->sufficient_improvement_factor);
    }

    template<typename Scalar>
    double NewtonMatrixSolver<Scalar>::update_solution_return_change_norm(Scalar* linear_system_solution)
    {
      double current_damping_factor = this->get_parameter_value(this->p_damping_factors).back();
      
      double solution_change_norm = 0.;
      for (int i = 0; i < this->problem_size; i++)
      {
        solution_change_norm += std::pow(std::abs(linear_system_solution[i]), 2.);
        this->sln_vector[i] += current_damping_factor * linear_system_solution[i];
      }
      return std::sqrt(solution_change_norm) * current_damping_factor;
    }

    template class HERMES_API NewtonMatrixSolver<double>;
    template class HERMES_API NewtonMatrixSolver<std::complex<double> >;
  }
}
