// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_RUNGE_KUTTA_H
#define __H2D_RUNGE_KUTTA_H

/// Creates an augmented weak formulation for the multi-stage Runge-Kutta problem.
/// The original discretized equation is M\dot{Y} = F(t, Y) where M is the mass
/// matrix, Y the coefficient vector, and F the (nonlinear) stationary residual.
/// Below, "stage_wf_left" and "stage_wf_right" refer to the left-hand side
/// and right-hand side of the equation, respectively.
void HERMES_API create_stage_wf(double current_time, double time_step, ButcherTable* bt, 
                     DiscreteProblem* dp, WeakForm* stage_wf_right, WeakForm* stage_wf_left);

bool HERMES_API rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                  scalar* coeff_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver,
                  double newton_tol, int newton_max_iter, bool verbose = false,
                  double newton_damping_coeff = 1.0, double newton_max_allowed_residual_norm = 1e6);

#endif
