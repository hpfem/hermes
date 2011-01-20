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
                                DiscreteProblem* dp, WeakForm* stage_wf_left,
                                WeakForm* stage_wf_right);

/// Takes a matrix M of size ndof times ndof, extends it (formally) to
/// a num_stages*ndof times num_stages*ndof matrix that has M in diagonal blocks and
/// zero everywhere else, and multiplies the new matrix with the vector stage_coeff_vec
/// which has length num_stages*ndof. The result is saved in vector_left which also
/// has length num_stages*ndof.
/// TODO: enable this for other types of matrices.
void HERMES_API multiply_as_diagonal_block_matrix(UMFPackMatrix* matrix_left, int num_stages,
                                                  scalar* stage_coeff_vec, scalar* vector_left);

// Perform one explicit or implicit time step using the Runge-Kutta method
// corresponding to a given Butcher's table. If err_vec != NULL then it will be 
// filled with an error vector calculated using the second B-row of the Butcher's
// table (the second B-row B2 must be nonzero in that case). The negative default 
// values for newton_tol and newton_max_iter are for linear problems.
// Many improvements are needed, a todo list is presented at the beginning of
// the crresponding .cpp file.
bool HERMES_API rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                             scalar* coeff_vec, scalar* err_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver,
                             bool verbose = false, double newton_tol = -1.0, int newton_max_iter = -1,
                             double newton_damping_coeff = 1.0, double newton_max_allowed_residual_norm = 1e6);

// This is a wrapper for the previous function if err_vec is not desired (adaptive time stepping 
// is not attempted). 
bool HERMES_API rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                             scalar* coeff_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver,
                             bool verbose = false, double newton_tol = -1.0, int newton_max_iter = -1,
                             double newton_damping_coeff = 1.0, double newton_max_allowed_residual_norm = 1e6);


#endif
