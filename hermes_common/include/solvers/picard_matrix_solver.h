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
#ifndef __HERMES_COMMON_PICARD_MATRIX_SOLVER_H_
#define __HERMES_COMMON_PICARD_MATRIX_SOLVER_H_

#include "solvers/nonlinear_matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    /// See H2D: PicardSolver.
    template<typename Scalar>
    class HERMES_API PicardMatrixSolver : public NonlinearMatrixSolver<Scalar>
    {
    public:
      PicardMatrixSolver();
      virtual ~PicardMatrixSolver() {};

#pragma region anderson-public
      /// Turn on / off the Anderson acceleration. By default it is off.
      void use_Anderson_acceleration(bool to_set);

      /// Set how many last vectors will be used for Anderson acceleration. See the details about the Anderson acceleration for 
      /// explanation of this parameter.
      void set_num_last_vector_used(int num);

      /// Set the Anderson beta coefficient. See the details about the Anderson acceleration for 
      /// explanation of this parameter.
      void set_anderson_beta(double beta);
#pragma endregion

    protected:
      virtual double update_solution_return_change_norm(Scalar* linear_system_solution);

      /// Initialization - called at the beginning of solving.
      virtual void init_solving(Scalar*& coeff_vec);
      
      /// Internal.
      virtual void deinit_solving();

      /// Solve the step's linear system.
      /// Overriden because of Anderson.
      virtual void solve_linear_system();

      /// Norm for convergence.
      virtual double calculate_residual_norm();

      virtual bool damping_factor_condition();

      /// Common constructors code.
      /// Internal setting of default values (see individual set methods).
      void init_picard();

      /// State querying helpers.
      inline std::string getClassName() const { return "PicardMatrixSolver"; }

#pragma region anderson-private
      // Anderson.
      int num_last_vectors_used;
      bool anderson_is_on;
      double anderson_beta;
      /// To store num_last_vectors_used last coefficient vectors.
      Scalar** previous_vectors;
      /// To store num_last_vectors_used - 1 Anderson coefficients.
      Scalar* anderson_coeffs;

      /// Initialization.
      void init_anderson();
      /// Deinitialization.
      void deinit_anderson();

      /// Handle the previous vectors.
      void handle_previous_vectors();
      /// Calcualte the coefficients.
      void calculate_anderson_coeffs();

      /// Previous solution (linear combination coming from Anderson).
      Scalar* previous_Anderson_sln_vector;

      int vec_in_memory;
#pragma endregion
    };
  }
}
#endif