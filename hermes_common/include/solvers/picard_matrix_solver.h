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
      virtual ~PicardMatrixSolver();

      /// Solve.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual void solve(Scalar* coeff_vec);

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

#pragma region damping-public
      /// Set the damping coefficient Gamma.
      /// Default value = 1.0.
      void set_damping_coefficient(double gamma);
#pragma endregion

    protected:
      /// Common constructors code.
      /// Internal setting of default values (see individual set methods).
      void init_picard();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "PicardMatrixSolver"; }

      /// Init - deinit one solving.
      virtual void init_solving(Scalar*& coeff_vec);
      void deinit_solving(Scalar* coeff_vec);

      /// Finalize solving (+deinit)
      /// For "good" finish.
      void finalize_solving(Scalar* coeff_vec);
      
      /// Calculate and store solution norm and solution change norm.
      void calculate_error(Scalar* coeff_vec);
      
      /// Initial iteratios is handled separately (though it is completely identical - this is just to reflect Newton solver).
      bool do_initial_step_return_finished(Scalar* coeff_vec);

      /// Act upon the convergence state.
      /// \return If the main loop in solve() should finalize after this.
      bool handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec);
      
      void solve_linear_system(Scalar* coeff_vec);

      /// Shortcut method for getting the current iteration.
      int get_current_iteration_number();
      
      /// Output info about the step.
      void step_info();

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
      void handle_previous_vectors(unsigned int& vec_in_memory);
      /// Calcualte the coefficients.
      void calculate_anderson_coeffs();

#pragma endregion

#pragma region damping-private
      /// Damping coefficient Gamma.
      /// Default value = 1.0.
      double gamma;
#pragma endregion

#pragma region OutputAttachable
      // For derived classes - read-only access.
      const OutputParameterUnsignedInt& iteration() const { return this->p_iteration; };
      const OutputParameterUnsignedInt& vec_in_memory() const { return this->p_vec_in_memory; };

    private:
      // Parameters for OutputAttachable mixin.
      OutputParameterUnsignedInt p_iteration;
      OutputParameterUnsignedInt p_vec_in_memory;
#pragma endregion

      friend class NonlinearConvergenceMeasurement<Scalar>;
    };
  }
}
#endif