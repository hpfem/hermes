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

#include "solver/nonlinear_solver.h"
#include "solver/nonlinear_convergence_measurement.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup userSolvingAPI
    /// Class for the Picard's method.<br>
    /// For details about the optionally applied Anderson acceleration, the following website<br>
    /// http://hpfem.org/hermes/hermes-tutorial/doc/_build/html/src/hermes2d/B-nonlinear/01-picard.html
    /// will give an overview.<br>
    /// Typical usage:<br>
    /// // Initialize Picard's solver.<br>
    /// // Here wf is Hermes2D::WeakForm<std::complex<double> >, space is Hermes2D::Space<std::complex<double> ><br>
    /// Hermes::Hermes2D::PicardSolver<std::complex<double> > picard_solver(&wf, &space);<br>
    /// <br>
    /// // Here we have an initial guess for the Picard's method - let us say that we already have a previous time level solution<br>
    /// Solution<std::complex<double> > prevTimeLevelSolution;<br>
    /// <br>
    /// ..... // Here we fill prevTimeLevelSolution.<br>
    /// <br>
    /// // Solve the linear problem.<br>
    /// try<br>
    /// {<br>
    ///&nbsp;// Call solve with the initial guess.<br>
    ///&nbsp;picard_solver.solve(&prevTimeLevelSolution);<br>
    /// <br>
    ///&nbsp;// Get the solution vector from the solver.<br>
    ///&nbsp;std::complex<double> * sln_vector = picard_solver.get_sln_vector();<br>
    /// <br>
    ///&nbsp;// Translate the solution vector into the previously initialized Solution<std::complex<double> > using the static method vector_to_solution.<br>
    ///&nbsp;Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(sln_vector, &space, &sln);<br>
    /// }<br>
    /// // All kinds of Exceptions may happen (Linear algebraic solver, some bad parameters, some data not initialized...)<br>
    /// catch(Hermes::Exceptions::Exception& e)<br>
    /// {<br>
    ///&nbsp;e.print_msg();<br>
    ///&nbsp;return -1;<br>
    /// }<br>
    /// // For illustrative purposes, otherwise one can just catch std::exception directly, as Hermes::Exceptions::Exception derive from it.<br>
    /// catch(std::exception& e)<br>
    /// {<br>
    ///&nbsp;std::cout << e.what();<br>
    ///&nbsp;return -1;<br>
    /// }<br>
    template<typename Scalar>
    class HERMES_API PicardSolver : public Hermes::Hermes2D::NonlinearSolver<Scalar>
    {
    public:
      PicardSolver();
      PicardSolver(DiscreteProblem<Scalar>* dp);
      PicardSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      PicardSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual ~PicardSolver();

      // See the base class for details, the following serves only for avoiding C++ name-hiding.
      using NonlinearSolver<Scalar>::solve;
      
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

    protected:
      /// Common constructors code.
      /// Internal setting of default values (see individual set methods).
      void init_picard();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "PicardSolver"; }

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