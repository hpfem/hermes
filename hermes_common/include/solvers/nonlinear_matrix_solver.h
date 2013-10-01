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
/*! \file nonlinear_matrix_solver.h
\brief General nonlinear solver functionality.
*/
#ifndef __HERMES_COMMON_NONLINEAR_MATRIX_SOLVER_H_
#define __HERMES_COMMON_NONLINEAR_MATRIX_SOLVER_H_

#include "precond.h"
#include "exceptions.h"
#include "algebra/cs_matrix.h"
#include "algebra/vector.h"
#include "mixins.h"
#include "algebra/algebra_mixins.h"
#include "solvers/nonlinear_convergence_measurement.h"
#include "solvers/matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    /// Convergence measurement strategies.
    /// Count of the types - for solver to hold arrays of such a length.
    const int NonlinearConvergenceMeasurementTypeCount = 7;

    /// This specifies the quantity that is compared to newton_tolerance (settable by set_tolerance()).
    enum NonlinearConvergenceMeasurementType
    {
      ResidualNormRelativeToInitial = 0x0001,
      ResidualNormRelativeToPrevious = 0x0002,
      ResidualNormRatioToInitial = 0x0004,
      ResidualNormRatioToPrevious = 0x0008,
      ResidualNormAbsolute = 0x0010,
      SolutionChangeAbsolute = 0x0020,
      SolutionChangeRelative = 0x0040
    };

    /// Nonlinear Convergence state.
    enum NonlinearConvergenceState
    {
      Converged,
      NotConverged,
      BelowMinDampingCoeff,
      AboveMaxAllowedResidualNorm,
      AboveMaxIterations,
      Error
    };

    template<typename Scalar> class HERMES_API NonlinearConvergenceMeasurement;
    
    /// \brief Base class for defining interface for nonlinear solvers.
    ///
    template <typename Scalar>
    class NonlinearMatrixSolver : 
      public virtual Hermes::Solvers::MatrixSolver<Scalar>,
      public virtual Hermes::Mixins::OutputAttachable,
      public virtual Hermes::Mixins::StateQueryable, 
      public virtual Hermes::Mixins::TimeMeasurable
    {
    public:
      NonlinearMatrixSolver();
      virtual ~NonlinearMatrixSolver();

      virtual void assemble_residual(Scalar* coeff_vec) = 0;
      virtual void assemble_jacobian(Scalar* coeff_vec) = 0;
      virtual void assemble(Scalar* coeff_vec) = 0;
      
      /// Set the maximum number of iterations, thus co-determine when to stop iterations.
      void set_max_allowed_iterations(int max_allowed_iterations);

      /// Clear (reset) all tolerances.
      virtual void clear_tolerances();

      /// Sets the maximum allowed norm of the residual during the calculation.
      /// Default: 1E9
      void set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set);

      /// Find out the convergence state.
      NonlinearConvergenceState get_convergence_state();

      /// Act upon the convergence state.
      /// \return If the main loop in solve() should finalize after this.
      virtual bool handle_convergence_state_return_finished(NonlinearConvergenceState state, Scalar* coeff_vec) = 0;

      /// Initialization - called at the beginning of solving.
      /// Very important e.g. for assigning DOFs before assembling.
      virtual void init_solving(Scalar*& coeff_vec);

      /// Set the residual norm tolerance for ending the Newton's loop.
      /// Default: this->set_tolerance(1e-8, ResidualNormAbsolute);
      /// \param[in] handleMultipleTolerancesAnd If true, multiple tolerances defined will have to be all fulfilled in order to proclaim
      /// solution as a correct one. If false, only one will be enough.
      void set_tolerance(double newton_tol, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd = false);
      
      /// Get the number of iterations.
      int get_num_iters() const;

    protected:
      /// Norm for convergence.
      double calculate_residual_norm();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NonlinearMatrixSolver"; }

      /// Shared code for constructors.
      void init_nonlinear();

      /// Shortcut method for getting the current iteration.
      virtual int get_current_iteration_number() = 0;

      /// Maximum allowed residual norm. If this number is exceeded, the methods solve() return 'false'.
      /// By default set to 1E6.
      /// Possible to change via method set_max_allowed_residual_norm().
      double max_allowed_residual_norm;

      /// Maximum number of iterations allowed.
      int max_allowed_iterations;

      /// There was no initial coefficient vector passed, so this instance had to create one
      /// and this serves as the identificator according to which it will be deleted.
      bool delete_coeff_vec;

      /// Tolerances for all NonlinearConvergenceMeasurementType numbered sequentially as the enum NonlinearConvergenceMeasurementType is.
      double tolerance[NonlinearConvergenceMeasurementTypeCount];

      /// info about set tolerances.
      bool tolerance_set[NonlinearConvergenceMeasurementTypeCount];

      /// If true, multiple tolerances defined will have to be all fulfilled in order to proclaim
      /// solution as a correct one. If false, only one will be enough.
      bool handleMultipleTolerancesAnd;

      int num_iters;

#pragma region OutputAttachable
      // For derived classes - read-only access.
      const OutputParameterDoubleVector& residual_norms() const { return this->p_residual_norms; };
      const OutputParameterDoubleVector& solution_norms() const { return this->p_solution_norms; };
      const OutputParameterDoubleVector& solution_change_norms() const { return this->p_solution_change_norms; };

      // Parameters for OutputAttachable mixin.
      OutputParameterDoubleVector p_residual_norms;
      OutputParameterDoubleVector p_solution_norms;
      OutputParameterDoubleVector p_solution_change_norms;
#pragma endregion
      friend class NonlinearConvergenceMeasurement<Scalar>;
    };
  }

  namespace Exceptions
  {
    class HERMES_API NonlinearException : public Hermes::Exceptions::Exception
    {
    public:
      NonlinearException(Hermes::Solvers::NonlinearConvergenceState convergenceState);

      Hermes::Solvers::NonlinearConvergenceState get_exception_state();

    protected:
      Hermes::Solvers::NonlinearConvergenceState convergenceState;
    };
  }
}
#endif