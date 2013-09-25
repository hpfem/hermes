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
/*! \file nonlinear_solver.cpp
\brief General nonlinear solver functionality.
*/
#include "solver/nonlinear_solver.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver() : Solver<Scalar>()
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(DiscreteProblem<Scalar>* dp) : Solver<Scalar>(dp)
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver<Scalar>(wf, space)
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver<Scalar>(wf, spaces)
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::init_nonlinear()
    {
      this->handleMultipleTolerancesAnd = false;
      this->max_allowed_iterations = 20;
      this->residual_as_function = false;
      this->max_allowed_residual_norm = 1E9;
      this->clear_tolerances();
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::~NonlinearSolver()
    {
    }

    template<typename Scalar>
    bool NonlinearSolver<Scalar>::isOkay() const
    {
      bool toleranceSet = false;
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        if(this->tolerance_set[i])
          toleranceSet = true;
      if(!toleranceSet)
      {
        throw Exceptions::Exception("No tolerance set in NonlinearSolver.");
        return false;
      }
      return Solver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_max_allowed_iterations(int max_allowed_iterations_)
    {
      if(max_allowed_iterations_ < 1)
        throw Exceptions::ValueException("max_allowed_iterations", max_allowed_iterations_, 1);
      this->max_allowed_iterations = max_allowed_iterations_;
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_tolerance(double tolerance_, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
    {
      this->handleMultipleTolerancesAnd = handleMultipleTolerancesAnd;

      if(tolerance_ < 0.0)
        throw Exceptions::ValueException("tolerance", tolerance_, 0.0);

      switch(toleranceType)
      {
      case ResidualNormRelativeToInitial:
        {
          this->tolerance[0] = tolerance_;
          this->tolerance_set[0] = true;
        }
        break;
      case ResidualNormRelativeToPrevious:
        {
          this->tolerance[1] = tolerance_;
          this->tolerance_set[1] = true;
        }
        break;
      case ResidualNormRatioToInitial:
        {
          this->tolerance[2] = tolerance_;
          this->tolerance_set[2] = true;
        }
        break;
      case ResidualNormRatioToPrevious:
        {
          this->tolerance[3] = tolerance_;
          this->tolerance_set[3] = true;
        }
        break;
      case ResidualNormAbsolute:
        {
          this->tolerance[4] = tolerance_;
          this->tolerance_set[4] = true;
        }
        break;
      case SolutionChangeAbsolute:
        {
          this->tolerance[5] = tolerance_;
          this->tolerance_set[5] = true;
        }
        break;
      case SolutionChangeRelative:
        {
          this->tolerance[6] = tolerance_;
          this->tolerance_set[6] = true;
        }
        break;
      default:
        throw Exceptions::Exception("Unknown NonlinearConvergenceMeasurementType in NonlinearSolver::set_tolerance.");
      }
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set)
    {
      if(max_allowed_residual_norm_to_set < 0.0)
        throw Exceptions::ValueException("max_allowed_residual_norm_to_set", max_allowed_residual_norm_to_set, 0.0);
      this->max_allowed_residual_norm = max_allowed_residual_norm_to_set;
    }

    template<typename Scalar>
    int NonlinearSolver<Scalar>::get_num_iters() const
    {
      return this->num_iters;
    }

    template<typename Scalar>
    NonlinearConvergenceState NonlinearSolver<Scalar>::get_convergence_state()
    {
      double residual_norm = this->get_parameter_value(p_residual_norms).back();

      if(residual_norm > this->max_allowed_residual_norm)
        return AboveMaxAllowedResidualNorm;

      if(this->get_current_iteration_number() >= this->max_allowed_iterations)
        return AboveMaxIterations;

      if(NonlinearConvergenceMeasurement<Scalar>::converged(this))
        return Converged;
      else
        return NotConverged;

      return Error;
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::clear_tolerances()
    {
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        this->tolerance[i] = std::numeric_limits<double>::max();
      memset(this->tolerance_set, 0, sizeof(bool)*NonlinearConvergenceMeasurementTypeCount);
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_residual_as_function()
    {
      this->residual_as_function = true;
    }

    template<typename Scalar>
    double NonlinearSolver<Scalar>::calculate_residual_norm()
    {
      // Measure the residual norm.
      if(residual_as_function)
      {
        // Prepare solutions for measuring residual norm.
        Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions;
        Hermes::vector<bool> dir_lift_false;
        for (unsigned int i = 0; i < this->dp->get_spaces().size(); i++) 
        {
          MeshFunctionSharedPtr<Scalar> sharedPtr(new Solution<Scalar>());
          solutions.push_back(sharedPtr);
          dir_lift_false.push_back(false);
        }

        Solution<Scalar>::vector_to_solutions(this->residual, this->dp->get_spaces(), solutions, dir_lift_false);

        // Calculate the norm.

        DefaultNormCalculator<Scalar, HERMES_L2_NORM> normCalculator(solutions.size());
        return normCalculator.calculate_norms(solutions);
      }
      else
      {
        // Calculate the l2-norm of residual vector, this is the traditional way.
        return get_l2_norm(this->residual);
      }
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_iterative_method(const char* iterative_method_name)
    {
#ifdef HAVE_AZTECOO
      if(Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != SOLVER_AZTECOO)
      {
        this->warn("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      else
      {
        this->iterative_method = (char*)iterative_method_name;
        dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar>*>(linear_solver)->set_solver(iterative_method_name);
      }
#else
      this->warn("Trying to set iterative method without AztecOO present.");
#endif
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_preconditioner(const char* preconditioner_name)
    {
#ifdef HAVE_AZTECOO
      if(Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != SOLVER_AZTECOO)
      {
        this->warn("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      else
      {
        dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar> *>(linear_solver)->set_precond(preconditioner_name);
        this->preconditioner = (char*)preconditioner_name;
      }
#else
      this->warn("Trying to set iterative method without AztecOO present.");
#endif
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_preconditioner(Hermes::Preconditioners::Precond<Scalar>* pc)
    {
#ifdef HAVE_AZTECOO
      if(Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != SOLVER_AZTECOO)
      {
        this->warn("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      else
      {
        dynamic_cast<Hermes::Solvers::AztecOOSolver<Scalar> *>(this->matrix_solver)->set_precond(pc);
        this->preconditioner = "Hermes::Preconditioners::Precond";
      }
#else
      this->warn("Trying to set iterative method without AztecOO present.");
#endif
    }

    template class HERMES_API NonlinearSolver<double>;
    template class HERMES_API NonlinearSolver<std::complex<double> >;
  }

  namespace Exceptions
  {
    NonlinearException::NonlinearException(Hermes2D::NonlinearConvergenceState convergenceState) : Exception("NonlinearException"), convergenceState(convergenceState)
    {
    }

    Hermes2D::NonlinearConvergenceState NonlinearException::get_exception_state()
    {
      return this->convergenceState;
    }
  }
}