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
#include "solvers/nonlinear_matrix_solver.h"
#include "common.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    NonlinearMatrixSolver<Scalar>::NonlinearMatrixSolver() : MatrixSolver<Scalar>()
    {
      this->init_nonlinear();
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::init_nonlinear()
    {
      this->handleMultipleTolerancesAnd = false;
      this->max_allowed_iterations = 20;
      this->max_allowed_residual_norm = 1E9;
      this->num_iters = 0;
      this->delete_coeff_vec = false;

      this->clear_tolerances();
    }

    template<typename Scalar>
    NonlinearMatrixSolver<Scalar>::~NonlinearMatrixSolver()
    {
    }

    template<typename Scalar>
    bool NonlinearMatrixSolver<Scalar>::isOkay() const
    {
      bool toleranceSet = false;
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        if(this->tolerance_set[i])
          toleranceSet = true;
      if(!toleranceSet)
      {
        throw Exceptions::Exception("No tolerance set in NonlinearMatrixSolver.");
        return false;
      }
      return true;
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_max_allowed_iterations(int max_allowed_iterations_)
    {
      if(max_allowed_iterations_ < 1)
        throw Exceptions::ValueException("max_allowed_iterations", max_allowed_iterations_, 1);
      this->max_allowed_iterations = max_allowed_iterations_;
    }

    template<typename Scalar>
    double NonlinearMatrixSolver<Scalar>::calculate_residual_norm()
    {
      return get_l2_norm(this->get_residual());
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::init_solving(Scalar*& coeff_vec)
    {
      this->check();
      this->tick();

      // Number of DOFs.
      assert(this->problem_size > 0);

      if(this->sln_vector != NULL)
      {
        delete [] this->sln_vector;
        this->sln_vector = NULL;
      }

      this->sln_vector = new Scalar[this->problem_size];

      if(coeff_vec == NULL)
        memset(this->sln_vector, 0, this->problem_size*sizeof(Scalar));
      else
        memcpy(this->sln_vector, coeff_vec, this->problem_size*sizeof(Scalar));

      this->delete_coeff_vec = false;
      if(coeff_vec == NULL)
      {
        coeff_vec = (Scalar*)calloc(this->problem_size, sizeof(Scalar));
        this->delete_coeff_vec = true;
      }

      this->on_initialization();
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_tolerance(double tolerance_, NonlinearConvergenceMeasurementType toleranceType, bool handleMultipleTolerancesAnd)
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
        throw Exceptions::Exception("Unknown NonlinearConvergenceMeasurementType in NonlinearMatrixSolver::set_tolerance.");
      }
    }

    template<typename Scalar>
    void NonlinearMatrixSolver<Scalar>::set_max_allowed_residual_norm(double max_allowed_residual_norm_to_set)
    {
      if(max_allowed_residual_norm_to_set < 0.0)
        throw Exceptions::ValueException("max_allowed_residual_norm_to_set", max_allowed_residual_norm_to_set, 0.0);
      this->max_allowed_residual_norm = max_allowed_residual_norm_to_set;
    }

    template<typename Scalar>
    int NonlinearMatrixSolver<Scalar>::get_num_iters() const
    {
      return this->num_iters;
    }

    template<typename Scalar>
    NonlinearConvergenceState NonlinearMatrixSolver<Scalar>::get_convergence_state()
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
    void NonlinearMatrixSolver<Scalar>::clear_tolerances()
    {
      for(int i = 0; i < NonlinearConvergenceMeasurementTypeCount; i++)
        this->tolerance[i] = std::numeric_limits<double>::max();
      memset(this->tolerance_set, 0, sizeof(bool)*NonlinearConvergenceMeasurementTypeCount);
    }

    template class HERMES_API NonlinearMatrixSolver<double>;
    template class HERMES_API NonlinearMatrixSolver<std::complex<double> >;
  }

  namespace Exceptions
  {
    NonlinearException::NonlinearException(Solvers::NonlinearConvergenceState convergenceState) : Exception("NonlinearException"), convergenceState(convergenceState)
    {
    }

    Solvers::NonlinearConvergenceState NonlinearException::get_exception_state()
    {
      return this->convergenceState;
    }
  }
}