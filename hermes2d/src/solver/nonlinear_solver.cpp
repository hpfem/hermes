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

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver() : Solver<Scalar>()
    {
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(DiscreteProblem<Scalar>* dp) : Solver<Scalar>(dp)
    {
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver<Scalar>(wf, space)
    {
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver<Scalar>(wf, spaces)
    {
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::~NonlinearSolver()
    {
    }
    
    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_max_allowed_iterations(int max_allowed_iterations_)
    {
      if(max_allowed_iterations_ < 1)
        throw Exceptions::ValueException("max_allowed_iterations", max_allowed_iterations_, 1);
      this->max_allowed_iterations = max_allowed_iterations_;
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

    template class HERMES_API NonlinearSolver<double>;
    template class HERMES_API NonlinearSolver<std::complex<double> >;
  }
}