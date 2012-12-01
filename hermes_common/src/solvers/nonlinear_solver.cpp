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
#include "nonlinear_solver.h"
#include "api.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    NonlinearSolver<Scalar>::NonlinearSolver(DiscreteProblemInterface<Scalar>* dp) : Hermes::Mixins::Loggable(true, NULL), dp(dp), sln_vector(NULL)
    {
    }

    template<typename Scalar>
    NonlinearSolver<Scalar>::~NonlinearSolver()
    {
      if(sln_vector != NULL)
        delete [] sln_vector;
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("NonlinearSolver<Scalar>::solve");
      return;
    }

    template<typename Scalar>
    Scalar *NonlinearSolver<Scalar>::get_sln_vector()
    {
      return sln_vector;
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_iterative_method(const char* iterative_method_name)
    {
      if(Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != SOLVER_AZTECOO)
      {
        this->warn("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      else
      {
        this->iterative_method = (char*)iterative_method_name;
      }
    }

    template<typename Scalar>
    void NonlinearSolver<Scalar>::set_preconditioner(const char* preconditioner_name)
    {
      if(Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType) != SOLVER_AZTECOO)
      {
        this->warn("Trying to set iterative method for a different solver than AztecOO.");
        return;
      }
      else
      {
        this->preconditioner = (char*)preconditioner_name;
      }
    }

    template class HERMES_API NonlinearSolver<double>;
    template class HERMES_API NonlinearSolver<std::complex<double> >;
  }
}