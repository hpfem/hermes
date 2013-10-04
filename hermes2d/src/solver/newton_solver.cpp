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
#include "solver/newton_solver.h"
#include "projections/ogprojection.h"
#include "hermes_common.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver() : Solver<Scalar>(), NewtonMatrixSolver<Scalar>()
    {
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(DiscreteProblem<Scalar>* dp) : Solver<Scalar>(dp), NewtonMatrixSolver<Scalar>()
    {
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver<Scalar>(wf, space), NewtonMatrixSolver<Scalar>()
    {
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::NewtonSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver<Scalar>(wf, spaces), NewtonMatrixSolver<Scalar>()
    {
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init()
    {
      this->dp->set_linear(false);
    }

    template<typename Scalar>
    NewtonSolver<Scalar>::~NewtonSolver()
    {
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      NewtonMatrixSolver<Scalar>::solve(coeff_vec);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::assemble_residual(bool store_previous_residual)
    {
      this->dp->assemble(this->sln_vector, this->get_residual());
      this->process_vector_output(this->get_residual(), this->get_current_iteration_number());
      this->get_residual()->change_sign();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::assemble_jacobian(bool store_previous_jacobian)
    {
      this->dp->assemble(this->sln_vector, this->get_jacobian());
      this->process_matrix_output(this->get_jacobian(), this->get_current_iteration_number()); 
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::assemble(bool store_previous_jacobian, bool store_previous_residual)
    {
      this->dp->assemble(this->sln_vector, this->get_jacobian(), this->get_residual());
      this->get_residual()->change_sign();
      this->process_vector_output(this->get_residual(), this->get_current_iteration_number());
      this->process_matrix_output(this->get_jacobian(), this->get_current_iteration_number());
    }

    template<typename Scalar>
    bool NewtonSolver<Scalar>::isOkay() const
    {
      return Solver<Scalar>::isOkay() && Hermes::Solvers::NewtonMatrixSolver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf)
    {
      Solver<Scalar>::set_weak_formulation(wf);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::init_solving(Scalar* coeff_vec)
    {
      this->problem_size = Space<Scalar>::assign_dofs(this->get_spaces());
      NewtonMatrixSolver<Scalar>::init_solving(coeff_vec);
    }

    template<typename Scalar>
    void NewtonSolver<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      Solver<Scalar>::set_spaces(spaces);
      this->jacobian_reusable = false;
    }

    template class HERMES_API NewtonSolver<double>;
    template class HERMES_API NewtonSolver<std::complex<double> >;
  }
}
