// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://www.hpfem.org/.
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
/*! \file linear_solver.cpp
\brief General linear solver functionality.
*/
#include "solver/linear_solver.h"
#include "solvers/matrix_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(bool force_use_direct_solver) : Solver<Scalar>(false), Hermes::Solvers::MatrixSolver<Scalar>(force_use_direct_solver)
    {
      this->dp = new DiscreteProblem<Scalar>(true);
      this->own_dp = true;
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(DiscreteProblem<Scalar>* dp, bool force_use_direct_solver) : Solver<Scalar>(dp), Hermes::Solvers::MatrixSolver<Scalar>(force_use_direct_solver)
    {
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakFormSharedPtr<Scalar> wf, SpaceSharedPtr<Scalar> space, bool force_use_direct_solver) : Solver<Scalar>(false), Hermes::Solvers::MatrixSolver<Scalar>(force_use_direct_solver)
    {
      this->dp = new DiscreteProblem<Scalar>(wf, space, true);
      this->own_dp = true;
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakFormSharedPtr<Scalar> wf, std::vector<SpaceSharedPtr<Scalar> > spaces, bool force_use_direct_solver) : Solver<Scalar>(false), Hermes::Solvers::MatrixSolver<Scalar>(force_use_direct_solver)
    {
      this->dp = new DiscreteProblem<Scalar>(wf, spaces, true);
      this->own_dp = true;
    }

    template<typename Scalar>
    LinearSolver<Scalar>::~LinearSolver()
    {
    }

    template<typename Scalar>
    Scalar* LinearSolver<Scalar>::get_sln_vector()
    {
      return this->sln_vector;
    }

    template<typename Scalar>
    bool LinearSolver<Scalar>::isOkay() const
    {
      return Solver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_weak_formulation(WeakFormSharedPtr<Scalar> wf)
    {
      Solver<Scalar>::set_weak_formulation(wf);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_spaces(std::vector<SpaceSharedPtr<Scalar> > spaces)
    {
      Solver<Scalar>::set_spaces(spaces);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Solvers::MatrixSolver<Scalar>::set_verbose_output(to_set);
      this->dp->set_verbose_output(to_set);
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->check();

      this->on_initialization();

      this->tick();

      // Extremely important.
      Space<Scalar>::assign_dofs(this->dp->get_spaces());

      // Assemble the residual always and the Matrix when necessary (nonconstant jacobian, not reusable, ...).
      if (this->jacobian_reusable && this->constant_jacobian)
      {
        this->info("\tLinearSolver: assembling... [reusing matrix, assembling rhs].");
        this->dp->assemble(coeff_vec, this->get_residual());
        this->linear_matrix_solver->set_reuse_scheme(Hermes::Solvers::HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
      }
      else
      {
        if (this->jacobian_reusable)
          this->info("\tLinearSolver: assembling... [re-assembling with a reusable matrix structure].");
        else
          this->info("\tLinearSolver: assembling... [assembling the matrix and rhs anew].");
        this->dp->assemble(coeff_vec, this->get_jacobian(), this->get_residual());
        this->linear_matrix_solver->set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
      }

      this->process_matrix_output(this->get_jacobian(), 1);
      this->process_vector_output(this->get_residual(), 1);

      this->tick();
      this->info("\tLinearSolver: assembling done in %s. Solving...", this->last_str().c_str());
      this->tick();

      // Solve, if the solver is iterative, give him the initial guess.
      this->linear_matrix_solver->solve(coeff_vec);

      this->sln_vector = this->linear_matrix_solver->get_sln_vector();

      this->on_finish();

      this->tick();
      this->info("\tLinearSolver: solving done in %s.", this->last_str().c_str());
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}