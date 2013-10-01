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
/*! \file linear_solver.cpp
\brief General linear solver functionality.
*/
#include "solver/linear_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(bool force_use_direct_solver) : Solver<Scalar>(force_use_direct_solver)
    {
      this->init_linear(force_use_direct_solver);
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(DiscreteProblem<Scalar>* dp, bool force_use_direct_solver) : Solver<Scalar>(dp, force_use_direct_solver)
    {
      this->init_linear(force_use_direct_solver);
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space, bool force_use_direct_solver) : Solver<Scalar>(wf, space, force_use_direct_solver)
    {
      this->init_linear(force_use_direct_solver);
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool force_use_direct_solver) : Solver<Scalar>(wf, spaces, force_use_direct_solver)
    {
      this->init_linear(force_use_direct_solver);
    }

    template<typename Scalar>
    LinearSolver<Scalar>::~LinearSolver()
    {
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::init_linear(bool force_use_direct_solver)
    {
      this->dp->set_linear();
    }

    template<typename Scalar>
    bool LinearSolver<Scalar>::isOkay() const
    {
      return Solver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf)
    {
      Solver<Scalar>::set_weak_formulation(wf);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      Solver<Scalar>::set_spaces(spaces);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->check();

      this->tick();

      this->on_initialization();

      this->info("\tLinearSolver: assembling...");

      // Extremely important.
      Space<Scalar>::assign_dofs(this->dp->get_spaces());

      // Assemble the residual always and the Matrix when necessary (nonconstant jacobian, not reusable, ...).
      if(this->jacobian_reusable && this->constant_jacobian)
      {
        this->info("\tLinearSolver: reusing Matrix, assembling RHS.");
        this->dp->assemble(coeff_vec, this->get_residual());
        this->linear_matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
      }
      else
      {
        if(this->jacobian_reusable)
          this->info("\tLinearSolver: recalculating a reusable Matrix.");
        else
          this->info("\tLinearSolver: calculating the Matrix.");
        this->dp->assemble(coeff_vec, this->get_jacobian(), this->get_residual());
        this->linear_matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
      }

      this->process_matrix_output(this->get_jacobian(), 1);
      this->process_vector_output(this->get_residual(), 1);

      this->info("\tLinearSolver: assembling done. Solving...");

      // Solve, if the solver is iterative, give him the initial guess.
      this->linear_matrix_solver->solve(coeff_vec);

      this->sln_vector = this->linear_matrix_solver->get_sln_vector();

      this->on_finish();

      this->tick();
      this->info("\tLinearSolver: done.");
      this->info("\tLinearSolver: solution duration: %f s.", this->last());
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}
