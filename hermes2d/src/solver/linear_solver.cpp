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
      delete matrix_solver;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::init_linear(bool force_use_direct_solver)
    {
      this->dp->set_linear();

      SparseMatrix<Scalar>* jacobian = create_matrix<Scalar>(force_use_direct_solver);
      Vector<Scalar>* residual = create_vector<Scalar>(force_use_direct_solver);
      this->matrix_solver = create_linear_solver<Scalar>(jacobian, residual, force_use_direct_solver);
    }

    template<typename Scalar>
    bool LinearSolver<Scalar>::isOkay() const
    {
      return this->dp->isOkay();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::conditionally_assemble(Scalar* coeff_vec, bool force_reuse_jacobian_values, bool assemble_residual)
    {
      if(this->jacobian_reusable)
      {
        if(this->constant_jacobian || force_reuse_jacobian_values)
        {
          this->info("\tSolver: reusing Jacobian.");
          if(assemble_residual)
            this->dp->assemble(coeff_vec, this->get_residual());
          this->matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
        }
        else
        {
          this->info("\tSolver: recalculating a reusable Jacobian.");
          this->matrix_solver->set_reuse_scheme(HERMES_REUSE_MATRIX_REORDERING);
          if(assemble_residual)
            this->dp->assemble(coeff_vec, this->get_jacobian(), this->get_residual());
          else
            this->dp->assemble(coeff_vec, this->get_jacobian());
        }
      }
      else
      {
        this->info("\tSolver: Calculating Jacobian.");
        if(assemble_residual)
          this->dp->assemble(coeff_vec, this->get_jacobian(), this->get_residual());
        else
          this->dp->assemble(coeff_vec, this->get_jacobian());
        this->matrix_solver->set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        this->jacobian_reusable = true;
      }
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      this->check();

      this->tick();

      this->on_initialization();

      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      this->info("\tLinear: assembling...");
      Space<Scalar>::assign_dofs(this->dp->get_spaces());

      // Assemble the residual always and the jacobian when necessary (nonconstant jacobian, not reusable, ...).
      this->conditionally_assemble();
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);

      this->process_matrix_output(this->get_jacobian(), 1);
      this->process_vector_output(this->get_residual(), 1);

      this->info("\tLinear: assembling done. Solving...");

      // Solve, if the solver is iterative, give him the initial guess.
      this->matrix_solver->solve(coeff_vec);
      this->handle_UMFPACK_reports();

      this->sln_vector = this->matrix_solver->get_sln_vector();

      this->on_finish();

      this->tick();
      this->info("\tLinear: done.");
      this->info("\tLinear: solution duration: %f s.", this->last());
    }

    template<typename Scalar>
    LinearMatrixSolver<Scalar>* LinearSolver<Scalar>::get_linear_solver()
    {
      return this->matrix_solver;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* LinearSolver<Scalar>::get_jacobian()
    {
      return this->matrix_solver->get_matrix();
    }

    template<typename Scalar>
    Vector<Scalar>* LinearSolver<Scalar>::get_residual()
    {
      return this->matrix_solver->get_rhs();
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}
