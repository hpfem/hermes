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
    LinearSolver<Scalar>::LinearSolver() : Solver<Scalar>()
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(DiscreteProblem<Scalar>* dp) : Solver<Scalar>(dp)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver<Scalar>(wf, space)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver<Scalar>(wf, spaces)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::~LinearSolver()
    {
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::init_linear()
    {
      this->dp->set_linear();
    }

    template<typename Scalar>
    bool LinearSolver<Scalar>::isOkay() const
    {
      return this->dp->isOkay();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve()
    {
      this->check();

      this->tick();

      this->on_initialization();

      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      Space<Scalar>::assign_dofs(this->dp->get_spaces());

      // Assemble the residual always and the jacobian when necessary (nonconstant jacobian, not reusable, ...).
      this->conditionally_assemble();
      if(this->report_cache_hits_and_misses)
          this->add_cache_hits_and_misses(this->dp);

      this->process_matrix_output(this->jacobian, 1);
      this->process_vector_output(this->residual, 1);

      this->matrix_solver->solve();

      this->sln_vector = this->matrix_solver->get_sln_vector();

      this->on_finish();
      
      this->tick();
      this->info("\tLinear solver solution duration: %f s.\n", this->last());
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}
