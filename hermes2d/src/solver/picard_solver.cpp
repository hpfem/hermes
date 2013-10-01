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
#include "solver/picard_solver.h"
#include "projections/ogprojection.h"
#include "exact_solution.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver() : Solver<Scalar>(), PicardMatrixSolver<Scalar>()
    {
      this->init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(DiscreteProblem<Scalar>* dp) : Solver<Scalar>(dp), PicardMatrixSolver<Scalar>()
    {
      this->init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver<Scalar>(wf, space), PicardMatrixSolver<Scalar>()
    {
      this->init();
    }

    template<typename Scalar>
    PicardSolver<Scalar>::PicardSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver<Scalar>(wf, spaces), PicardMatrixSolver<Scalar>()
    {
      this->init();
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init()
    {
      this->dp->set_linear(false, false);
      assert(this->matrix_solver);
    }

    template<typename Scalar>
    PicardSolver<Scalar>::~PicardSolver()
    {
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      PicardMatrixSolver<Scalar>::solve(coeff_vec);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::assemble_residual(Scalar* coeff_vec)
    {
      this->dp->assemble(coeff_vec, this->get_residual());
      this->process_vector_output(this->get_residual(), this->get_current_iteration_number());
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::assemble_jacobian(Scalar* coeff_vec)
    {
      this->dp->assemble(coeff_vec, this->get_jacobian());
      this->process_matrix_output(this->get_jacobian(), this->get_current_iteration_number()); 
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::assemble(Scalar* coeff_vec)
    {
      this->dp->assemble(coeff_vec, this->get_jacobian(), this->get_residual());
      this->process_vector_output(this->get_residual(), this->get_current_iteration_number());
      this->process_matrix_output(this->get_jacobian(), this->get_current_iteration_number());
    }

    template<typename Scalar>
    int PicardSolver<Scalar>::get_dimension()
    {
      return this->get_jacobian()->get_size();
    }

    template<typename Scalar>
    LinearMatrixSolver<Scalar>* PicardSolver<Scalar>::get_linear_solver()
    {
      return this->matrix_solver;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::isOkay() const
    {
      return Solver<Scalar>::isOkay() && Hermes::Solvers::PicardMatrixSolver<Scalar>::isOkay();
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_step_end()
    {
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);
      this->handle_UMFPACK_reports();

      return true;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_initialization()
    {
      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      // UMFPACK reporting.
      if(this->do_UMFPACK_reporting)
        memset(this->UMFPACK_reporting_data, 0, 3 * sizeof(double));
      return true;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_initial_step_end()
    {
      if(this->report_cache_hits_and_misses)
        this->add_cache_hits_and_misses(this->dp);
      this->handle_UMFPACK_reports();
      return true;
    }

    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}