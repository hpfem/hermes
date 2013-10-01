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
    bool PicardSolver<Scalar>::isOkay() const
    {
      return Solver<Scalar>::isOkay() && Hermes::Solvers::PicardMatrixSolver<Scalar>::isOkay();
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf)
    {
      Solver<Scalar>::set_weak_formulation(wf);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::init_solving(Scalar*& coeff_vec)
    {
      this->problem_size = Space<Scalar>::assign_dofs(this->get_spaces());
      PicardMatrixSolver<Scalar>::init_solving(coeff_vec);
    }

    template<typename Scalar>
    void PicardSolver<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      Solver<Scalar>::set_spaces(spaces);
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_step_end()
    {
      return true;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_initialization()
    {
      return true;
    }

    template<typename Scalar>
    bool PicardSolver<Scalar>::on_initial_step_end()
    {
      return true;
    }

    template class HERMES_API PicardSolver<double>;
    template class HERMES_API PicardSolver<std::complex<double> >;
  }
}