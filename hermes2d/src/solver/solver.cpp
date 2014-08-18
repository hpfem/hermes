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
/*! \file solver.cpp
\brief General solver functionality.
*/
#include "solver/solver.h"
#include "projections/ogprojection.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Solver<Scalar>::Solver(bool initialize_discrete_problem)
    {
      if (initialize_discrete_problem)
      {
        this->dp = new DiscreteProblem<Scalar>();
        own_dp = true;
      }
      else
        own_dp = false;
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(DiscreteProblem<Scalar>* dp) : dp(dp), own_dp(false)
    {
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakFormSharedPtr<Scalar> wf, SpaceSharedPtr<Scalar> space) : dp(new DiscreteProblem<Scalar>(wf, space)), own_dp(true)
    {
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakFormSharedPtr<Scalar> wf, std::vector<SpaceSharedPtr<Scalar> > spaces) : dp(new DiscreteProblem<Scalar>(wf, spaces)), own_dp(true)
    {
    }

    template<typename Scalar>
    Solver<Scalar>::~Solver()
    {
      if (own_dp)
        delete this->dp;
      else
      {
        this->dp->set_matrix(nullptr);
        this->dp->set_rhs(nullptr);
      }
    }

    template<typename Scalar>
    void Solver<Scalar>::solve()
    {
      this->solve(nullptr);
    }

    template<typename Scalar>
    void Solver<Scalar>::solve(MeshFunctionSharedPtr<Scalar> initial_guess)
    {
      if (this->dp->get_spaces().size() != 1)
        throw Hermes::Exceptions::ValueException("dp->get_spaces().size()", this->dp->get_spaces().size(), 1);
      Scalar* coeff_vec = new Scalar[Space<Scalar>::get_num_dofs(this->dp->get_spaces())];
      OGProjection<Scalar>::project_global(this->dp->get_spaces()[0], initial_guess, coeff_vec);
      this->solve(coeff_vec);
      delete[] coeff_vec;
    }

    template<typename Scalar>
    void Solver<Scalar>::solve(std::vector<MeshFunctionSharedPtr<Scalar> > initial_guess)
    {
      Scalar* coeff_vec = new Scalar[Space<Scalar>::get_num_dofs(this->dp->get_spaces())];
      OGProjection<Scalar>::project_global(this->dp->get_spaces(), initial_guess, coeff_vec);
      this->solve(coeff_vec);
      delete[] coeff_vec;
    }

    template<typename Scalar>
    bool Solver<Scalar>::isOkay() const
    {
      return this->dp->isOkay();
    }

    template<typename Scalar>
    void Solver<Scalar>::set_time(double time)
    {
      Space<Scalar>::update_essential_bc_values(this->dp->get_spaces(), time);
      this->dp->wf->set_current_time(time);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_weak_formulation(WeakFormSharedPtr<Scalar> wf)
    {
      this->dp->set_weak_formulation(wf);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_time_step(double time_step)
    {
      this->dp->wf->set_current_time_step(time_step);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_spaces(std::vector<SpaceSharedPtr<Scalar> > spaces)
    {
      this->dp->set_spaces(spaces);
    }

    template<typename Scalar>
    std::vector<SpaceSharedPtr<Scalar> > Solver<Scalar>::get_spaces()
    {
      return this->dp->get_spaces();
    }

    template class HERMES_API Solver < double > ;
    template class HERMES_API Solver < std::complex<double> > ;
  }
}
