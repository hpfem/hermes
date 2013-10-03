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
    Solver<Scalar>::Solver(bool force_use_direct_solver) : dp(new DiscreteProblem<Scalar>()), own_dp(true)
    {
      this->init(force_use_direct_solver);
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(DiscreteProblem<Scalar>* dp, bool force_use_direct_solver) :  dp(dp), own_dp(false)
    {
      this->init(force_use_direct_solver);
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space, bool force_use_direct_solver) :  dp(new DiscreteProblem<Scalar>(wf, space)), own_dp(true)
    {
      this->init(force_use_direct_solver);
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool force_use_direct_solver) :  dp(new DiscreteProblem<Scalar>(wf, spaces)), own_dp(true)
    {
      this->init(force_use_direct_solver);
    }
    
    template<typename Scalar>
    Solver<Scalar>::~Solver()
    {
      if(own_dp)
        delete this->dp;
      else
      {
        this->dp->set_matrix(nullptr);
        this->dp->set_rhs(nullptr);
      }
    }

    template<typename Scalar>
    void Solver<Scalar>::init(bool force_use_direct_solver)
    {
    }

    template<typename Scalar>
    void Solver<Scalar>::solve()
    {
      this->solve(nullptr);
    }

    template<typename Scalar>
    void Solver<Scalar>::solve(MeshFunctionSharedPtr<Scalar>& initial_guess)
    {
      if(this->dp->get_spaces().size() != 1)
        throw Hermes::Exceptions::ValueException("dp->get_spaces().size()", this->dp->get_spaces().size(), 1);
      Scalar* coeff_vec = new Scalar[Space<Scalar>::get_num_dofs(this->dp->get_spaces())];
      OGProjection<Scalar>::project_global(this->dp->get_spaces()[0], initial_guess, coeff_vec);
      this->solve(coeff_vec);
      delete [] coeff_vec;
    }

    template<typename Scalar>
    void Solver<Scalar>::solve(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& initial_guess)
    {
      Scalar* coeff_vec = new Scalar[Space<Scalar>::get_num_dofs(this->dp->get_spaces())];
      OGProjection<Scalar>::project_global(this->dp->get_spaces(), initial_guess, coeff_vec);
      this->solve(coeff_vec);
      delete [] coeff_vec;
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
    void Solver<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf)
    {
      this->dp->set_weak_formulation(wf);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_time_step(double time_step)
    {
      this->dp->wf->set_current_time_step(time_step);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      this->dp->set_spaces(spaces);
    }
    
    template<typename Scalar>
    Hermes::vector<SpaceSharedPtr<Scalar> >& Solver<Scalar>::get_spaces()
    {
      return this->dp->get_spaces();
    }

    template<typename Scalar>
    void Solver<Scalar>::set_report_cache_hits_and_misses(bool to_set)
    {
      DiscreteProblemCacheSettings::set_report_cache_hits_and_misses(to_set);
      this->dp->set_report_cache_hits_and_misses(to_set);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_do_not_use_cache(bool to_set)
    {
      DiscreteProblemCacheSettings::set_do_not_use_cache(to_set);
      this->dp->set_do_not_use_cache(to_set);
    }

    template<typename Scalar>
    void Solver<Scalar>::free_cache()
    {
      this->dp->cache.free();
    }

    template class HERMES_API Solver<double>;
    template class HERMES_API Solver<std::complex<double> >;
  }
}
