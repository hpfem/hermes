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

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Solver<Scalar>::Solver() : dp(new DiscreteProblem<Scalar>()), own_dp(true)
    {
      this->init();
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(DiscreteProblem<Scalar>* dp) : dp(dp), own_dp(false)
    {
      this->init();
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : dp(new DiscreteProblem<Scalar>(wf, space)), own_dp(true)
    {
      this->init();
    }

    template<typename Scalar>
    Solver<Scalar>::Solver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : dp(new DiscreteProblem<Scalar>(wf, spaces)), own_dp(true)
    {
      this->init();
    }
    
    template<typename Scalar>
    Solver<Scalar>::~Solver()
    {
      delete jacobian;
      delete residual;
      delete matrix_solver;
      if(own_dp)
        delete this->dp;
      else
      {
        this->dp->set_matrix(NULL);
        this->dp->set_rhs(NULL);
      }
    }

    template<typename Scalar>
    void Solver<Scalar>::init()
    {
      this->jacobian = create_matrix<Scalar>();
      this->residual = create_vector<Scalar>();
      this->matrix_solver = create_linear_solver<Scalar>(this->jacobian, this->residual);
      this->set_verbose_output(true);
      this->sln_vector = NULL;

      this->jacobian_reusable = false;
      this->constant_jacobian = false;
    }

    template<typename Scalar>
    void Solver<Scalar>::set_jacobian_constant(bool to_set)
    {
      this->constant_jacobian = to_set;
      if(!to_set)
        this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void Solver<Scalar>::set_do_not_use_cache(bool to_set)
    {
      DiscreteProblemCacheSettings::set_do_not_use_cache(to_set);
      this->dp->set_do_not_use_cache(to_set);
    }

    template<typename Scalar>
    void Solver<Scalar>::set_report_cache_hits_and_misses(bool to_set)
    {
      DiscreteProblemCacheSettings::set_report_cache_hits_and_misses(to_set);
      this->dp->set_report_cache_hits_and_misses(to_set);
    }
      
    template<typename Scalar>
    void Solver<Scalar>::keep_element_values(int marker, typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side)
    {
      this->dp->selectiveAssembler.state_reuse_kept[dimension][equation_side][marker] = true;
    }

    template<typename Scalar>
    bool Solver<Scalar>::isOkay() const
    {
      if(this->dp->get_weak_formulation() == NULL)
        return false;
      if(this->dp->get_spaces().size() == 0)
        return false;
      return true;
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
      this->jacobian_reusable = false;
    }

    template<typename Scalar>
    void Solver<Scalar>::set_time_step(double time_step)
    {
      this->dp->wf->set_current_time_step(time_step);
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* Solver<Scalar>::get_jacobian()
    {
      return this->jacobian;
    }

    template<typename Scalar>
    Vector<Scalar>* Solver<Scalar>::get_residual()
    {
      return this->residual;
    }

    template<typename Scalar>
    void Solver<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      this->dp->set_spaces(spaces);
      this->jacobian_reusable = false;
    }
    
    template<typename Scalar>
    Hermes::vector<SpaceSharedPtr<Scalar> >& Solver<Scalar>::get_spaces()
    {
      return this->dp->get_spaces();
    }

    template<typename Scalar>
    Scalar *Solver<Scalar>::get_sln_vector()
    {
      return sln_vector;
    }

    template<typename Scalar>
    void Solver<Scalar>::free_cache()
    {
      this->dp->cache.free();
    }

    template<typename Scalar>
    bool Solver<Scalar>::reuse_jacobian_values()
    {
      if(this->constant_jacobian)
        return true;
      else
        return false;
    }

    template<typename Scalar>
    void Solver<Scalar>::conditionally_assemble(Scalar* coeff_vec, bool force_reuse_jacobian_values, bool assemble_residual)
    {
      if(this->jacobian_reusable)
      {
        if(this->reuse_jacobian_values() || force_reuse_jacobian_values)
        {
          if(assemble_residual)
            this->dp->assemble(coeff_vec, this->residual);
          this->matrix_solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
        }
        else
        {
          this->matrix_solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING_AND_SCALING);
          if(assemble_residual)
            this->dp->assemble(coeff_vec, this->jacobian, this->residual);
          else
            this->dp->assemble(coeff_vec, this->jacobian);
        }
      }
      else
      {
        if(assemble_residual)
          this->dp->assemble(coeff_vec, this->jacobian, this->residual);
        else
          this->dp->assemble(coeff_vec, this->jacobian);
        this->matrix_solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
        this->jacobian_reusable = true;
      }
    }

    template class HERMES_API Solver<double>;
    template class HERMES_API Solver<std::complex<double> >;
  }
}
