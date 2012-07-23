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
/*! \file nonlinear_solver.cpp
\brief General nonlinear solver functionality.
*/
#include "linear_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(DiscreteProblemLinear<Scalar>* dp) : dp(dp), sln_vector(NULL), own_dp(false)
    {
      this->init();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space) : dp(new DiscreteProblemLinear<Scalar>(wf, space)), sln_vector(NULL), own_dp(true)
    {
      this->init();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces) : dp(new DiscreteProblemLinear<Scalar>(wf, spaces)), sln_vector(NULL), own_dp(true)
    {
      this->init();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::init()
    {
      this->jacobian = create_matrix<Scalar>();
      this->residual = create_vector<Scalar>();
      this->matrix_solver = create_linear_solver<Scalar>(this->jacobian, this->residual);
    }
    
    template<typename Scalar>
    void LinearSolver<Scalar>::setTime(double time)
    {
      Hermes::vector<Space<Scalar>*> spaces;
      for(unsigned int i = 0; i < this->dp->get_spaces().size(); i++)
        spaces.push_back(const_cast<Space<Scalar>*>(this->dp->get_space(i)));

      Space<Scalar>::update_essential_bc_values(spaces, time);
      const_cast<WeakForm<Scalar>*>(this->dp->wf)->set_current_time(time);
    }
      
    template<typename Scalar>
    void LinearSolver<Scalar>::setTimeStep(double timeStep)
    {
      const_cast<WeakForm<Scalar>*>(this->dp->wf)->set_current_time_step(timeStep);
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_spaces(Hermes::vector<const Space<Scalar>*> spaces)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_spaces(spaces);
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::set_space(const Space<Scalar>* space)
    {
      static_cast<DiscreteProblem<Scalar>*>(this->dp)->set_space(space);
    }
    
    template<typename Scalar>
    Hermes::vector<const Space<Scalar>*> LinearSolver<Scalar>::get_spaces() const
    {
      return static_cast<DiscreteProblem<Scalar>*>(this->dp)->get_spaces();
    }
    
    template<typename Scalar>
    LinearSolver<Scalar>::~LinearSolver()
    {
      delete jacobian;
      delete residual;
      delete matrix_solver;
      if(own_dp)
        delete this->dp;
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve()
    {
      dp->assemble(this->jacobian, this->residual);

      this->matrix_solver->solve();

      this->sln_vector = matrix_solver->get_sln_vector();
    }

    template<typename Scalar>
    Scalar *LinearSolver<Scalar>::get_sln_vector()
    {
      return sln_vector;
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}