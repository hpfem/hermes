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
/*! \file solver.h
\brief General solver functionality.
*/
#ifndef __H2D_SOLVER_H_
#define __H2D_SOLVER_H_

#include "discrete_problem.h"
#include "global.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup userSolvingAPI User solving API
     * \brief Collection of classes that provide the top-level solving capabilities.
    */

    template <typename Scalar>
    class Solver: 
      public virtual Hermes::Mixins::TimeMeasurable,
      public Hermes::Mixins::SettableComputationTime, 
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public virtual Hermes::Mixins::OutputAttachable,
      public Hermes::Algebra::Mixins::MatrixRhsOutput<Scalar>, 
      public Hermes::Mixins::IntegrableWithGlobalOrder, 
      public virtual Hermes::Mixins::StateQueryable, 
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>
    {
    public:
      Solver(bool initialize_discrete_problem = true);
      Solver(DiscreteProblem<Scalar>* dp);
      Solver(WeakFormSharedPtr<Scalar> wf, SpaceSharedPtr<Scalar> space);
      Solver(WeakFormSharedPtr<Scalar> wf, std::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual ~Solver();

      /// Basic solve method.
      virtual void solve();

      /// Basic solve method.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual void solve(Scalar* coeff_vec) = 0;

      /// Get sln vector.
      virtual Scalar* get_sln_vector() = 0;

      /// Solve.
      /// \param[in] initial_guess Solution to start from (which is projected to obtain the initial coefficient vector.
      virtual void solve(MeshFunctionSharedPtr<Scalar> initial_guess);

      /// Solve.
      /// \param[in] initial_guess Solutions to start from (which is projected to obtain the initial coefficient vector.
      virtual void solve(std::vector<MeshFunctionSharedPtr<Scalar> >& initial_guess);

      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// SettableSpaces helper.
      virtual void set_spaces(std::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual std::vector<SpaceSharedPtr<Scalar> > get_spaces();

      /// DiscreteProblemWeakForm helper.
      virtual void set_weak_formulation(WeakFormSharedPtr<Scalar> wf);

    protected:
      virtual bool isOkay() const;
      
      ///< FE problem being solved.
      DiscreteProblem<Scalar>* dp;

      /// This instance owns its DP.
      bool own_dp;

      template<typename Scalar, typename SolverType> friend class AdaptSolver;
    };
  }
}
#endif