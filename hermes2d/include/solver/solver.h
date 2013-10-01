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
      public Hermes::Hermes2D::Mixins::DiscreteProblemCacheSettings,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>
    {
    public:
      Solver(bool force_use_direct_solver = false);
      Solver(DiscreteProblem<Scalar>* dp, bool force_use_direct_solver = false);
      Solver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space, bool force_use_direct_solver = false);
      Solver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool force_use_direct_solver = false);
      virtual ~Solver();

      /// Basic solve method.
      virtual void solve();

      /// Basic solve method.
      /// \param[in] coeff_vec initiall guess as a vector of coefficients wrt. basis functions.
      virtual void solve(Scalar* coeff_vec) = 0;

      /// Solve.
      /// \param[in] initial_guess Solution to start from (which is projected to obtain the initial coefficient vector.
      virtual void solve(MeshFunctionSharedPtr<Scalar>& initial_guess);

      /// Solve.
      /// \param[in] initial_guess Solutions to start from (which is projected to obtain the initial coefficient vector.
      virtual void solve(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& initial_guess);

      /// See DiscreteProblemCacheSettings in mixins2d.h for details.
      virtual void free_cache();

      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// SettableSpaces helper.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual Hermes::vector<SpaceSharedPtr<Scalar> >& get_spaces();

      /// DiscreteProblemWeakForm helper.
      virtual void set_weak_formulation(WeakForm<Scalar>* wf);

      /// If the cache should not be used for any reason.
      virtual void set_do_not_use_cache(bool to_set = true);
      
      /// \TODO This is not used now.
      /// Report cache hits and misses.
      virtual void set_report_cache_hits_and_misses(bool to_set = true);
    protected:
      virtual bool isOkay() const;
      
      ///< FE problem being solved.
      DiscreteProblem<Scalar>* dp;

      /// This instance owns its DP.
      const bool own_dp;
      
    private:
      void init(bool force_use_direct_solver);
    };
  }
}
#endif