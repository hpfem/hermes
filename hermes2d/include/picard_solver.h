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
/*! \file solver_picard.h
\brief Picard's method.
*/
#ifndef __H2D_SOLVER_PICARD_H_
#define __H2D_SOLVER_PICARD_H_

#include "global.h"
#include "projections/ogprojection.h"
#include "discrete_problem.h"
#include "views/scalar_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup userSolvingAPI
    /// Class for the Picard's method.
    template<typename Scalar>
    class HERMES_API PicardSolver : public NonlinearSolver<Scalar>, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, public Hermes::Mixins::OutputAttachable
    {
    public:
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Solution<Scalar>* sln_prev_iter);
      PicardSolver(DiscreteProblemLinear<Scalar>* dp, Hermes::vector<Solution<Scalar>* > slns_prev_iter);
      PicardSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space, Solution<Scalar>* sln_prev_iter);
      PicardSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces, Solution<Scalar>* sln_prev_iter);
      PicardSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space, Hermes::vector<Solution<Scalar>* > slns_prev_iter);
      PicardSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>* > slns_prev_iter);
      ~PicardSolver();
      /// Sets the attribute verbose_output for the inner Newton's loop to the paramater passed.
      void set_verbose_output_linear_solver(bool verbose_output_to_set);

      /// Solve with default tolerances.
      virtual void solve();

      /// set time information for time-dependent problems.
      virtual void setTime(double time);
      virtual void setTimeStep(double timeStep);

      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;
    
      void set_picard_tol(double tol);
      void set_picard_max_iter(int max_iter);
      void set_num_last_vector_used(int num);
      void set_anderson_beta(double beta);
    private:
      void init();
      Hermes::vector<Solution<Scalar>* > slns_prev_iter;
      bool verbose_output_linear_solver;

      /// This instance owns its DP.
      const bool own_dp;

      double tol;
      int max_iter;
      int num_last_vectors_used;
      double beta;
    };
  }
}
#endif
