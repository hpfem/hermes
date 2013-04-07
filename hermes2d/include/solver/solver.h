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

namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup userSolvingAPI User solving API
     * \brief Collection of classes that provide the top-level solving capabilities.
    */
    template <typename Scalar>
    class Solver: 
      public Hermes::Mixins::Loggable, 
      public Hermes::Mixins::TimeMeasurable, 
      public Hermes::Mixins::SettableComputationTime, 
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public Hermes::Mixins::OutputAttachable, 
      public Hermes::Hermes2D::Mixins::MatrixRhsOutput<Scalar>, 
      public Hermes::Mixins::IntegrableWithGlobalOrder, 
      public Hermes::Hermes2D::Mixins::StateQueryable, 
      public Hermes::Hermes2D::Mixins::DiscreteProblemCacheSettings,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>
    {
    public:
      Solver();
      Solver(DiscreteProblem<Scalar>* dp);
      Solver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      Solver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual ~Solver();

      void init();

      void keep_matrix_volume_values(int marker, MatrixFormVol<Scalar>* form = NULL);
      void keep_rhs_volume_values(int marker, VectorFormVol<Scalar>* form = NULL);
      void keep_matrix_surface_values(int marker, MatrixFormSurf<Scalar>* form = NULL);
      void keep_rhs_surface_values(int marker, VectorFormSurf<Scalar>* form = NULL);

      virtual bool isOkay() const;
      
      /// See DiscreteProblemCacheSettings in mixins2d.h for details.
      virtual void free_cache();

      /// Return the solution vector.
      virtual Scalar *get_sln_vector();
      
      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// SettableSpaces helpers.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual void set_space(SpaceSharedPtr<Scalar>& space);
      virtual Hermes::vector<SpaceSharedPtr<Scalar> >& get_spaces();

      virtual void set_weak_formulation(WeakForm<Scalar>* wf);

      /// Get the Jacobian.
      SparseMatrix<Scalar>* get_jacobian();

      /// Get the Residual.
      Vector<Scalar>* get_residual();
    protected:
      ///< FE problem being solved.
      DiscreteProblem<Scalar>* dp;

      /// The solution vector.
      Scalar* sln_vector;

      /// Jacobian.
      SparseMatrix<Scalar>* jacobian;

      /// Residual.
      Vector<Scalar>* residual;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* matrix_solver;
      
      /// This instance owns its DP.
      const bool own_dp;
    };
  }
}
#endif
