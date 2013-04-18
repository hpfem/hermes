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

      /// Method setting that the element values on an internal marker will be kept if reusable from a previous solution.
      /// IMPORTANT: so far this method is implemented (works) only on the matrix and volumetric markers.
      /// It is for a discussion if anything else is needed.
      /// \param[in] marker The INTERNAL marker specifying the elements where matrix-vector entries will be reused.
      /// \param[in] dimension Specifying either surface (1d), or volumetric (2d) integrals, thus determining the meaning of marker.
      /// \param[in] equation_side Specifying either matrix or right-hand side.
      void keep_element_values(int marker, typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side);

      /// See DiscreteProblemCacheSettings in mixins2d.h for details.
      virtual void free_cache();

      /// Return the solution vector.
      virtual Scalar *get_sln_vector();
      
      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// SettableSpaces helper.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      virtual Hermes::vector<SpaceSharedPtr<Scalar> >& get_spaces();

      /// DiscreteProblemWeakForm helper.
      virtual void set_weak_formulation(WeakForm<Scalar>* wf);

      /// Sets the jacobian to be constant, i.e. reused whenever possible.
      void set_jacobian_constant(bool to_set = true);

      /// Get the Jacobian.
      SparseMatrix<Scalar>* get_jacobian();

      /// Get the Residual.
      Vector<Scalar>* get_residual();

    protected:
      /// Handle the jacobian re-calculation and re-usage of a previous one.
      void conditionally_assemble(Scalar* coeff_vec = NULL, bool force_reuse_jacobian_values = false, bool assemble_residual = true);

      /// Internal checking.
      virtual bool isOkay() const;
      
      /// Jacobian can be reused if possible.
      bool constant_jacobian;

      /// Jacobian is ready to be reused if desirable.
      bool jacobian_reusable;
      
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

      /// For deciding if the jacobian is constant at this point.
      virtual bool reuse_jacobian_values();

    private:
      void init();
    };
  }
}
#endif