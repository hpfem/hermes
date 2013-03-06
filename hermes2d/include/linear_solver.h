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
/*! \file nonlinear_solver.h
\brief General nonlinear solver functionality.
*/
#ifndef __HERMES_COMMON_LINEAR_SOLVER_H_
#define __HERMES_COMMON_LINEAR_SOLVER_H_

#include "discrete_problem_linear.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup userSolvingAPI User solving API
     * \brief Collection of classes that provide the top-level solving capabilities.
    */

    /// @ingroup userSolvingAPI
    /// Class for solving linear problems.<br>
    /// Typical usage:<br>
    /// // Initialize linear solver.<br>
    /// // Here wf is Hermes2D::WeakForm<double>, space is Hermes2D::Space<double><br>
    /// Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, &space);<br>
    /// <br>
    /// // Solve the linear problem.<br>
    /// try<br>
    /// {<br>
    ///&nbsp;// Just call solve().<br>
    ///&nbsp;linear_solver.solve();<br>
    /// <br>
    ///&nbsp;// Get the solution vector from the solver.<br>
    ///&nbsp;double* sln_vector = linear_solver.get_sln_vector();<br>
    /// <br>
    ///&nbsp;// Translate the solution vector into the previously initialized Solution<double> using the static method vector_to_solution.<br>
    ///&nbsp;Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, &space, &sln);<br>
    /// }<br>
    /// // All kinds of Exceptions may happen (Linear algebraic solver, some bad parameters, some data not initialized...)<br>
    /// catch(Hermes::Exceptions::Exception& e)<br>
    /// {<br>
    ///&nbsp;e.print_msg();<br>
    ///&nbsp;return -1;<br>
    /// }<br>
    /// // For illustrative purposes, otherwise one can just catch std::exception directly, as Hermes::Exceptions::Exception derive from it.<br>
    /// catch(std::exception& e)<br>
    /// {<br>
    ///&nbsp;std::cout << e.what();<br>
    ///&nbsp;return -1;<br>
    /// }<br>
    template <typename Scalar>
    class LinearSolver : public Hermes::Mixins::Loggable, public Hermes::Mixins::TimeMeasurable, public Hermes::Mixins::SettableComputationTime, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, public Hermes::Mixins::OutputAttachable, public Hermes::Hermes2D::Mixins::MatrixRhsOutput<Scalar>, public Hermes::Hermes2D::Mixins::StateQueryable
    {
    public:
      LinearSolver();
      LinearSolver(DiscreteProblemLinear<Scalar>* dp);
      LinearSolver(const WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space);
      LinearSolver(const WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      void init();

      ~LinearSolver();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "LinearSolver"; }

      /// Basic solve method.
      virtual void solve();

      Scalar *get_sln_vector();
      
      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual void set_space(SpaceSharedPtr<Scalar> space);
      virtual Hermes::vector<SpaceSharedPtr<Scalar> > get_spaces() const;
      
      /// Set the weak forms.
      void set_weak_formulation(const WeakForm<Scalar>* wf);

      /// Get the Jacobian.
      SparseMatrix<Scalar>* get_jacobian();

      /// Get the Residual.
      Vector<Scalar>* get_residual();
    protected:
      DiscreteProblemLinear<Scalar>* dp; ///< FE problem being solved.

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
