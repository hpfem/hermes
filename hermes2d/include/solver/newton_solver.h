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
/*! \file solver_newton.h
\brief Newton's method.
*/
#ifndef __H2D_SOLVER_NEWTON_H_
#define __H2D_SOLVER_NEWTON_H_

#include "solvers/newton_matrix_solver.h"
#include "solver.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup userSolvingAPI
    /// Class for Newton's method.<br>
    /// Typical usage:<br>
    /// // Initialize Newton's solver.<br>
    /// // Here wf is Hermes2D::WeakForm<double>, space is Hermes2D::Space<double><br>
    /// Hermes::Hermes2D::NewtonSolver<double> newton_solver(&wf, &space);<br>
    /// Set a whole bunch of parameters according to your liking.<br>
    /// See the class documentation for all possible parameters.<br>
    /// newton_solver.set_tolerance(1e-6);<br>
    /// newton_solver.set_max_allowed_iterations(15);<br>
    /// newton_solver.set_max_allowed_residual_norm(1e6);<br>
    /// newton_solver.set_min_allowed_damping_coeff(1e-3);<br>
    /// <br>
    /// // Solve the linear problem.<br>
    /// try<br>
    /// {<br>
    ///&nbsp;// Just call solve().<br>
    ///&nbsp;newton_solver.solve();<br>
    /// <br>
    ///&nbsp;// Get the solution vector from the solver.<br>
    ///&nbsp;double* sln_vector = newton_solver.get_sln_vector();<br>
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
    ///&nbsp;std::cout << e.what(); <br>
    ///&nbsp;return -1;<br>
    /// }<br>
    template<typename Scalar>
    class HERMES_API NewtonSolver : 
      public Hermes::Hermes2D::Solver<Scalar>,
      public Hermes::Solvers::NewtonMatrixSolver<Scalar>
    {
    public:
      NewtonSolver();
      NewtonSolver(DiscreteProblem<Scalar>* dp);
      NewtonSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      NewtonSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      void init();
      virtual ~NewtonSolver();

      // See the base class for details, the following serves only for avoiding C++ name-hiding.
      using Solver<Scalar>::solve;
      
      /// Basic solve method - in linear solvers it serves only as an initial guess for iterative solvers.
      /// \param[in] coeff_vec initiall guess.
      virtual void solve(Scalar* coeff_vec);
      
      /// DiscreteProblemWeakForm helper.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);

      /// DiscreteProblemWeakForm helper.
      virtual void set_weak_formulation(WeakForm<Scalar>* wf);

      virtual void assemble_residual(bool store_previous_residual);
      virtual void assemble_jacobian(bool store_previous_jacobian);
      virtual void assemble(bool store_previous_jacobian, bool store_previous_residual);

      /// Initialization - called at the beginning of solving.
      virtual void init_solving(Scalar* coeff_vec);

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "NewtonSolver"; }
    };
  }
}
#endif
