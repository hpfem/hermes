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
    /// Class for the Picard's method.<br>
    /// For details about the optionally applied Anderson acceleration, the following website<br>
    /// http://hpfem.org/hermes/hermes-tutorial/doc/_build/html/src/hermes2d/B-nonlinear/01-picard.html
    /// will give an overview.<br>
    /// Typical usage:<br>
    /// // Initialize Picard's solver.<br>
    /// // Here wf is Hermes2D::WeakForm<std::complex<double> >, space is Hermes2D::Space<std::complex<double> ><br>
    /// Hermes::Hermes2D::PicardSolver<std::complex<double> > picard_solver(&wf, &space);<br>
    /// <br>
    /// // Here we have an initial guess for the Picard's method - let us say that we already have a previous time level solution<br>
    /// Solution<std::complex<double> > prevTimeLevelSolution;<br>
    /// <br>
    /// ..... // Here we fill prevTimeLevelSolution.<br>
    /// <br>
    /// // Solve the linear problem.<br>
    /// try<br>
    /// {<br>
    ///&nbsp;// Call solve with the initial guess.<br>
    ///&nbsp;picard_solver.solve(&prevTimeLevelSolution);<br>
    /// <br>
    ///&nbsp;// Get the solution vector from the solver.<br>
    ///&nbsp;std::complex<double> * sln_vector = picard_solver.get_sln_vector();<br>
    /// <br>
    ///&nbsp;// Translate the solution vector into the previously initialized Solution<std::complex<double> > using the static method vector_to_solution.<br>
    ///&nbsp;Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(sln_vector, &space, &sln);<br>
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
    template<typename Scalar>
    class HERMES_API PicardSolver : public Solvers::NonlinearSolver<Scalar>, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, public Hermes::Mixins::OutputAttachable, public Hermes::Hermes2D::Mixins::MatrixRhsOutput<Scalar>, public Hermes::Hermes2D::Mixins::StateQueryable
    {
    public:
      PicardSolver();
      PicardSolver(DiscreteProblem<Scalar>* dp);
      PicardSolver(const WeakForm<Scalar>* wf, const Space<Scalar>* space);
      PicardSolver(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar>*> spaces);
      ~PicardSolver();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "PicardSolver"; }

      /// Sets the attribute verbose_output for the inner Newton's loop to the paramater passed.
      void set_verbose_output_linear_solver(bool verbose_output_to_set);

      /// Solve.
      /// \param[in] coeff_vec Ceofficient vector to start from.
      virtual void solve(Scalar* coeff_vec = NULL);

      /// Solve.
      /// \param[in] initial_guess Solution to start from (which is projected to obtain the initial coefficient vector.
      void solve(Solution<Scalar>* initial_guess);

      /// Solve.
      /// \param[in] initial_guess Solutions to start from (which is projected to obtain the initial coefficient vector.
      void solve(Hermes::vector<Solution<Scalar>*> initial_guess);

      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// Overridden Mixins::SettableSpaces methods.
      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;

      /// Turn on / off the Anderson acceleration. By default it is off.
      void use_Anderson_acceleration(bool to_set);
    
      /// Set the relative tolerance, thus co-determine when to stop Picard's iterations.
      void set_picard_tol(double tol);
      /// Set the maximum number of Picard's iterations, thus co-determine when to stop Picard's iterations.
      void set_picard_max_iter(int max_iter);
      /// Set how many last vectors will be used for Anderson acceleration. See the details about the Anderson acceleration for 
      /// explanation of this parameter.
      void set_num_last_vector_used(int num);
      /// Set the Anderson beta coefficient. See the details about the Anderson acceleration for 
      /// explanation of this parameter.
      void set_anderson_beta(double beta);

      /// Set the weak forms.
      void set_weak_formulation(const WeakForm<Scalar>* wf);
    protected:
      void init();
      
      static void calculate_anderson_coeffs(Scalar** previous_vectors, Scalar* anderson_coeffs, int num_last_vectors_used, int ndof);
      
      bool verbose_output_linear_solver;

      /// Matrix.
      SparseMatrix<Scalar>* matrix;

      /// Right-hand side.
      Vector<Scalar>* rhs;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* linear_solver;

      /// This instance owns its DP.
      const bool own_dp;

      double tol;
      int max_iter;
      int num_last_vectors_used;
      bool anderson_is_on;
      double anderson_beta;
    };
  }
}
#endif
