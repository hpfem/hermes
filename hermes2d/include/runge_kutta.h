// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_RUNGE_KUTTA_H
#define __H2D_RUNGE_KUTTA_H

#include "global.h"
#include "weakform/weakform.h"
#include "function/filter.h"
#include "exceptions.h"
#include "mixins2d.h"
namespace Hermes
{
  namespace Hermes2D
  {
    // TODO LIST:
    //
    // (1) With explicit and diagonally implicit methods, the matrix is treated
    //     in the same way as with fully implicit ones. To make this more
    //     efficient, with explicit and diagonally implicit methods one should
    //     first only solve for the upper left block, then eliminate all blocks
    //     under it, then solve for block at position 22, eliminate all blocks
    //     under it, etc. Currently this is not done and everything is left to
    //     the matrix solver.
    //
    // (2) In example 03-timedep-adapt-space-and-time with implicit Euler
    //     method, Newton's method takes much longer than in 01-timedep-adapt-space-only
    //     (that also uses implicit Euler method). This means that the initial guess for
    //     the K_vector should be improved (currently it is zero).
    //
    // (3) At the end of rk_time_step_newton(), the previous time level solution is
    //     projected onto the space of the new time-level solution so that
    //     it can be added to the stages. This projection is slow so we should
    //     find a way to do this differently. In any case, the projection
    //     is not necessary when no adaptivity in space takes place and the
    //     two spaces are the same (but it is done anyway).
    //
    // (4) We do not take advantage of the fact that all blocks in the
    //     Jacobian matrix have the same structure. Thus it is enough to
    //     assemble the matrix M (one block) and copy the sparsity structure
    //     into all remaining nonzero blocks (and diagonal blocks). Right
    //     now, the sparsity structure is created expensively in each block
    //     again.
    //
    // (5) If space does not change, the sparsity does not change. Right now
    //     we discard everything at the end of every time step, we should not
    //     do it.
    //
    // (6) If the problem does not depend explicitly on time, then all the blocks
    //     in the Jacobian matrix of the stationary residual are the same up
    //     to a multiplicative constant. Thus they do not have to be aassembled
    //     from scratch.
    //
    //
    // (7) In practice, Butcher's tables are being transformed to the
    //     Jordan canonical form (I think) for better performance. This
    //     can be found, I think, in newer Butcher's papers or presentation
    //     (he has them online), and possibly in his book.


    /// @defgroup userSolvingAPI Top-level solving

    /// @ingroup userSolvingAPI
    template<typename Scalar>
    class HERMES_API RungeKutta : public Hermes::Mixins::Loggable, public Hermes::Mixins::IntegrableWithGlobalOrder, public Hermes::Mixins::SettableComputationTime, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>
    {
    public:
      /// Constructor.
      /// Parameter start_from_zero_K_vector: if set to true, the last K_vector will NOT be used
      /// as an initial guess for the Newton's method, instead zero vector will be used.
      RungeKutta(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces, ButcherTable* bt);

      /// Constructor for one equation.
      RungeKutta(const WeakForm<Scalar>* wf, const Space<Scalar>* space, ButcherTable* bt);

      /// Projections will be local (projection-based).
      void use_local_projections();

      void set_start_from_zero_K_vector();
      void set_residual_as_solutions();
      void set_block_diagonal_jacobian();

      /// Destructor.
      ~RungeKutta();

      // Perform one explicit or implicit time step using the Runge-Kutta method
      // corresponding to a given Butcher's table. If err_vec != NULL then it will be
      // filled with an error vector calculated using the second B-row of the Butcher's
      // table (the second B-row B2 must be nonzero in that case). The negative default
      // values for newton_tol and newton_max_iter are for linear problems.
      // Many improvements are needed, a todo list is presented at the beginning of
      // the corresponding .cpp file.
      // freeze_jacobian... if true then the Jacobian is not recalculated in each
      //                    iteration of the Newton's method.
      // block_diagonal_jacobian... if true then the tensor product block Jacobian is
      //                            reduced to just the diagonal blocks.
      void rk_time_step_newton(Hermes::vector<Solution<Scalar>*> slns_time_prev, Hermes::vector<Solution<Scalar>*> slns_time_new, Hermes::vector<Solution<Scalar>*> error_fns);
      void rk_time_step_newton(Solution<Scalar>* slns_time_prev, Solution<Scalar>* slns_time_new, Solution<Scalar>* error_fn);

      // This is a wrapper for the previous function if error_fn is not provided
      // (adaptive time stepping is not wanted).
      void rk_time_step_newton(Hermes::vector<Solution<Scalar>*> slns_time_prev, Hermes::vector<Solution<Scalar>*> slns_time_new);
      void rk_time_step_newton(Solution<Scalar>* sln_time_prev, Solution<Scalar>* sln_time_new);

      void set_freeze_jacobian();
      void set_newton_tol(double newton_tol);
      void set_newton_max_iter(int newton_max_iter);
      void set_newton_damping_coeff(double newton_damping_coeff);
      void set_newton_max_allowed_residual_norm(double newton_max_allowed_residual_norm);

      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;

      /**
       \fn  void RungeKutta::set_filters_to_reinit(Hermes::vector<Filter<Scalar>*> filters_to_reinit);

       \brief Sets the filters to reinitialize.

       \author  Lk
       \date  10/29/2011

       \param[in]  filters_to_reinit the filters to reinitialize.
       */

      void set_filters_to_reinit(Hermes::vector<Filter<Scalar>*> filters_to_reinit);

    protected:
      /// Initialization of DiscreteProblems.
      void init();

      /// Takes a matrix M of size ndof times ndof, extends it (formally) to
      /// a num_stages*ndof times num_stages*ndof matrix that has M in diagonal blocks and
      /// zero everywhere else, and multiplies the new matrix with the vector stage_coeff_vec
      /// which has length num_stages*ndof. The result is saved in vector_left which also
      /// has length num_stages*ndof.
      /// TODO: enable this for other types of matrices.
      void multiply_as_diagonal_block_matrix(SparseMatrix<Scalar>* matrix_left, int num_stages,
        Scalar* stage_coeff_vec, Scalar* vector_left);

      /// Creates an augmented weak formulation for the multi-stage Runge-Kutta problem.
      /// The original discretized equation is M\dot{Y} = F(t, Y) where M is the mass
      /// matrix, Y the coefficient vector, and F the (nonlinear) stationary residual.
      /// Below, "stage_wf_left" and "stage_wf_right" refer to the left-hand side M\dot{Y}
      /// and right-hand side F(t, Y) of the above equation, respectively.
      void create_stage_wf(unsigned int size, bool block_diagonal_jacobian);
      
      /// Updates the augmented weak formulation.
      void update_stage_wf(Hermes::vector<Solution<Scalar>*> slns_time_prev);

      // Prepare u_ext_vec.
      void prepare_u_ext_vec();

      /// Matrix for the time derivative part of the equation (left-hand side).
      SparseMatrix<Scalar>* matrix_left;

      /// Matrix and vector for the rest (right-hand side).
      SparseMatrix<Scalar>* matrix_right;
      Vector<Scalar>* vector_right;

      /// Matrix solver.
      LinearMatrixSolver<Scalar>* solver;

      /// Weak formulation.
      const WeakForm<Scalar>* wf;

      /// Space instances for all equations in the system.
      Hermes::vector<const Space<Scalar>*> spaces;
      Hermes::vector<Space<Scalar>*> spaces_mutable;

      /// ButcherTable.
      ButcherTable* bt;

      /// Number of stages.
      unsigned int num_stages;

      /// Multistage weak formulation.
      // For the main part equation (written on the right),
      /// size num_stages*ndof times num_stages*ndof.
      WeakForm<Scalar> stage_wf_right;
      DiscreteProblem<Scalar>* stage_dp_right;

      /// For the matrix M (size ndof times ndof).
      WeakForm<Scalar> stage_wf_left;
      DiscreteProblem<Scalar>* stage_dp_left;

      bool start_from_zero_K_vector;
      bool block_diagonal_jacobian;
      bool residual_as_vector;

      /// Number of previous calls to rk_time_step_newton().
      unsigned int iteration;

      bool do_global_projections;
      bool freeze_jacobian;
      double newton_tol;
      int newton_max_iter;
      double newton_damping_coeff;
      double newton_max_allowed_residual_norm;
      
      Hermes::vector<Solution<Scalar>*> residuals_vector;

      /// Vector K_vector of length num_stages * ndof. will represent
      /// the 'K_i' vectors in the usual R-K notation.
      Scalar* K_vector;

      /// Vector u_ext_vec will represent h \sum_{j = 1}^s a_{ij} K_i.
      Scalar* u_ext_vec;

      /// Vector for the left part of the residual.
      Scalar* vector_left;
      
      ///< The filters to reinitialize in every Newton's loop
      Hermes::vector<Filter<Scalar>*> filters_to_reinit;
    };
  }
}
#endif