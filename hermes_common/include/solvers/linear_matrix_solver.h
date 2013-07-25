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
/*! \file linear_solver.h
\brief General linear/iterative solver functionality.
*/
#ifndef __HERMES_COMMON_SOLVER_H_
#define __HERMES_COMMON_SOLVER_H_

#include "precond.h"
#include "exceptions.h"
#include "mixins.h"

using namespace Hermes::Algebra;
using namespace Hermes::Solvers;

/// \brief General namespace for the Hermes library.
namespace Hermes
{
  /// \brief Namespace for linear / nonlinear / iterative solvers.
  namespace Solvers
  {
    using namespace Hermes::Preconditioners;
    /// \ brief Options for matrix factorization reuse.
    ///
    /// Reusing the information computed during previous solution of a similar problem
    /// significantly improves efficiency of the solver.
    ///
    /// <b>Usage:</b>
    /// Each solver which allows factorization reuse should perform complete factorization
    /// from scratch for the first time it is invoked, keep the precomputed structures
    /// according to the current factorization reuse stratregy and use them for next
    /// factorization.
    ///
    /// <b>Enabled solvers:</b>
    ///&nbsp;-\c SuperLU - performs reordering, scaling and factorization separately. When the
    ///&nbsp;              multithreaded version is used, scaling is performed during the factorization
    ///&nbsp;              phase (if neccessary) and thus \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING
    ///&nbsp;              and \c HERMES_REUSE_MATRIX_REORDERING have the same effect.
    ///&nbsp;-\c UMFPack - like the MT version of SuperLU, performs scaling and factorization in one step.
    ///&nbsp;              \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING has thus the same effect as
    ///&nbsp;              \c HERMES_REUSE_MATRIX_REORDERING (saves the preceding symbolic analysis step).
    ///&nbsp;-\c MUMPS   - If \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING is set, scaling is performed
    ///&nbsp;              during analysis and only factorization is repeated during each solve.
    ///&nbsp;              If \c HERMES_REUSE_MATRIX_REORDERING is set, scaling is performed during
    ///&nbsp;              the factorization phase. This may be less efficient, but more reliable,
    ///&nbsp;              especially for highly non-symmetric matrices.
    ///&nbsp;-\c AMESOS  - Behaves like UMFPack.
    ///&nbsp;-\c PETSc   - Factorization reuse applies to the construction of PETSc preconditioner.
    ///&nbsp;              Both \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING and
    ///&nbsp;              \c HERMES_REUSE_MATRIX_REORDERING allow to reuse the non-zero pattern of the
    ///&nbsp;              previously created preconditioner, \c HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY
    ///&nbsp;              indicates that the preconditioner may be reused completely for future solves.
    ///
    /// <b>Typical scenario:</b>
    /// When \c rhsonly was set to \c true for the assembly phase,
    /// \c HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY should be set for the following solution phase.
    enum MatrixStructureReuseScheme
    {
      HERMES_CREATE_STRUCTURE_FROM_SCRATCH,              ///< Perform new factorization (operator creation), don't reuse
      ///< existing factorization data.
      HERMES_REUSE_MATRIX_REORDERING,             ///< Factorize matrix (create operatoer) with the same sparsity
      ///< pattern as the one already factorized.
      HERMES_REUSE_MATRIX_REORDERING_AND_SCALING, ///< Factorize matrix with the same sparsity
      ///< pattern and similar numerical values
      ///< as the one already factorized.
      HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY       ///< Completely reuse the already performed
      ///< factorization / operator.
    };

    /// \brief Abstract class for defining solver interface.
    ///
    ///\todo Adjust interface to support faster update of matrix and rhs
    ///
    /// Typical usage is through the function create_linear_solver(Matrix<Scalar>* matrix, Vector<Scalar>* rhs, bool use_direct_solver).
    template <typename Scalar>
    class LinearMatrixSolver : public Hermes::Mixins::Loggable, public Hermes::Mixins::TimeMeasurable
    {
    public:
      LinearMatrixSolver(MatrixStructureReuseScheme reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

      virtual ~LinearMatrixSolver();

      /// Solve.
      /// @return true on succes
      virtual void solve() = 0;

      /// Solve.
      /// @return true on succes
      /// \param[in] initial guess.
      virtual void solve(Scalar* initial_guess) = 0;

      /// Get solution vector.
      /// @return solution vector ( #sln )
      Scalar *get_sln_vector();

      /// Get time spent on solving.
      /// @return time spent on solving in secs ( #time )
      double get_time();

      /// Get size of matrix
      virtual int get_matrix_size() = 0;

      /// Get factorization scheme.
      virtual MatrixStructureReuseScheme get_used_reuse_scheme() const;

      /// Set factorization scheme.
      /// @param[in] reuse_scheme factoriztion scheme to set
      virtual void set_reuse_scheme(MatrixStructureReuseScheme reuse_scheme);

      /// Set factorization scheme to default.
      virtual void set_reuse_scheme();
      
      /// Set matrix ordering in the case of a system of PDEs.
      virtual void use_node_wise_ordering(unsigned int num_pdes);
      virtual void use_equations_wise_ordering();

    protected:
      /// Factorization scheme
      MatrixStructureReuseScheme reuse_scheme;

      /// Solution vector.
      Scalar *sln;

      ///< Time spent on solving (in secs).
      double time;
      
      /// Number of equations in a system of PDEs.
      unsigned int n_eq;
      
      bool node_wise_ordering;
    };

    /// \brief Base class for defining interface for direct linear solvers.
    /// Internal, though utilizable for defining interfaces to other algebraic packages.
    template <typename Scalar>
    class DirectSolver : public LinearMatrixSolver<Scalar>
    {
    public:
      DirectSolver(MatrixStructureReuseScheme reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
      virtual void solve() = 0;
      virtual void solve(Scalar* initial_guess);
    };

    /// \brief Abstract middle-class for solvers that work in a loop of a kind (iterative, multigrid, ...)
    template <typename Scalar>
    class LoopSolver : public LinearMatrixSolver<Scalar>
    {
    public:
      LoopSolver(MatrixStructureReuseScheme reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

      /// Various tolerances.
      /// Not necessarily supported by all iterative solvers used.
      enum ToleranceType
      {
        AbsoluteTolerance = 0,
        RelativeTolerance = 1,
        DivergenceTolerance = 2
      };

      /// Get the number of iterations performed.
      virtual int get_num_iters() = 0;
      
      /// Get the final residual.
      virtual double get_residual() = 0;

      /// Set the convergence tolerance.
      /// @param[in] tol - the tolerance to set
      virtual void set_tolerance(double tol);

      /// Set the convergence tolerance.
      /// @param[in] tolerance - the tolerance to set
      /// @param[in] toleranceType - the tolerance to set
      virtual void set_tolerance(double tolerance, ToleranceType toleranceType);

      /// Set maximum number of iterations to perform.
      /// @param[in] iters - number of iterations
      virtual void set_max_iters(int iters);

    protected:
      /// Maximum number of iterations.
      int max_iters;
      /// Convergence tolerance.
      double tolerance;
      /// Convergence tolerance type.
      /// See the enum.
      ToleranceType toleranceType;
    };
    

    /// \brief  Abstract class for defining interface for iterative solvers.
    /// Internal, though utilizable for defining interfaces to other algebraic packages.
    template <typename Scalar>
    class IterSolver : public LoopSolver<Scalar>
    {
    public:
      IterSolver(MatrixStructureReuseScheme reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH);

      /// Set preconditioner.
      virtual void set_precond(Precond<Scalar> *pc) = 0;

    protected:
      /// Whether the solver is preconditioned.
      bool precond_yes;
    };
    
    /// \brief  Abstract class for defining interface for Algebraic Multigrid solvers.
    /// Internal, though utilizable for defining interfaces to other algebraic packages.
    template <typename Scalar>
    class AMGSolver : public LoopSolver<Scalar>
    {
    public:
      AMGSolver(MatrixStructureReuseScheme reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    };

    /// \brief Function returning a solver according to the users's choice.
    /// @param[in] matrix matrix
    /// @param[in] rhs right hand side vector
    /// @return created linear solver
    template<typename Scalar>
    HERMES_API LinearMatrixSolver<Scalar>*
      create_linear_solver(Matrix<Scalar>* matrix, Vector<Scalar>* rhs, bool use_direct_solver = false);
  }
}
/*@}*/ // End of documentation group Solvers.
#endif