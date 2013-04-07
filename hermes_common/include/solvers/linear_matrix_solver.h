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
    ///&nbsp;              previously created preconditioner, \c HERMES_REUSE_FACTORIZATION_COMPLETELY
    ///&nbsp;              indicates that the preconditioner may be reused completely for future solves.
    ///
    /// <b>Typical scenario:</b>
    /// When \c rhsonly was set to \c true for the assembly phase,
    /// \c HERMES_REUSE_FACTORIZATION_COMPLETELY should be set for the following solution phase.
    enum FactorizationScheme
    {
      HERMES_FACTORIZE_FROM_SCRATCH,              ///< Perform new factorization, don't reuse
      ///< existing factorization data.
      HERMES_REUSE_MATRIX_REORDERING,             ///< Factorize matrix with the same sparsity
      ///< pattern as the one already factorized.
      HERMES_REUSE_MATRIX_REORDERING_AND_SCALING, ///< Factorize matrix with the same sparsity
      ///< pattern and similar numerical values
      ///< as the one already factorized.
      HERMES_REUSE_FACTORIZATION_COMPLETELY       ///< Completely reuse the already performed
      ///< factorization.
    };

    /// \brief Abstract class for defining solver interface.
    ///
    ///\todo Adjust interface to support faster update of matrix and rhs
    ///
    /// Typical usage is through the function create_linear_solver(Matrix<Scalar>* matrix, Vector<Scalar>* rhs).
    template <typename Scalar>
    class LinearMatrixSolver : public Hermes::Mixins::Loggable, public Hermes::Mixins::TimeMeasurable
    {
    public:
      LinearMatrixSolver();

      virtual ~LinearMatrixSolver();

      /// Solve.
      /// @return true on succes
      virtual bool solve() = 0;

      /// Get solution vector.
      /// @return solution vector ( #sln )
      Scalar *get_sln_vector();

      /// @return #error
      int get_error();
      /// Get time spent on solving.
      /// @return time spent on solving in secs ( #time )
      double get_time();
      /// Get size of matrix
      virtual int get_matrix_size() = 0;

      /// Set factorization scheme.
      /// @param[in] reuse_scheme factoriztion scheme to set
      virtual void set_factorization_scheme(FactorizationScheme reuse_scheme) { };

      /// Set factorization scheme to default.
      virtual void set_factorization_scheme();

    protected:
      /// Solution vector.
      Scalar *sln;
      /// \todo document (not sure what it do)
      int error;
      double time;  ///< Time spent on solving (in secs).
    };

    /// \brief Base class for defining interface for direct linear solvers.
    /// Internal, though utilizable for defining interfaces to other algebraic packages.
    template <typename Scalar>
    class DirectSolver : public LinearMatrixSolver<Scalar>
    {
    public:
      DirectSolver(unsigned int factorization_scheme = HERMES_FACTORIZE_FROM_SCRATCH)
        : LinearMatrixSolver<Scalar>(), factorization_scheme(factorization_scheme) {};

    protected:
      virtual void set_factorization_scheme(FactorizationScheme reuse_scheme);

      unsigned int factorization_scheme;
    };

    /// \brief  Abstract class for defining interface for iterative solvers.
    /// Internal, though utilizable for defining interfaces to other algebraic packages.
    template <typename Scalar>
    class IterSolver : public LinearMatrixSolver<Scalar>
    {
    public:
      IterSolver() : LinearMatrixSolver<Scalar>(), max_iters(10000), tolerance(1e-8), precond_yes(false) {};

      virtual int get_num_iters() = 0;
      virtual double get_residual() = 0;

      /// Set the convergence tolerance.
      /// @param[in] tol - the tolerance to set
      void set_tolerance(double tol);

      /// Set maximum number of iterations to perform.
      /// @param[in] iters - number of iterations
      void set_max_iters(int iters);

      virtual void set_precond(const char *name) = 0;

      virtual void set_precond(Precond<Scalar> *pc) = 0;

    protected:
      int max_iters;          ///< Maximum number of iterations.
      double tolerance;       ///< Convergence tolerance.
      bool precond_yes;
    };

    /// \brief Function returning a solver according to the users's choice.
    /// @param[in] matrix matrix
    /// @param[in] rhs right hand side vector
    /// @return created linear solver
    template<typename Scalar>
    HERMES_API LinearMatrixSolver<Scalar>*
      create_linear_solver(Matrix<Scalar>* matrix, Vector<Scalar>* rhs);
  }
}
/*@}*/ // End of documentation group Solvers.
#endif