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
#include "dpinterface.h"

#ifdef HAVE_TEUCHOS
#include <Teuchos_RefCountPtr.hpp>
#endif

using namespace Hermes::Algebra;

/// @defgroup solvers Solvers
///
///\todo description
///
/*@{*/ // Beginning of documentation group Solvers.

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
    ///   -\c SuperLU - performs reordering, scaling and factorization separately. When the 
    ///                 multithreaded version is used, scaling is performed during the factorization
    ///                 phase (if neccessary) and thus \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING
    ///                 and \c HERMES_REUSE_MATRIX_REORDERING have the same effect.
    ///   -\c UMFPack - like the MT version of SuperLU, performs scaling and factorization in one step.
    ///                 \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING has thus the same effect as
    ///                 \c HERMES_REUSE_MATRIX_REORDERING (saves the preceding symbolic analysis step).
    ///   -\c MUMPS   - If \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING is set, scaling is performed
    ///                 during analysis and only factorization is repeated during each solve.
    ///                 If \c HERMES_REUSE_MATRIX_REORDERING is set, scaling is performed during 
    ///                 the factorization phase. This may be less efficient, but more reliable,
    ///                 especially for highly non-symmetric matrices.
    ///   -\c AMESOS  - Behaves like UMFPack.
    ///   -\c PETSc   - Factorization reuse applies to the construction of PETSc preconditioner.
    ///                 Both \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING and 
    ///                 \c HERMES_REUSE_MATRIX_REORDERING allow to reuse the non-zero pattern of the
    ///                 previously created preconditioner, \c HERMES_REUSE_FACTORIZATION_COMPLETELY
    ///                 indicates that the preconditioner may be reused completely for future solves.
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
    template <typename Scalar>
    class LinearSolver {
    public:
      LinearSolver() { sln = NULL; time = -1.0; }

      virtual ~LinearSolver() { if (sln != NULL) delete [] sln; }

      /// Solve.
      /// @return true on succes
      virtual bool solve() = 0;

      /// Get solution vector.
      /// @return solution vector ( #sln )
      Scalar *get_sln_vector() { return sln; }

      /// @return #error
      int get_error() { return error; }
      /// Get time spent on solving.
      /// @return time spent on solving in secs ( #time )
      double get_time() { return time; }

      /// Set factorization scheme.
      /// @param[in] reuse_scheme factoriztion scheme to set
      virtual void set_factorization_scheme(FactorizationScheme reuse_scheme) { };

      /// Set factorization scheme to default.
      virtual void set_factorization_scheme() 
      {
        set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY); 
      }

    protected:
      /// Solution vector.
      Scalar *sln;
      /// \todo document (not sure what it do)
      int error;
      double time;  ///< Time spent on solving (in secs).
    };


    /// \brief Base class for defining interface for direct linear solvers.
    ///
    template <typename Scalar>
    class DirectSolver : virtual public LinearSolver<Scalar>
    {
    public:
      DirectSolver(unsigned int factorization_scheme = HERMES_FACTORIZE_FROM_SCRATCH) 
        : LinearSolver<Scalar>(), factorization_scheme(factorization_scheme) {};

    protected:
      virtual void set_factorization_scheme(FactorizationScheme reuse_scheme) { 
        factorization_scheme = reuse_scheme;
      }

      unsigned int factorization_scheme;
    };

    /// \brief  Abstract class for defining interface for iterative solvers.
    ///
    template <typename Scalar>
    class IterSolver : public virtual LinearSolver<Scalar>
    {
    public:
      IterSolver() : LinearSolver<Scalar>(), max_iters(10000), tolerance(1e-8), precond_yes(false) {};

      virtual int get_num_iters() = 0;
      virtual double get_residual() = 0;

      /// Set the convergence tolerance.
      /// @param[in] tol - the tolerance to set
      void set_tolerance(double tol)
      {
        this->tolerance = tol;
      }

      /// Set maximum number of iterations to perform.
      /// @param[in] iters - number of iterations
      void set_max_iters(int iters)
      {
        this->max_iters = iters; 
      }

      virtual void set_precond(const char *name) = 0;

#ifdef HAVE_TEUCHOS
      virtual void set_precond(Teuchos::RCP<Precond<Scalar> > &pc) = 0;
#else
      virtual void set_precond(Precond<Scalar> *pc) = 0;
#endif            

    protected:    
      int max_iters;          ///< Maximum number of iterations.
      double tolerance;       ///< Convergence tolerance.
      bool precond_yes;
    };

    /// \brief Function returning a solver according to the users's choice.
    /// @param[in] matrix_solver_type the choice of solver, an element of enum Hermes::MatrixSolverType.
    /// @param[in] matrix matrix
    /// @param[in] rhs right hand side vector
    /// @return created linear solver
    template<typename Scalar> 
    HERMES_API LinearSolver<Scalar>*  
      create_linear_solver(Hermes::MatrixSolverType matrix_solver_type, Matrix<Scalar>* matrix, Vector<Scalar>* rhs);

  }
}
/*@}*/ // End of documentation group Solvers.
#endif
