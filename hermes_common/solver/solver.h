// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "precond.h"
                        
#ifdef HAVE_TEUCHOS
  #include <Teuchos_RefCountPtr.hpp>
#endif

/// @defgroup solvers Solvers
///
/// TODO: description
///
/*@{*/ // Beginning of documentation group Solvers.

/// Options for matrix factorization reuse.
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
///   -\c Pardiso - performs reordering and scaling in one step, during the symbolic analysis. Both 
///                 options \c HERMES_REUSE_MATRIX_REORDERING_AND_SCALING and 
///                 \c HERMES_REUSE_MATRIX_REORDERING have the same meaning, i.e. results of the 
///                 symbolic phase will be reused and L, U will be computed anew.
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
///
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

class DiscreteProblem;

/// Abstract class for defining solver interface.
///
/// TODO: Adjust interface to support faster update of matrix and rhs
///
class Solver {
public:
  Solver() { sln = NULL; time = -1.0; }
  virtual ~Solver() { if (sln != NULL) delete [] sln; }

  virtual bool solve() = 0;
  scalar *get_solution() { return sln; }

  int get_error() { return error; }
  double get_time() { return time; }
  
  virtual void set_factorization_scheme(FactorizationScheme reuse_scheme) { };
  virtual void set_factorization_scheme() {
    set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY); 
  }

protected:
  scalar *sln;
  int error;
  double time;  ///< time spent on solving (in secs)
};


/// Abstract class for defining interface for linear solvers.
///
class LinearSolver : public Solver 
{
  public:
    LinearSolver(unsigned int factorization_scheme = HERMES_FACTORIZE_FROM_SCRATCH) 
      : Solver(), factorization_scheme(factorization_scheme) {};
    
  protected:
    virtual void set_factorization_scheme(FactorizationScheme reuse_scheme) { 
      factorization_scheme = reuse_scheme;
    }
        
    unsigned int factorization_scheme;
};

/// Abstract class for defining interface for nonlinear solvers.
///
class NonlinearSolver : public Solver {
  public:
    NonlinearSolver() : Solver() { dp = NULL; }
    NonlinearSolver(DiscreteProblem* dp) : Solver() { this->dp = dp; }
    
  protected:
    DiscreteProblem* dp;        // FE problem being solved (not NULL in case of using
    // NonlinearProblem(DiscreteProblem *) ctor
};

/// Abstract class for defining interface for iterative solvers.
///
class IterSolver : public Solver
{
  public:
    IterSolver() : Solver(), max_iters(10000), tolerance(1e-8), precond_yes(false) {};
    
    virtual int get_num_iters() = 0;
    virtual double get_residual() = 0;
    
    /// Set the convergence tolerance
    /// @param[in] tol - the tolerance to set
    void set_tolerance(double tol) { this->tolerance = tol; }
    /// Set maximum number of iterations to perform
    /// @param[in] iters - number of iterations
    void set_max_iters(int iters) { this->max_iters = iters; }
    
    virtual void set_precond(const char *name) = 0;
    #ifdef HAVE_TEUCHOS
      virtual void set_precond(Teuchos::RCP<Precond> &pc) = 0;
    #else
      virtual void set_precond(Precond *pc) = 0;
    #endif            
      
  protected:    
    int max_iters;          ///< Maximum number of iterations.
    double tolerance;       ///< Convergence tolerance.
    bool precond_yes;
};

/*@}*/ // End of documentation group Solvers.

#endif
