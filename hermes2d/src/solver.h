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

#ifndef __H2D_SOLVER_H
#define __H2D_SOLVER_H

class DiscreteProblem;
class LinearProblem;

/// \brief Abstract interface to sparse linear solvers.
///
///  Solver is an abstract class defining the interface to all linear solvers
///  used by Hermes2D. A concrete derived class (UmfpackSolver, PardisoSolver...)
///  is instantiated by the user and passed to the DiscreteProblem class to solve the
///  discrete problem. The user never directly calls any of the methods of this
///  class.
///
///  The linear solver can be direct or iterative, although the analyze() and
///  factorize() methods are clearly designed for direct solvers. These
///  methods are optional and do not have to be implemented by iterative solvers.
///
///  Another feature geared towards direct solvers is the possibility to store
///  auxiliary data (usually factorization data) between the calls to analyze(),
///  factorize() and solve() in an arbitrary data structure called "context".
///  This way the factorization can be reused to solve for different right hand
///  sides and/or different, but structurally identical matrices.
///
class Solver
{
protected:

  /// Must return true if the solvers expects compressed row (CSR) format.
  /// Otherwise DiscreteProblem assumes the compressed column (CSC) format.
  virtual bool is_row_oriented() = 0;

  /// Must return true if the solver is capable of solving structurally
  /// symmetric matrices, ie., when only the upper half of the matrix is
  /// stored in the CSR or CSC arrays.
  virtual bool handles_symmetry() = 0;


  /// Creates a new data block containing (optional) factorization data.
  /// This method is called by DiscreteProblem on its creation.
  virtual void* new_context(bool sym) { return NULL; }

  /// Frees the data block created by new_context(). Called by DiscreteProblem on destruction.
  virtual void free_context(void *ctx) {}


  /// After the sparse structure of the matrix is calculated, DiscreteProblem calls
  /// this function to give the solver a chance to analyze the matrix and store
  /// the results for reuse in the execution context. If the solver does not
  /// support the reuse of structural analysis, this method does not have to be
  /// implemented.  \return true on success, false otherwise.
  virtual bool analyze(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym) { return false; }

  /// Called by DiscreteProblem after the stiffness matrix has been assembled.
  /// Direct solvers should implement this function and store the result
  /// of the factorization in the execution context, so that it can be used
  /// many times by solve() for different right hand sides.
  /// \return true on success, false otherwise.
  virtual bool factorize(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym) { return false; }

  /// Called by DiscreteProblem when the user requests the solution of the linear system.
  /// Direct solvers will want to use the matrix factorization stored in "ctx".
  /// Iterative solvers will probably solve the system from scratch in this call,
  /// but can expect "vec" to contain the initial guess for the solution vector.
  /// \param ctx execution context
  /// \param n [in] number of rows and columns of the matrix
  /// \param Ap, Ai, Ax [in] the matrix in CSR or CSC format
  /// \param sym [in] true if the matrix is symmetric and only the upper half is stored.
  /// \param RHS [in] right hand side vector
  /// \param vec [in/out] solution vector
  /// \return true on success, false otherwise.
  virtual bool solve(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym,
                     scalar* RHS, scalar* vec) = 0;

  /// Must free all auxiliary data created in the calls to analyze() and/or factorize().
  virtual void free_data(void* ctx) {};


  friend class DiscreteProblem;
  friend class LinearProblem;

};


#endif
