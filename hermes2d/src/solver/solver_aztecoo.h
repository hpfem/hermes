// This file is part of Hermes2D
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

#ifndef __H2D_SOLVER_AZTECOO_H_
#define __H2D_SOLVER_AZTECOO_H_

#include "solver_epetra.h"
#include "itersolver.h"
#include "precond_ifpack.h"

#ifdef HAVE_AZTECOO
#include <AztecOO.h>
#endif

/// Encapsulation of AztecOO linear solver
///
/// @ingroup solvers
class AztecOOSolver : public IterSolver
{
public:
	AztecOOSolver(EpetraMatrix &m, EpetraVector &rhs);
	virtual ~AztecOOSolver();

	virtual bool solve();

	int get_num_iters();

	/// Set the type of the solver
	/// @param[in] solver - name of the solver [ gmres | cg | cgs | tfqmr | bicgstab ]
	void set_solver(const char *solver);
	/// Set the convergence tolerance
	/// @param[in] tol - the tolerance to set
	void set_tolerance(double tol) { this->tolerance = tol; }
	/// Set maximum number of iterations to perform
	/// @param[in] iters - number of iterations
	void set_max_iters(int iters) { this->max_iters = iters; }

	/// Set Aztec internal preconditioner
	/// @param[in] name - name of the preconditioner [ none | jacobi | neumann | least-squares ]
	void set_precond(const char *name);
	/// Set preconditioner from IFPACK
	/// @param[in] pc - IFPACK preconditioner
	void set_precond(Precond *pc) { this->pc = pc; }

protected:
#ifdef WITH_TRILINOS
	AztecOO aztec;					/// instance of aztec solver
#endif
	EpetraMatrix &m;
	EpetraVector &rhs;
	Precond *pc;

	int max_iters;					/// maximum number of iterations
	double tolerance;				/// convergence tolerance
};

#endif
