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

#ifndef _AZTECOO_SOLVER_H_
#define _AZTECOO_SOLVER_H_

#include "epetra.h"
#include "../solver.h"
#include "../precond/ifpack.h"

#ifdef HAVE_AZTECOO
#include <AztecOO.h>
#endif

/// Encapsulation of AztecOO linear solver
///
/// @ingroup solvers
class AztecOOSolver : public LinearSolver {
public:
	AztecOOSolver(EpetraMatrix *m, EpetraVector *rhs);
	AztecOOSolver(LinearProblem *lp);
	virtual ~AztecOOSolver();

	virtual bool solve();

	int get_num_iters();
	double get_residual();

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

	/// Option setting function
	void set_option(int option, int value);

	/// Parameter setting function
	void set_param(int param, double value);

protected:
#ifdef WITH_TRILINOS
	AztecOO aztec;					/// instance of aztec solver
#endif
	EpetraMatrix *m;
	EpetraVector *rhs;
	Precond *pc;

	int max_iters;					/// maximum number of iterations
	double tolerance;				/// convergence tolerance
};

#endif
