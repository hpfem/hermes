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

#include "../h3dconfig.h"
#include "amesos.h"
#include "../linear_problem.h"
#include <common/callstack.h>
#include <common/timer.h>

#ifdef HAVE_AMESOS
#include <Amesos_ConfigDefs.h>
#endif

#ifdef HAVE_AMESOS
Amesos AmesosSolver::factory;
#endif

// Amesos solver ///////////////////////////////////////////////////////////////////////////////////

AmesosSolver::AmesosSolver(const char *solver_type, EpetraMatrix *m, EpetraVector *rhs)
	: LinearSolver(), m(m), rhs(rhs)
{
	_F_
#ifdef HAVE_AMESOS
	solver = factory.Create(solver_type, problem);
	assert(solver != NULL);
#else
	warning("hermes3d was not built with AMESOS support.");
	exit(128);
#endif
}

AmesosSolver::AmesosSolver(const char *solver_type, LinearProblem *lp)
	: LinearSolver(lp)
{
	_F_
#ifdef HAVE_AMESOS
	solver = factory.Create(solver_type, problem);
	assert(solver != NULL);
	m = new EpetraMatrix;
	rhs = new EpetraVector;
#else
	warning("hermes3d was not built with AMESOS support.");
	exit(128);
#endif
}

AmesosSolver::~AmesosSolver()
{
	_F_
#ifdef HAVE_AMESOS
	if (lp != NULL) {
		delete m;
		delete rhs;
	}
	delete solver;
#endif
}

bool AmesosSolver::is_available(const char *name)
{
	_F_
#ifdef HAVE_AMESOS
	return factory.Query(name);
#else
	return false;
#endif
}

void AmesosSolver::set_use_transpose(bool use_transpose)
{
	_F_
#ifdef HAVE_AMESOS
	solver->SetUseTranspose(use_transpose);
#endif
}

bool AmesosSolver::use_transpose()
{
	_F_
#ifdef HAVE_AMESOS
	return solver->UseTranspose();
#else
	return false;
#endif
}

bool AmesosSolver::solve()
{
	_F_
#ifdef HAVE_AMESOS
	assert(m != NULL);
	assert(rhs != NULL);

	if (lp != NULL)
		lp->assemble(m, rhs);
	assert(m->size == rhs->size);

	Timer tmr;
	tmr.start();

	Epetra_Vector x(*rhs->std_map);

	problem.SetOperator(m->mat);
	problem.SetRHS(rhs->vec);
	problem.SetLHS(&x);

	if ((error = solver->SymbolicFactorization()) != 0) return false;
	if ((error = solver->NumericFactorization()) != 0) return false;
	if ((error = solver->Solve()) != 0) return false;

	tmr.stop();
	time = tmr.get_seconds();

	delete [] sln;
	sln = new scalar[m->size]; MEM_CHECK(sln);
	// copy the solution into sln vector
	memset(sln, 0, m->size * sizeof(scalar));
	for (int i = 0; i < m->size; i++) sln[i] = x[i];

	return true;
#else
	return false;
#endif
}
