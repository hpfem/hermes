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

#include "config.h"
#include "solver_aztecoo.h"
#ifdef HAVE_KOMPLEX
#include <Komplex_LinearProblem.h>
#endif

#define AZTECOO_NOT_COMPILED "hermes2d was not built with AztecOO support."

// AztecOO solver //////////////////////////////////////////////////////////////////////////////////

AztecOOSolver::AztecOOSolver(EpetraMatrix &m, EpetraVector &rhs)
	: IterSolver(), m(m), rhs(rhs)
{
#ifdef HAVE_AZTECOO
	// set default values
	max_iters = 10000;
	tolerance = 10e-8;
	pc = NULL;
#else
	error(AZTECOO_NOT_COMPILED);
#endif
}

AztecOOSolver::~AztecOOSolver()
{
#ifdef HAVE_AZTECOO
#endif
}

void AztecOOSolver::set_solver(const char *name)
{
#ifdef HAVE_AZTECOO
	int az_solver;
	if (strcasecmp(name, "gmres") == 0) az_solver = AZ_gmres;
	else if (strcasecmp(name, "cg") == 0) az_solver = AZ_cg;
	else if (strcasecmp(name, "cgs") == 0) az_solver = AZ_cgs;
	else if (strcasecmp(name, "tfqmr") == 0) az_solver = AZ_tfqmr;
	else if (strcasecmp(name, "bicgstab") == 0) az_solver = AZ_bicgstab;
	else az_solver = AZ_gmres;

	aztec.SetAztecOption(AZ_solver, az_solver);
#endif
}

void AztecOOSolver::set_precond(const char *name)
{
#ifdef HAVE_AZTECOO
	int az_precond;
	if (strcasecmp(name, "none") == 0) az_precond = AZ_none;
	else if (strcasecmp(name, "jacobi") == 0) az_precond = AZ_Jacobi;
	else if (strcasecmp(name, "neumann") == 0) az_precond = AZ_Neumann;
	else if (strcasecmp(name, "least-squares") == 0) az_precond = AZ_ls;
	else az_precond = AZ_none;

	aztec.SetAztecOption(AZ_precond, az_precond);
#endif
}

bool AztecOOSolver::solve()
{
#ifdef HAVE_AZTECOO
	assert(m.size == rhs.size);

	// no output
	aztec.SetAztecOption(AZ_output, AZ_none);	// AZ_all | AZ_warnings | AZ_last | AZ_summary

#ifndef H2D_COMPLEX
	// setup the problem
	aztec.SetUserMatrix(m.mat);
	aztec.SetRHS(rhs.vec);
	Epetra_Vector x(*rhs.std_map);
	aztec.SetLHS(&x);

	if (pc != NULL) {
		Epetra_Operator *op = pc->get_obj();
		assert(op != NULL);		// can work only with Epetra_Operators
		aztec.SetPrecOperator(op);
	}

	// solve it
	aztec.Iterate(max_iters, tolerance);

	delete [] sln;
	sln = new scalar[m.size];
	memset(sln, 0, m.size * sizeof(scalar));

	// copy the solution into sln vector
	for (int i = 0; i < m.size; i++) sln[i] = x[i];
#else
	double c0r = 1.0, c0i = 0.0;
	double c1r = 0.0, c1i = 1.0;

	Epetra_Vector xr(*rhs.std_map);
	Epetra_Vector xi(*rhs.std_map);

	Komplex_LinearProblem kp(c0r, c0i, *m.mat, c1r, c1i, *m.mat_im, xr, xi, *rhs.vec, *rhs.vec_im);
	Epetra_LinearProblem *lp = kp.KomplexProblem();
	aztec.SetProblem(*lp);

	// solve it
	aztec.Iterate(max_iters, tolerance);

	kp.ExtractSolution(xr, xi);

	delete [] sln;
	sln = new scalar[m.size];
	memset(sln, 0, m.size * sizeof(scalar));

	// copy the solution into sln vector
	for (int i = 0; i < m.size; i++) sln[i] = scalar(xr[i], xi[i]);
#endif
	return true;
#else
	return false;
#endif
}

int AztecOOSolver::get_num_iters()
{
#ifdef HAVE_AZTECOO
	return aztec.NumIters();
#else
	return -1;
#endif
}
