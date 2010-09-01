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

#ifdef WITH_UMFPACK
extern "C" {
#include <umfpack.h>
}
#endif

#include "umfpack.h"
#include "../linear_problem.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/utils.h>
#include <common/callstack.h>
#include <common/timer.h>

UMFPackMatrix::UMFPackMatrix() {
	_F_
	Ap = NULL;
	Ai = NULL;
	Ax = NULL;
}

UMFPackMatrix::~UMFPackMatrix() {
	_F_
	free();
}

void UMFPackMatrix::alloc() {
	_F_
	assert(pages != NULL);

	// initialize the arrays Ap and Ai
	Ap = new int [size + 1];
	MEM_CHECK(Ap);
	int aisize = get_num_indices();
	Ai = new int [aisize];
	MEM_CHECK(Ai);

	// sort the indices and remove duplicities, insert into Ai
	int i, pos = 0;
	for (i = 0; i < size; i++) {
		Ap[i] = pos;
		pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
	}
	Ap[i] = pos;

	delete [] pages;
	pages = NULL;

	Ax = new scalar [Ap[size]];
	MEM_CHECK(Ax);
	memset(Ax, 0, sizeof(scalar) * Ap[size]);
}

void UMFPackMatrix::free() {
	_F_
	delete [] Ap; Ap = NULL;
	delete [] Ai; Ai = NULL;
	delete [] Ax; Ax = NULL;
}

scalar UMFPackMatrix::get(int m, int n)
{
	_F_

	// bin search the value
	register int lo = Ap[m], hi = Ap[m + 1], mid;
	while (1) {
		mid = (lo + hi) >> 1;

		if (n < Ai[mid]) hi = mid - 1;
		else if (n > Ai[mid]) lo = mid + 1;
		else break;

		if (lo > hi) return 0.0;		// entry not set -> e.i. it is zero
	}

	return Ax[mid];
}

void UMFPackMatrix::zero() {
	_F_
	memset(Ax, 0, sizeof(scalar) * Ap[size]);
}

void UMFPackMatrix::add(int m, int n, scalar v) {
	_F_
	if (v != 0.0 && m != H3D_DIRICHLET_DOF && n != H3D_DIRICHLET_DOF)		// ignore dirichlet DOFs
		insert_value(Ai + Ap[n], Ax + Ap[n], Ap[n + 1] - Ap[n], m, v);
}

void UMFPackMatrix::add(int m, int n, scalar **mat, int *rows, int *cols) {
	_F_
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)			// cols
			add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool UMFPackMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
	_F_
	switch (fmt) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
			for (int j = 0; j < size; j++)
				for (int i = Ap[j]; i < Ap[j + 1]; i++)
					fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
			fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

			return true;

		case DF_HERMES_BIN: {
			hermes_fwrite("H2DX\001\000\000\000", 1, 8, file);
			int ssize = sizeof(scalar);
			int nnz = Ap[size];
			hermes_fwrite(&ssize, sizeof(int), 1, file);
			hermes_fwrite(&size, sizeof(int), 1, file);
			hermes_fwrite(&nnz, sizeof(int), 1, file);
			hermes_fwrite(Ap, sizeof(int), size + 1, file);
			hermes_fwrite(Ai, sizeof(int), nnz, file);
			hermes_fwrite(Ax, sizeof(scalar), nnz, file);
			return true;
		}

		case DF_PLAIN_ASCII:
			EXIT(H3D_ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
}

int UMFPackMatrix::get_matrix_size() const {
	_F_
	assert(Ap != NULL);
	return (sizeof(int) + sizeof(scalar)) * (Ap[size] + size);
}

double UMFPackMatrix::get_fill_in() const {
	_F_
	return Ap[size] / (double) (size * size);
}

void UMFPackMatrix::insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value) {
	_F_
	if (idx >= 0) {
		register int lo = 0, hi = Alen - 1, mid;

		while (1) {
			mid = (lo + hi) >> 1;

			if (idx < Ai[mid]) hi = mid - 1;
			else if (idx > Ai[mid]) lo = mid + 1;
			else break;

			if (lo > hi) EXIT("Sparse matrix entry not found.");
		}

		Ax[mid] += value;
	}
}


// UMFPackVector ///////

UMFPackVector::UMFPackVector() {
	_F_
	v = NULL;
	size = 0;
}

UMFPackVector::~UMFPackVector() {
	_F_
	free();
}

void UMFPackVector::alloc(int n) {
	_F_
	free();
	size = n;
	v = new scalar [n];
	MEM_CHECK(v);
	zero();
}

void UMFPackVector::zero() {
	_F_
	memset(v, 0, size * sizeof(scalar));
}

void UMFPackVector::free() {
	_F_
	delete [] v;
	v = NULL;
	size = 0;
}

void UMFPackVector::set(int idx, scalar y) {
	_F_
	if (idx >= 0) v[idx] = y;
}

void UMFPackVector::add(int idx, scalar y) {
	_F_
	if (idx >= 0) v[idx] += y;
}

void UMFPackVector::add(int n, int *idx, scalar *y) {
	_F_
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) v[idx[i]] += y[i];
}

bool UMFPackVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
	_F_
	switch (fmt) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
			for (int i = 0; i < size; i++)
				fprintf(file, SCALAR_FMT "\n", SCALAR(v[i]));
			fprintf(file, " ];\n");
			return true;

		case DF_HERMES_BIN: {
			hermes_fwrite("H2DR\001\000\000\000", 1, 8, file);
			int ssize = sizeof(scalar);
			hermes_fwrite(&ssize, sizeof(int), 1, file);
			hermes_fwrite(&size, sizeof(int), 1, file);
			hermes_fwrite(v, sizeof(scalar), size, file);
			return true;
		}

		case DF_PLAIN_ASCII:
			EXIT(H3D_ERR_NOT_IMPLEMENTED);
			return false;

		default:
			return false;
	}
}

// UMFPack solver //////

#ifndef H3D_COMPLEX
// real case
#define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)		umfpack_di_symbolic(m, n, Ap, Ai, Ax, S, C, I)
#define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)			umfpack_di_numeric(Ap, Ai, Ax, S, N, C, I)
#define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I)	umfpack_di_solve(sys, Ap, Ai, Ax, X, B, N, C, I)
#define umfpack_free_symbolic							umfpack_di_free_symbolic
#define umfpack_free_numeric							umfpack_di_free_numeric
#define umfpack_defaults								umfpack_di_defaults
#else
// macros for calling complex UMFPACK in packed-complex mode
#define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)		umfpack_zi_symbolic(m, n, Ap, Ai, (double *) (Ax), NULL, S, C, I)
#define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)			umfpack_zi_numeric(Ap, Ai, (double *) (Ax), NULL, S, N, C, I)
#define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I)	umfpack_zi_solve(sys, Ap, Ai, (double *) (Ax), NULL, (double *) (X), NULL, (double *) (B), NULL, N, C, I)
#define umfpack_free_symbolic							umfpack_di_free_symbolic
#define umfpack_free_numeric							umfpack_zi_free_numeric
#define umfpack_defaults								umfpack_zi_defaults
#endif


UMFPackLinearSolver::UMFPackLinearSolver(UMFPackMatrix *m, UMFPackVector *rhs)
	: LinearSolver(), m(m), rhs(rhs)
{
	_F_
#ifdef WITH_UMFPACK
#else
	error("hermes3d was not built with UMFPACK support.");
#endif
}

UMFPackLinearSolver::UMFPackLinearSolver(LinearProblem *lp)
	: LinearSolver(lp)
{
	_F_
#ifdef WITH_UMFPACK
	m = new UMFPackMatrix;
	rhs = new UMFPackVector;
#else
	error("hermes3d was not built with UMFPACK support.");
#endif
}

UMFPackLinearSolver::~UMFPackLinearSolver() {
	_F_
#ifdef WITH_UMFPACK
	if (lp != NULL) {
		delete m;
		delete rhs;
	}
#endif
}

#ifdef WITH_UMFPACK

static void check_status(const char *fn_name, int status) {
	_F_
	switch (status) {
		case UMFPACK_OK: break;
		case UMFPACK_WARNING_singular_matrix:       warning("%s: singular matrix!", fn_name); break;
		case UMFPACK_ERROR_out_of_memory:           warning("%s: out of memory!", fn_name); break;
		case UMFPACK_ERROR_argument_missing:        warning("%s: argument missing", fn_name); break;
		case UMFPACK_ERROR_invalid_Symbolic_object: warning("%s: invalid Symbolic object", fn_name); break;
		case UMFPACK_ERROR_invalid_Numeric_object:  warning("%s: invalid Numeric object", fn_name); break;
		case UMFPACK_ERROR_different_pattern:       warning("%s: different pattern", fn_name); break;
		case UMFPACK_ERROR_invalid_system:          warning("%s: invalid system", fn_name); break;
		case UMFPACK_ERROR_n_nonpositive:           warning("%s: n nonpositive", fn_name); break;
		case UMFPACK_ERROR_invalid_matrix:          warning("%s: invalid matrix", fn_name); break;
		case UMFPACK_ERROR_internal_error:          warning("%s: internal error", fn_name); break;
		default:                                    warning("%s: unknown error (%d)", fn_name, status); break;
	}
}

#endif

bool UMFPackLinearSolver::solve() {
	_F_
#ifdef WITH_UMFPACK
	assert(m != NULL);
	assert(rhs != NULL);

	if (lp != NULL)
		lp->assemble(m, rhs);
	assert(m->size == rhs->size);

	Timer tmr;
	tmr.start();

	void *symbolic, *numeric;
	int status;

	status = umfpack_symbolic(m->size, m->size, m->Ap, m->Ai, m->Ax, &symbolic, NULL, NULL);
	if (status != UMFPACK_OK) {
		check_status("umfpack_di_symbolic", status);
		return false;
	}
	if (symbolic == NULL) EXIT("umfpack_di_symbolic error: symbolic == NULL");

	status = umfpack_numeric(m->Ap, m->Ai, m->Ax, symbolic, &numeric, NULL, NULL);
	if (status != UMFPACK_OK) {
		check_status("umfpack_di_numeric", status);
		return false;
	}
	if (numeric == NULL) EXIT("umfpack_di_numeric error: numeric == NULL");

	delete [] sln;
	sln = new scalar[m->size];
	MEM_CHECK(sln);
	memset(sln, 0, m->size * sizeof(scalar));

	status = umfpack_solve(UMFPACK_A, m->Ap, m->Ai, m->Ax, sln, rhs->v, numeric, NULL, NULL);
	if (status != UMFPACK_OK) {
		check_status("umfpack_di_solve", status);
		return false;
	}

	tmr.stop();
	time = tmr.get_seconds();

	umfpack_free_symbolic(&symbolic);
	umfpack_free_numeric(&numeric);

	return true;
#else
	return false;
#endif
}
