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
#include "pardiso.h"
#include "../linear_problem.h"

#include <common/trace.h>
#include <common/error.h>
#include <common/utils.h>
#include <common/callstack.h>

#ifdef __cplusplus
extern "C" {
#endif

extern int pardisoinit_(void *, int *, int *);

extern int
    pardiso_(void *, int *, int *, int *, int *, int *,
             scalar *, int *, int *, int *, int *, int *, int *, scalar *, scalar *, int *);

#define PARDISOINIT pardisoinit_
#define PARDISO pardiso_

#ifdef __cplusplus
}
#endif

PardisoMatrix::PardisoMatrix() {
	_F_
	Ap = NULL;
	Ai = NULL;
	Ax = NULL;
}

PardisoMatrix::~PardisoMatrix() {
	_F_
	free();
}

void PardisoMatrix::pre_add_ij(int row, int col) {
	_F_
	SparseMatrix::pre_add_ij(col, row);
}

void PardisoMatrix::alloc() {
	_F_
	assert(pages != NULL);

	// initialize the arrays Ap and Ai
	Ap = new int[size + 1];
	MEM_CHECK(Ap);
	int aisize = get_num_indices();
	Ai = new int[aisize];
	MEM_CHECK(Ai);

	// sort the indices and remove duplicities, insert into Ai
	int i, pos = 0;
	for (i = 0; i < size; i++) {
		Ap[i] = pos;
		pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
	}
	Ap[i] = pos;

	delete[] pages;
	pages = NULL;

	Ax = new scalar[Ap[size]];
	MEM_CHECK(Ax);
	zero();
}

void PardisoMatrix::free() {
	_F_
	delete [] Ap; Ap = NULL;
	delete [] Ai; Ai = NULL;
	delete [] Ax; Ax = NULL;
}

scalar PardisoMatrix::get(int m, int n)
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

void PardisoMatrix::zero() {
	_F_
    memset(Ax, 0, sizeof(scalar) * Ap[size]);
}

void PardisoMatrix::add(int m, int n, scalar v) {
	_F_
	insert_value(Ai + Ap[m], Ax + Ap[m], Ap[m + 1] - Ap[m], n, v);
}

void PardisoMatrix::add(int m, int n, scalar **mat, int *rows, int *cols) {
	_F_
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)			// cols
			if (mat[i][j] != 0.0 && rows[i] != H3D_DIRICHLET_DOF && cols[j] != H3D_DIRICHLET_DOF)		// ignore dirichlet DOFs
				add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool PardisoMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
	_F_
	switch (fmt) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
			for (int j = 0; j < size; j++)
				for (int i = Ap[j]; i < Ap[j + 1]; i++)
					fprintf(file, "%d %d " SCALAR_FMT ";\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
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

int PardisoMatrix::get_matrix_size() const {
	_F_
	assert(Ap != NULL);
	return (sizeof(int) + sizeof(scalar)) * (Ap[size] + size);
}

double PardisoMatrix::get_fill_in() const {
	_F_
	return Ap[size] / (double) (size * size);
}

void PardisoMatrix::insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value) {
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


// PardisoVector ///////

PardisoVector::PardisoVector() {
	_F_
	v = NULL;
	size = 0;
}

PardisoVector::~PardisoVector() {
	_F_
	free();
}

void PardisoVector::alloc(int n) {
	_F_
	free();
	v = new scalar[n];
	MEM_CHECK(v);
	size = n;
	zero();
}

void PardisoVector::zero() {
	_F_
	memset(v, 0, size * sizeof(scalar));
}

void PardisoVector::free() {
	_F_
	delete [] v;
	v = NULL;
	size = 0;
}

void PardisoVector::set(int idx, scalar y) {
	_F_
	if (idx >= 0) v[idx] = y;
}

void PardisoVector::add(int idx, scalar y) {
	_F_
	if (idx >= 0) v[idx] += y;
}

void PardisoVector::add(int n, int *idx, scalar *y) {
	_F_
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) v[idx[i]] += y[i];
}

bool PardisoVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
	_F_
	switch (fmt) {
		case DF_MATLAB_SPARSE:
			fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
			for (int i = 0; i < this->size; i++)
				fprintf(file, SCALAR_FMT ";\n", SCALAR(v[i]));
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

// PARDISO solver //////

PardisoLinearSolver::PardisoLinearSolver(PardisoMatrix *m, PardisoVector *rhs)
	: LinearSolver(), m(m), rhs(rhs)
{
	_F_
#ifdef WITH_PARDISO
#else
	warning("hermes3d was not built with Pardiso support.");
	exit(128);
#endif
}

PardisoLinearSolver::PardisoLinearSolver(LinearProblem *lp)
	: LinearSolver(lp)
{
	_F_
#ifdef WITH_PARDISO
	m = new PardisoMatrix;
	rhs = new PardisoVector;
#else
	warning("hermes3d was not built with Pardiso support.");
	exit(128);
#endif
}

PardisoLinearSolver::~PardisoLinearSolver() {
	_F_
#ifdef WITH_PARDISO
	if (lp != NULL) {
		delete m;
		delete rhs;
	}
#endif
}

bool PardisoLinearSolver::solve() {
	_F_
#ifdef WITH_PARDISO
	assert(m != NULL);
	assert(rhs != NULL);

	if (lp != NULL)
		lp->assemble(m, rhs);
	assert(m->size == rhs->size);

	bool res = true;
	int n = m->size;

	try {
		// Numbers of processors, value of OMP_NUM_THREADS
		int num_procs;
		char *var = getenv("OMP_NUM_THREADS");
		if (var != NULL) sscanf(var, "%d", &num_procs);
		else num_procs = 1;

		int mtype = 11;		// Real unsymmetric matrix
		int nrhs = 1;		// Number of right hand sides
		int nnz = m->Ap[n];	// The number of nonzero elements

		// Internal solver memory pointer pt,
		// 32-bit: int pt[64]; 64-bit: long int pt[64]
		// or void *pt[64] should be OK on both architectures
		void *pt[64];
		// Pardiso control parameters.
		int iparm[64];
		int maxfct, mnum, phase, err, msglvl;
		// Auxiliary variables.
		scalar ddum; // Double dummy
		int idum; // Integer dummy.

		iparm[2] = num_procs;

		maxfct = 1;		// Maximum number of numerical factorizations.
		mnum = 1;		// Which factorization to use.
		msglvl = 0;		// Do not print statistical information
		err = 0;		// Initialize error flag

		// Convert matrix from 0-based C-notation to Fortran 1-based notation.
		for (int i = 0; i < n + 1; i++) m->Ap[i] += 1;
		for (int i = 0; i < nnz; i++) m->Ai[i] += 1;

		Timer tmr;
		tmr.start();

		// Setup Pardiso control parameters.
		PARDISOINIT(pt, &mtype, iparm);

		// .. Reordering and Symbolic Factorization. This step also allocates
		// all memory that is necessary for the factorization.
		phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err);
		if (err != 0) {
			// ERROR during symbolic factorization: err
			throw ERR_FAILURE;
		}

		// .. Numerical factorization.
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err);
		if (err != 0) {
			// ERROR during numerical factorization: err
			throw ERR_FAILURE;
		}

		// .. Back substitution and iterative refinement.
		delete [] sln;
		sln = new scalar[m->size];
		MEM_CHECK(sln);
		memset(sln, 0, (m->size) * sizeof(scalar));

		phase = 33;
		iparm[7] = 1; // Max numbers of iterative refinement steps.
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, rhs->v, sln, &err);
		if (err != 0) {
			// ERROR during solution: err
			throw ERR_FAILURE;
		}

		tmr.stop();
		time = tmr.get_seconds();

		//  Convert matrix back to 0-based C-notation.
		for (int i = 0; i < n + 1; i++) m->Ap[i] -= 1;
		for (int i = 0; i < nnz; i++) m->Ai[i] -= 1;

		// .. Termination and release of memory.
		phase = -1; // Release internal memory.
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err);
	}
	catch (int e) {
		error = e;
		res = false;
	}

	return res;
#else
	return false;
#endif
}
