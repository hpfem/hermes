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


#include "petsc.h"
#include "../linear_problem.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>
#include <common/timer.h>

#define H3D_PETSC_NOT_COMPILED	"hermes3d was not built with PETSc support."


PetscMatrix::PetscMatrix() {
	_F_
#ifdef WITH_PETSC
	inited = false;
#else
	error(H3D_PETSC_NOT_COMPILED);
#endif
}

PetscMatrix::~PetscMatrix() {
	_F_
	free();
}

void PetscMatrix::alloc() {
	_F_
#ifdef WITH_PETSC
	assert(pages != NULL);

	// calc nnz
	int *nnz = new int[size];
	MEM_CHECK(nnz);

	// fill in nnz
	int aisize = get_num_indices();
	int *ai = new int[aisize];
	MEM_CHECK(ai);

	// sort the indices and remove duplicities, insert into ai
	int pos = 0;
	for (int i = 0; i < size; i++) {
		nnz[i] = sort_and_store_indices(pages[i], ai + pos, ai + aisize);
		pos += nnz[i];
	}
	delete [] pages; pages = NULL;
	delete [] ai;

	//
	MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 0, nnz, &matrix);
//	MatSetOption(matrix, MAT_ROW_ORIENTED);
//	MatSetOption(matrix, MAT_ROWS_SORTED);

	delete [] nnz;

	inited = true;
#endif
}

void PetscMatrix::free() {
	_F_
#ifdef WITH_PETSC
	if (inited) MatDestroy(matrix);
	inited = false;
#endif
}

void PetscMatrix::finish()
{
	_F_
#ifdef WITH_PETSC
	MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
#endif
}

scalar PetscMatrix::get(int m, int n)
{
	_F_
	scalar v = 0.0;
#ifdef WITH_PETSC
	MatGetValues(matrix, 1, &m, 1, &n, &v);
#endif
	return v;
}

void PetscMatrix::zero() {
	_F_
#ifdef WITH_PETSC
	MatZeroEntries(matrix);
#endif
}

void PetscMatrix::add(int m, int n, scalar v) {
	_F_
#ifdef WITH_PETSC
	if (v != 0.0 && m != H3D_DIRICHLET_DOF && n != H3D_DIRICHLET_DOF)		// ignore "dirichlet DOF"
		MatSetValue(matrix, m, n, (PetscScalar) v, ADD_VALUES);
#endif
}

void PetscMatrix::add(int m, int n, scalar **mat, int *rows, int *cols) {
	_F_
#ifdef WITH_PETSC
	// TODO: pass in just the block of the matrix without H3D_DIRICHLET_DOFs (so that can use MatSetValues directly without checking
	// row and cols for -1)
	for (int i = 0; i < m; i++)				// rows
		for (int j = 0; j < n; j++)			// cols
			add(rows[i], cols[j], mat[i][j]);
#endif
}

bool PetscMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat) {
	_F_
#ifdef WITH_PETSC
#endif
	return false;
}

int PetscMatrix::get_matrix_size() const {
	_F_
	return 0;
}

double PetscMatrix::get_fill_in() const {
	_F_
	return 0;
}

// PETSc vector //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscVector::PetscVector() {
	_F_
#ifdef WITH_PETSC
	inited = false;
#else
	error(H3D_PETSC_NOT_COMPILED);
#endif
}

PetscVector::~PetscVector() {
	_F_
	free();
}

void PetscVector::alloc(int n) {
	_F_
#ifdef WITH_PETSC
	free();
	size = n;
	VecCreateSeq(PETSC_COMM_SELF, size, &vec);
	inited = true;
#endif
}

void PetscVector::free() {
	_F_
#ifdef WITH_PETSC
	if (inited) VecDestroy(vec);
	inited = false;
#endif
}

void PetscVector::finish()
{
	_F_
#ifdef WITH_PETSC
	VecAssemblyBegin(vec);
	VecAssemblyEnd(vec);
#endif
}

scalar PetscVector::get(int idx) {
	_F_
	scalar y = 0;
#ifdef WITH_PETSC
	VecGetValues(vec, 1, &idx, &y);
#endif
	return y;
}

void PetscVector::extract(scalar *v) const {
	_F_
#ifdef WITH_PETSC
	int *idx = new int [size];
	for (int i = 0; i < size; i++) idx[i] = i;
	VecGetValues(vec, size, idx, (PetscScalar *) v);
	delete [] idx;
#endif
}

void PetscVector::zero() {
	_F_
#ifdef WITH_PETSC
	VecZeroEntries(vec);
#endif
}

void PetscVector::set(int idx, scalar y) {
	_F_
#ifdef WITH_PETSC
	if (idx >= 0) VecSetValue(vec, idx, (PetscScalar) y, INSERT_VALUES);
#endif
}

void PetscVector::add(int idx, scalar y) {
	_F_
#ifdef WITH_PETSC
	if (idx >= 0) VecSetValue(vec, idx, (PetscScalar) y, ADD_VALUES);
#endif
}

void PetscVector::add(int n, int *idx, scalar *y) {
	_F_
#ifdef WITH_PETSC
	for (int i = 0; i < n; i++)
		if (idx[i] >= 0) VecSetValue(vec, idx[i], (PetscScalar) y[i], ADD_VALUES);
#endif
}

bool PetscVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat) {
	_F_
#ifdef WITH_PETSC
#endif
	return false;
}

// PETSc linear solver ///////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscLinearSolver::PetscLinearSolver(PetscMatrix *mat, PetscVector *rhs)
	: LinearSolver(), m(mat), rhs(rhs)
{
	_F_
#ifdef WITH_PETSC
#else
	warning(H3D_PETSC_NOT_COMPILED);
	exit(128);
#endif
}

PetscLinearSolver::PetscLinearSolver(LinearProblem *lp)
	: LinearSolver(lp)
{
#ifdef WITH_PETSC
	m = new PetscMatrix;
	rhs = new PetscVector;
#else
	warning(H3D_PETSC_NOT_COMPILED);
	exit(128);
#endif
}


PetscLinearSolver::~PetscLinearSolver() {
	_F_
#ifdef WITH_PETSC
	if (lp != NULL) {
		delete m;
		delete rhs;
	}
#endif
}

bool PetscLinearSolver::solve() {
	_F_
#ifdef WITH_PETSC
	assert(m != NULL);
	assert(rhs != NULL);

	if (lp != NULL)
		lp->assemble(m, rhs);
	assert(m->size == rhs->size);

	PetscErrorCode ec;
	KSP ksp;
	Vec x;

	Timer tmr;
	tmr.start();

	KSPCreate(PETSC_COMM_WORLD, &ksp);

	KSPSetOperators(ksp, m->matrix, m->matrix, DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);
	VecDuplicate(rhs->vec, &x);

	ec = KSPSolve(ksp, rhs->vec, x);
	if (ec) return false;

	tmr.stop();
	time = tmr.get_seconds();

	// allocate memory for solution vector
	delete [] sln;
	sln = new scalar [m->size];
	MEM_CHECK(sln);
	memset(sln, 0, m->size * sizeof(scalar));

	// index map vector (basic serial code uses the map sln[i] = x[i] for all dofs.
	int *idx = new int [m->size];
	MEM_CHECK(idx);
	for (int i = 0; i < m->size; i++) idx[i] = i;

	// copy solution to the output solution vector
	VecGetValues(x, m->size, idx, (PetscScalar *) sln);
	delete [] idx;

	KSPDestroy(ksp);
	VecDestroy(x);

	return true;
#else
	return false;
#endif
}
