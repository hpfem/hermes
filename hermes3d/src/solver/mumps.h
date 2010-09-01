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

#ifndef _MUMPS_SOLVER_H_
#define _MUMPS_SOLVER_H_

#include "../solver.h"
#include "../matrix.h"

#ifdef WITH_MUMPS
extern "C" {
	#include <mumps_c_types.h>
#ifndef H3D_COMPLEX
	#include <dmumps_c.h>
#else
	#include <zmumps_c.h>
#endif
}
#else
struct ZMUMPS_COMPLEX {
	double r, i;
};
#endif


class MumpsMatrix : public SparseMatrix {
public:
	MumpsMatrix();
	virtual ~MumpsMatrix();

	virtual void alloc();
	virtual void free();
	virtual scalar get(int m, int n);
	virtual void zero();
	virtual void add(int m, int n, scalar v);
	virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
	virtual int get_matrix_size() const;
	virtual double get_fill_in() const;

protected:
	// MUMPS specific data structures for storing matrix, rhs
	int nnz;				// number of non-zero elements
	int *irn;				// row indices
	int *jcn;				// column indices
#ifndef H3D_COMPLEX
	scalar *a;				// matrix entries
#else
	ZMUMPS_COMPLEX *a;
#endif
	int *ap;
	int *ai;

	friend class MumpsSolver;
};


class MumpsVector : public Vector {
public:
	MumpsVector();
	virtual ~MumpsVector();

	virtual void alloc(int ndofs);
	virtual void free();
#ifndef H3D_COMPLEX
	virtual scalar get(int idx) { return v[idx]; }
#else
	virtual scalar get(int idx) { return std::complex<double>(v[idx].r, v[idx].i); }
#endif
	virtual void extract(scalar *v) const { memcpy(v, this->v, size * sizeof(scalar)); }
	virtual void zero();
	virtual void set(int idx, scalar y);
	virtual void add(int idx, scalar y);
	virtual void add(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
#ifndef H3D_COMPLEX
	scalar *v;
#else
	ZMUMPS_COMPLEX *v;
#endif

	friend class MumpsSolver;
};


/// Encapsulation of MUMPS linear solver
///
/// @ingroup solvers
class MumpsSolver : public LinearSolver {
public:
	MumpsSolver(MumpsMatrix *m, MumpsVector *rhs);
	MumpsSolver(LinearProblem *lp);
	virtual ~MumpsSolver();

	virtual bool solve();

protected:
	MumpsMatrix *m;
	MumpsVector *rhs;
};

#endif
