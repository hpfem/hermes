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


#ifndef _PARDISO_SOLVER_H_
#define _PARDISO_SOLVER_H_

#include "../solver.h"
#include "../matrix.h"

class PardisoMatrix : public SparseMatrix {
public:
	PardisoMatrix();
	virtual ~PardisoMatrix();

	virtual void pre_add_ij(int row, int col);
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
	// PARDISO specific data structures for storing matrix, rhs
	int *Ap;
	int *Ai;
	scalar *Ax;

	static void insert_value(int *Ai, scalar *Ax, int Alen, int idx, scalar value);

	friend class PardisoLinearSolver;
};

class PardisoVector : public Vector {
public:
	PardisoVector();
	virtual ~PardisoVector();

	virtual void alloc(int ndofs);
	virtual void free();
	virtual scalar get(int idx) { return v[idx]; }
	virtual void zero();
	virtual void set(int idx, scalar y);
	virtual void add(int idx, scalar y);
	virtual void add(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
	scalar *v;

	friend class PardisoLinearSolver;
};

/// Encapsulation of PARDISO linear solver
///
/// @ingroup solvers
class PardisoLinearSolver : public LinearSolver {
public:
	PardisoLinearSolver(PardisoMatrix *m, PardisoVector *rhs);
	PardisoLinearSolver(LinearProblem *lp);
	virtual ~PardisoLinearSolver();

	virtual bool solve();

protected:
	PardisoMatrix *m;
	PardisoVector *rhs;
};

#endif /* _PARDISO_SOLVER_H_*/
