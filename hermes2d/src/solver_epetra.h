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

#ifndef __H2D_SOLVER_EPETRA_H_
#define __H2D_SOLVER_EPETRA_H_

#include "config.h"
#include "itersolver.h"
#include "matrix_old.h"

#ifdef HAVE_EPETRA
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#endif


class EpetraMatrix : public SparseMatrix
{
public:
	EpetraMatrix();
#ifdef HAVE_EPETRA
	EpetraMatrix(Epetra_RowMatrix &mat);
#endif
	virtual ~EpetraMatrix();

	virtual void prealloc(int n);
	virtual void pre_add_ij(int row, int col);
	virtual void finish();

	virtual void alloc();
	virtual void free();
	virtual scalar get(int m, int n);
	virtual int get_num_row_entries(int row);
	virtual void extract_row_copy(int row, int len, int &n_entries, double *vals, int *idxs);
	virtual void zero();
	virtual void add(int m, int n, scalar v);
	virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
	virtual int get_matrix_size() const;

protected:
#ifdef HAVE_EPETRA
	Epetra_BlockMap *std_map;
	Epetra_CrsGraph *grph;
	Epetra_CrsMatrix *mat;
#ifdef H2D_COMPLEX
	Epetra_CrsMatrix *mat_im;		// imaginary part of the matrix, mat holds the real part
#endif
	bool owner;
#endif

	friend class AmesosSolver;
	friend class AztecOOSolver;
	friend class NoxSolver;
	friend class IfpackPrecond;
	friend class MlPrecond;
};


class EpetraVector : public _Vector
{
public:
	EpetraVector();
#ifdef HAVE_EPETRA
	EpetraVector(const Epetra_Vector &v);
#endif
	virtual ~EpetraVector();

	virtual void alloc(int ndofs);
	virtual void free();
#ifdef HAVE_EPETRA
	virtual scalar get(int idx) { return (*vec)[idx]; }
#ifndef H2D_COMPLEX
	virtual void extract(scalar *v) const { vec->ExtractCopy(v); }
#else
	virtual void extract(scalar *v) const { }
#endif
#else
	virtual scalar get(int idx) { return 0.0; }
	virtual void extract(scalar *v) const { }
#endif
	virtual void zero();
	virtual void set(int idx, scalar y);
	virtual void add(int idx, scalar y);
	virtual void add(int n, int *idx, scalar *y);
	virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
#ifdef HAVE_EPETRA
	Epetra_BlockMap *std_map;
	Epetra_Vector *vec;
#ifdef H2D_COMPLEX
	Epetra_Vector *vec_im;		// imaginary part of the vector, vec holds the real part
#endif
	bool owner;
#endif

	friend class AztecOOSolver;
	friend class NoxSolver;
};

#endif
