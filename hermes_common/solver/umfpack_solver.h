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

#ifndef _UMFPACK_SOLVER_H_
#define _UMFPACK_SOLVER_H_

#include "solver.h"
#include "../matrix.h"

class HERMES_API UMFPackMatrix : public SparseMatrix {
public:
  UMFPackMatrix();
  UMFPackMatrix(int size);
  virtual ~UMFPackMatrix();

  virtual void alloc();
  virtual void free();
  virtual scalar get(int m, int n);
  virtual void zero();
  virtual void add(int m, int n, scalar v);
  virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
  virtual int get_matrix_size() const;
  int get_nnz() {return this->nnz;}
  virtual double get_fill_in() const;

  // Applies the matrix to vector_in and saves result to vector_out.
  void multiply(scalar* vector_in, scalar* vector_out);
  // Creates matrix in CSC format using size, nnz, and the three arrays.
  void create(int size, int nnz, int* ap, int* ai, scalar* ax);
  // Exposes pointers to the CSC arrays.
  int *get_Ap() {
      return this->Ap;
  }
  int *get_Ai() {
      return this->Ai;
  }
  scalar *get_Ax() {
      return this->Ax;
  }

protected:
  // UMFPack specific data structures for storing the system matrix (CSC format).
  scalar *Ax;   // Matrix entries (column-wise).
  int *Ai;      // Row indices of values in Ax.
  int *Ap;      // Index to Ax/Ai, where each column starts.
  int nnz;      // Number of non-zero entries (= Ap[size]).

  friend class UMFPackLinearSolver;
};

class HERMES_API UMFPackVector : public Vector {
public:
  UMFPackVector();
  UMFPackVector(int size);
  virtual ~UMFPackVector();

  virtual void alloc(int ndofs);
  virtual void free();
  virtual scalar get(int idx) { return v[idx]; }
  virtual void extract(scalar *v) const { memcpy(v, this->v, size * sizeof(scalar)); }
  virtual void zero();
  virtual void set(int idx, scalar y);
  virtual void add(int idx, scalar y);
  virtual void add(int n, int *idx, scalar *y);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

  scalar *get_c_array() {
      return this->v;
  }

protected:
  //UMFPack specific data structures for storing the rhs.
  scalar *v;
  friend class UMFPackLinearSolver;
};


/// Encapsulation of UMFPACK linear solver
///
/// @ingroup solvers
class HERMES_API UMFPackLinearSolver : public LinearSolver {
public:
  UMFPackLinearSolver(UMFPackMatrix *m, UMFPackVector *rhs);
  virtual ~UMFPackLinearSolver();

  virtual bool solve();
    
protected:
  UMFPackMatrix *m;
  UMFPackVector *rhs;
  
  // Reusable factorization information (A denotes matrix represented by the pointer 'm').
  void *symbolic; // Reordering of matrix A to reduce fill-in during factorization.
  void *numeric;  // LU factorization of matrix A.
  
  bool setup_factorization();
  void free_factorization_data();
};

#endif
