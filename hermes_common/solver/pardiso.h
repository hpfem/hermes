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

#include "solver.h"
#include "../matrix.h"

class PardisoMatrix : public SparseMatrix 
{
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
  virtual int get_nnz() const;
  virtual double get_fill_in() const;

protected:
  // PARDISO specific data structures for storing the system matrix (CSR format).
  scalar *Ax;   // Matrix entries (row-wise). 
  int *Ai;      // Column indices of values in Ax.
  int *Ap;      // Index to Ax/Ai, where each row starts.
  int nnz;      // Number of non-zero entries (= Ap[size]).
  
  friend class PardisoLinearSolver;
};

class PardisoVector : public Vector 
{
public:
  PardisoVector();
  virtual ~PardisoVector();

  virtual void alloc(int ndofs);
  virtual void free();
  virtual scalar get(int idx) { return v[idx]; }
  virtual void zero();
  virtual void set(int idx, scalar y);
  virtual void add(int idx, scalar y);
  virtual void extract(scalar *v) const;
  virtual void add(int n, int *idx, scalar *y);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
  // PARDISO specific data structures for storing the rhs.
  scalar *v;

  friend class PardisoLinearSolver;
};

/// Encapsulation of PARDISO linear solver
///
/// @ingroup solvers
class HERMES_API PardisoLinearSolver : public LinearSolver
{
public:
  PardisoLinearSolver(PardisoMatrix *m, PardisoVector *rhs);
  virtual ~PardisoLinearSolver();

  virtual bool solve();

protected:
  PardisoMatrix *m;
  PardisoVector *rhs;
  
  bool inited;
  bool setup_factorization();
  bool check_status();
  
  // Variables needed on input to the Pardiso driver routine.
  int solver;     // Type of the solver (sparse direct or multi-recursive iterative).
  int mtype;      // Type of the matrix.
  int maxfct;     // Maximum number of numerical factorizations.
  int mnum;       // Which factorization to use.
  int msglvl;     // Controls output of statistical information.
  int err;        // Error flag.
  int nrhs;       // Number of right hand sides.
  int num_procs;  // Numbers of processors, value of OMP_NUM_THREADS.
  int phase;      // Calculation phase (symbolic factorization/numerical factorization/solution).
  
  // Internal solver memory pointer pt,
  //  32-bit: int pt[64]; 
  //  64-bit: long int pt[64]
  // void *pt[64] should be OK on both architectures.
  void *pt[64];
  
  // Pardiso control parameters. Consult Pardiso manual for interpretation of individual entries.
  int    iparm[64];
  double dparm[64];
};

#endif /* _PARDISO_SOLVER_H_*/
