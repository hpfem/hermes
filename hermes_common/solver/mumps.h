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

#include "solver.h"
#include "../matrix.h"

#ifdef WITH_MUMPS
  extern "C" {
    #include <mumps_c_types.h>
  #if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    #include <dmumps_c.h>
    typedef scalar mumps_scalar;
    #define MUMPS_SCALAR(a) SCALAR(a)
    #define MUMPS_STRUCT    DMUMPS_STRUC_C
  #else
    #include <zmumps_c.h>
    typedef ZMUMPS_COMPLEX mumps_scalar;
    #define MUMPS_SCALAR(a) a.r, a.i  
    #define MUMPS_STRUCT    ZMUMPS_STRUC_C
  #endif
  }
  
  #ifdef WITH_MPI
    #include <mpi.h>
  #endif
  
#else
  #if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    typedef scalar mumps_scalar;
    #define MUMPS_SCALAR(a) SCALAR(a)
  #else
    typedef struct { double r, i; } mumps_scalar;
    #define MUMPS_SCALAR(a) a.r, a.i
  #endif
#endif


class MumpsMatrix : public SparseMatrix 
{
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
  // MUMPS specific data structures for storing the system matrix (CSC format).
  int nnz;          // Number of non-zero elements. 
  int *irn;         // Row indices.
  int *jcn;         // Column indices.
  mumps_scalar *Ax; // Matrix entries (column-wise).
  int *Ai;          // Row indices of values in Ax.
  int *Ap;          // Index to Ax/Ai, where each column starts.

  friend class MumpsSolver;
};


class MumpsVector : public Vector {
public:
  MumpsVector();
  virtual ~MumpsVector();

  virtual void alloc(int ndofs);
  virtual void free();
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  virtual scalar get(int idx) { return v[idx]; }
#else
  virtual scalar get(int idx) { return cplx(v[idx].r, v[idx].i); }
#endif
  virtual void extract(scalar *v) const { memcpy(v, this->v, size * sizeof(scalar)); }
  virtual void zero();
  virtual void set(int idx, scalar y);
  virtual void add(int idx, scalar y);
  virtual void add(int n, int *idx, scalar *y);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
  // MUMPS specific data structures for storing the rhs.
  mumps_scalar *v;

  friend class MumpsSolver;
};


/// Encapsulation of MUMPS linear solver
///
/// @ingroup solvers
class HERMES_API MumpsSolver : public LinearSolver {
public:
  MumpsSolver(MumpsMatrix *m, MumpsVector *rhs);
  virtual ~MumpsSolver();

  virtual bool solve();

protected:
  MumpsMatrix *m;
  MumpsVector *rhs;
  
  bool setup_factorization();

#ifdef WITH_MUMPS
  MUMPS_STRUCT  param;
  
  bool check_status();
  
  /// (Re)initialize a MUMPS instance.
  bool reinit();
  bool inited;
#endif
};

#endif
