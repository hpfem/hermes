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

#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

#include "../matrix.h"
#include "solver.h"

#ifdef WITH_PETSC
  #include <petsc.h>
  #include <petscmat.h>
  #include <petscvec.h>
  #include <petscksp.h>
#endif

/// Wrapper of PETSc matrix, to store matrices used with PETSc in its native format
///
class PetscMatrix : public SparseMatrix {
public:
  PetscMatrix();
  virtual ~PetscMatrix();

  virtual void alloc();
  virtual void free();
  virtual void finish();
  virtual scalar get(int m, int n);
  virtual void zero();
  virtual void add(int m, int n, scalar v);
  virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
  virtual int get_matrix_size() const;
  virtual double get_fill_in() const;

protected:
#ifdef WITH_PETSC
  Mat matrix;
#endif
  bool inited;

  friend class PetscLinearSolver;
};

/// Wrapper of PETSc vector, to store vectors used with PETSc in its native format
///
class PetscVector : public Vector {
public:
  PetscVector();
  virtual ~PetscVector();

  virtual void alloc(int ndofs);
  virtual void free();
  virtual void finish();
  virtual scalar get(int idx);
  virtual void extract(scalar *v) const;
  virtual void zero();
  virtual void set(int idx, scalar y);
  virtual void add(int idx, scalar y);
  virtual void add(int n, int *idx, scalar *y);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

protected:
#ifdef WITH_PETSC
  Vec vec;
#endif
  bool inited;

  friend class PetscLinearSolver;
};

/// Encapsulation of PETSc linear solver
///
/// @ingroup solvers
class HERMES_API PetscLinearSolver : public LinearSolver {
public:
  PetscLinearSolver(PetscMatrix *mat, PetscVector *rhs);
  virtual ~PetscLinearSolver();

  virtual bool solve();

protected:
  PetscMatrix *m;
  PetscVector *rhs;
};

#endif
