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

#ifndef _SUPERLU_SOLVER_H_
#define _SUPERLU_SOLVER_H_

#include "solver.h"
#include "../matrix.h"

#ifdef WITH_SUPERLU
  #if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    #ifdef SLU_MT
      #include <pdsp_defs.h>  
    #else
      #include <slu_ddefs.h>
    #endif
    
    #define SLU_DTYPE               SLU_D
    #define SLU_CREATE_DENSE_MATRIX dCreate_Dense_Matrix
    #define SLU_CREATE_CSC_MATRIX   dCreate_CompCol_Matrix
    #define SLU_PRINT_CSC_MATRIX    dPrint_CompCol_Matrix
    #define SLU_SOLVER_DRIVER       dgssvx
    #define SLU_SCALAR_MALLOC       doubleMalloc
    
    typedef scalar slu_scalar;
    #define SUPERLU_SCALAR(a) SCALAR(a)
  #else
    #ifdef SLU_MT
      #include <pzsp_defs.h>  
    #else
      #include <slu_zdefs.h>
    #endif
    
    #define SLU_DTYPE               SLU_Z
    #define SLU_CREATE_DENSE_MATRIX zCreate_Dense_Matrix
    #define SLU_CREATE_CSC_MATRIX   zCreate_CompCol_Matrix
    #define SLU_PRINT_CSC_MATRIX    zPrint_CompCol_Matrix
    #define SLU_SOLVER_DRIVER       zgssvx
    #define SLU_SCALAR_MALLOC       doublecomplexMalloc
    
    typedef doublecomplex slu_scalar;
    #define SUPERLU_SCALAR(a) a.r, a.i
  #endif
#else
  #if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    typedef scalar slu_scalar;
    #define SUPERLU_SCALAR(a) SCALAR(a)
  #else
    typedef struct  { double r, i; } slu_scalar;
    #define SUPERLU_SCALAR(a) a.r, a.i
  #endif
#endif


class SuperLUMatrix : public SparseMatrix {
public:
  SuperLUMatrix();
  virtual ~SuperLUMatrix();

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
  // SUPERLU specific data structures for storing the matrix (CSC format).
  slu_scalar *Ax; // Matrix entries (column-wise).
  int *Ai;        // Row indices of values in Ax.
  int *Ap;        // Index to Ax/Ai, where each column starts.
  int nnz;        // Number of non-zero entries (= Ap[size]).
  
  friend class SuperLUSolver;
};


class SuperLUVector : public Vector {
public:
  SuperLUVector();
  virtual ~SuperLUVector();

  virtual void alloc(int ndofs);
  virtual void free();
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
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
  // SUPERLU specific data structures for storing the rhs.
  slu_scalar *v;     // Vector entries.

friend class SuperLUSolver;
};


/// Encapsulation of SUPERLU linear solver
///
/// @ingroup solvers
class HERMES_API SuperLUSolver : public LinearSolver {
public:
  SuperLUSolver(SuperLUMatrix *m, SuperLUVector *rhs);
  virtual ~SuperLUSolver();

  virtual bool solve();
  
protected:
  SuperLUMatrix *m;       
  SuperLUVector *rhs;
  
  bool has_A, has_B;            // Have the native SuperLU matrices been created?
  bool inited;                  // Have the factorization structures been allocated?
  bool A_changed;               // Indicates that the system matrix has been changed
                                // internally during factorization or externally by
                                // the user.
                                
  bool check_status(int info);  // Check the status returned from the solver routine.
  
  // Deep copies of matrix and rhs data vectors (they may be changed by the solver driver,
  // hence we need a copy so that the original SuperLUMatrix/Vector is preserved).
  int *local_Ai, *local_Ap;
  slu_scalar *local_Ax, *local_rhs;
  
  bool prepare_factorization_structures();
  void free_factorization_structures();
  void free_matrix();
  void free_rhs();
  
#ifdef WITH_SUPERLU  
  SuperMatrix A, B;             // Native SuperLU representations of 'm' and 'rhs'.
  SuperMatrix L, U;             // L/U factors of A.
  double *R, *C;                // Row/column scaling factors of A.
  int *perm_r;                  // Row permutations from partial pivoting.
  int *perm_c;                  // Column permutations to reduce fill-in (=> matrix Pc)
  int *etree;                   // Elimination tree of Pc'*A'*A*Pc.
  char equed[1];                // Form of equilibration that was done on A.
  superlu_options_t options;    // Structure holding the input options.
#endif 
};

#endif
