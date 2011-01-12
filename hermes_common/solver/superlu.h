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
            
      #define SLU_GSEQU         dgsequ
      #define SLU_LAQGS         dlaqgs
      #define SLU_SP_COLORDER   sp_colorder
      #define SLU_GSTRF         pdgstrf
      #define SLU_PIVOT_GROWTH  dPivotGrowth
      #define SLU_LANGS         dlangs
      #define SLU_GSCON         dgscon
      #define SLU_GSTRS         dgstrs
      #define SLU_GSRFS         dgsrfs
      #define SLU_LAMCH_        dlamch_
      #define SLU_QUERY_SPACE   superlu_dQuerySpace
      
      #define SLU_MULT(a,b)     a *= b
      
    #else
      #include <slu_ddefs.h>
    #endif
    
    #define SLU_DTYPE               SLU_D
    #define SLU_CREATE_DENSE_MATRIX dCreate_Dense_Matrix
    #define SLU_CREATE_CSC_MATRIX   dCreate_CompCol_Matrix
    #define SLU_PRINT_CSC_MATRIX    dPrint_CompCol_Matrix
    #define SLU_SCALAR_MALLOC       doubleMalloc
    
    #ifndef SLU_MT
      #define SLU_SOLVER_DRIVER     dgssvx
    #endif
    
    typedef scalar slu_scalar;
    #define SUPERLU_SCALAR(a) SCALAR(a)
  #else
    #ifdef SLU_MT
      #include <pzsp_defs.h>
      
      #define SLU_GSEQU         zgsequ
      #define SLU_LAQGS         zlaqgs
      #define SLU_SP_COLORDER   sp_colorder
      #define SLU_GSTRF         pzgstrf
      #define SLU_PIVOT_GROWTH  zPivotGrowth
      #define SLU_LANGS         zlangs
      #define SLU_GSCON         zgscon
      #define SLU_GSTRS         zgstrs
      #define SLU_GSRFS         zgsrfs
      #define SLU_LAMCH_        dlamch_
      #define SLU_QUERY_SPACE   superlu_zQuerySpace
      
      #define SLU_MULT(a,b)     zd_mult(&a, &a, b)
      
    #else
      #include <slu_zdefs.h>
    #endif
    
    #define SLU_DTYPE               SLU_Z
    #define SLU_CREATE_DENSE_MATRIX zCreate_Dense_Matrix
    #define SLU_CREATE_CSC_MATRIX   zCreate_CompCol_Matrix
    #define SLU_PRINT_CSC_MATRIX    zPrint_CompCol_Matrix
    #define SLU_SCALAR_MALLOC       doublecomplexMalloc
    
    #ifndef SLU_MT
      #define SLU_SOLVER_DRIVER     zgssvx
    #endif
    
    typedef doublecomplex slu_scalar;
    #define SUPERLU_SCALAR(a) a.r, a.i
  #endif
  
  #ifdef SLU_MT
    typedef superlumt_options_t       slu_options_t;
    typedef Gstat_t                   slu_stat_t;
    typedef superlu_memusage_t        slu_memusage_t;
    #define SLU_DESTROY_L             Destroy_SuperNode_SCP 
    #define SLU_DESTROY_U             Destroy_CompCol_NCP
    #define SLU_INIT_STAT(stat_ptr)   StatAlloc(m->size, options.nprocs,\
                                                options.panel_size, options.relax,\
                                                stat_ptr);\
                                      StatInit(m->size, options.nprocs, stat_ptr)  
    #define SLU_PRINT_STAT(stat_ptr)  PrintStat(stat_ptr)    
    
    void slu_mt_solver_driver(slu_options_t *options, SuperMatrix *A, 
                              int *perm_c, int *perm_r, SuperMatrix *AC,
                              equed_t *equed, double *R, double *C,
                              SuperMatrix *L, SuperMatrix *U,
                              SuperMatrix *B, SuperMatrix *X, 
                              double *recip_pivot_growth, double *rcond, 
                              double *ferr, double *berr, 
                              slu_stat_t *stat, slu_memusage_t *memusage,
                              int *info);        
  #else
    typedef superlu_options_t         slu_options_t;
    typedef SuperLUStat_t             slu_stat_t;
    typedef mem_usage_t               slu_memusage_t;
    #define SLU_DESTROY_L             Destroy_SuperNode_Matrix
    #define SLU_DESTROY_U             Destroy_CompCol_Matrix
    #define SLU_INIT_STAT(stat_ptr)   StatInit(stat_ptr)
    #define SLU_PRINT_STAT(stat_ptr)  StatPrint(stat_ptr)
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
  virtual void add_to_diagonal(scalar v);
  virtual void add(int m, int n, scalar **mat, int *rows, int *cols);
  virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
  virtual int get_matrix_size() const;
  virtual int get_nnz() const;
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
  virtual void change_sign();
  virtual void set(int idx, scalar y);
  virtual void add(int idx, scalar y);
  virtual void add(int n, int *idx, scalar *y);
  virtual void add_vector(Vector* vec) {
    assert(this->length() == vec->length());
    for (int i = 0; i < this->length(); i++) this->add(i, vec->get(i));
  };
  virtual void add_vector(scalar* vec) {
    for (int i = 0; i < this->length(); i++) this->add(i, vec[i]);
  };
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
  
  bool setup_factorization();
  void free_factorization_data();
  void free_matrix();
  void free_rhs();
  
#ifdef WITH_SUPERLU  
  SuperMatrix A, B;             // Native SuperLU representations of 'm' and 'rhs'.
  SuperMatrix L, U;             // L/U factors of A.
  double *R, *C;                // Row/column scaling factors of A.
  int *perm_r;                  // Row permutations from partial pivoting.
  int *perm_c;                  // Column permutations to reduce fill-in (=> matrix Pc)
  int *etree;                   // Elimination tree of Pc'*A'*A*Pc.
  slu_options_t options;        // Structure holding the input options for the solver.
  
  
  #ifndef SLU_MT
    char equed[1];              // Form of equilibration that was done on A.
  #else  
    equed_t equed;              // Form of equilibration that was done on A.
    SuperMatrix AC;             // Matrix A permuted by perm_c.
  #endif  
#endif 
};

#endif
