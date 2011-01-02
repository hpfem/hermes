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

#include "superlu.h"
#include "../trace.h"
#include "../error.h"
#include "../utils.h"
#include "../callstack.h"

// Binary search for the location of a particular CSC/CSR matrix entry.
//
// Typically, we search for the index into Ax that corresponds to a given 
// row (CSC) or column (CSR) ('idx') among indices of nonzero values in 
// a particular column (CSC) or row (CSR) ('Ai').
//
static int find_position(int *Ai, int Alen, int idx)
{
  _F_
  assert (idx >= 0);
  
  register int lo = 0, hi = Alen - 1, mid;
  
  while (1) 
  {
    mid = (lo + hi) >> 1;
    
    if (idx < Ai[mid]) hi = mid - 1;
    else if (idx > Ai[mid]) lo = mid + 1;
    else break;
    
    // Sparse matrix entry not found (raise an error when trying to add 
    // value to this position, return 0 when obtaining value there).
    if (lo > hi) mid = -1;
  }
  return mid;
}

SuperLUMatrix::SuperLUMatrix()
{
  _F_
  size = 0; nnz = 0;
  Ax = NULL;
  Ap = NULL;
  Ai = NULL;
}

SuperLUMatrix::~SuperLUMatrix()
{
  _F_
  this->free();
}

void SuperLUMatrix::alloc()
{
  _F_
  assert(pages != NULL);
  assert(size > 0);
  
  // Initialize the arrays Ap and Ai.
  Ap = new int [size + 1];
  MEM_CHECK(Ap);
  int aisize = get_num_indices();
  Ai = new int [aisize];
  MEM_CHECK(Ai);
  
  // sort the indices and remove duplicities, insert into Ai
  int i, pos = 0;
  for (i = 0; i < size; i++) {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;
  
  delete [] pages;
  pages = NULL;
  
  nnz = Ap[size];

  Ax = new slu_scalar [nnz];
  memset(Ax, 0, sizeof(slu_scalar) * nnz);
}

void SuperLUMatrix::free()
{
  _F_
  nnz = 0;
  delete [] Ap; Ap = NULL;
  delete [] Ai; Ai = NULL;
  delete [] Ax; Ax = NULL;
}

scalar SuperLUMatrix::get(int m, int n)
{
  _F_
  // Find m-th row in the n-th column.
  int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
  // Return 0 if the entry has not been found.
  if (mid < 0) return 0.0;
  // Otherwise, add offset to the n-th column and return the value.
  if (mid >= 0) mid += Ap[n];
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  return Ax[mid];
#else
  return cplx(Ax[mid].r, Ax[mid].i);
#endif
}

void SuperLUMatrix::zero()
{
  _F_
  memset(Ax, 0, sizeof(slu_scalar) * nnz);
}

void SuperLUMatrix::add(int m, int n, scalar v)
{
  _F_
  if (v != 0.0 && m >= 0 && n >= 0) // ignore dirichlet DOFs
  {   
    // Find m-th row in the n-th column.
    int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
    // Make sure we are adding to an existing non-zero entry.
    if (pos < 0) 
      error("Sparse matrix entry not found");
    // Add offset to the n-th column.
    pos += Ap[n];
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    Ax[pos] += v;
#else
    Ax[pos].r += v.real();
    Ax[pos].i += v.imag();
#endif
  }
}

void SuperLUMatrix::add(int m, int n, scalar **mat, int *rows, int *cols)
{
  _F_
  for (int i = 0; i < m; i++)       // rows
    for (int j = 0; j < n; j++)     // cols
      add(rows[i], cols[j], mat[i][j]);
}

/// Save matrix and right-hand side to a file.
///
bool SuperLUMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  // TODO
  switch (fmt) 
  {      
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
      for (int j = 0; j < size; j++)
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)          
          fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, SUPERLU_SCALAR(Ax[i]));
#else          
        fprintf(file, "%d %d %lf+%lfi\n", Ai[i] + 1, j + 1, SUPERLU_SCALAR(Ax[i]));
#endif          
      fprintf(file, "];\n%s = spconvert(temp);\n", var_name);
      
      return true;
      
    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DX\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(&nnz, sizeof(int), 1, file);
      hermes_fwrite(Ap, sizeof(int), size + 1, file);
      hermes_fwrite(Ai, sizeof(int), nnz, file);
      hermes_fwrite(Ax, sizeof(slu_scalar), nnz, file);
      return true;
    }
    
    default:
      return false;
  }
}

int SuperLUMatrix::get_matrix_size() const
{
  return size;
}

/* THIS WAS WRONG
int SuperLUMatrix::get_matrix_size() const
{
  _F_
  //           Ax               Ai                 Ap                      nnz
  return (sizeof(scalar) + sizeof(int)) * nnz + sizeof(int)*(size+1) + sizeof(int);
}
*/

int SuperLUMatrix::get_nnz() const
{
  return nnz;
}

double SuperLUMatrix::get_fill_in() const
{
  _F_
  return nnz / (double) (size * size);
}

// SuperLUVector /////////////////////////////////////////////////////////////////////////////////////

SuperLUVector::SuperLUVector()
{
  _F_
  v = NULL;
  size = 0;
}

SuperLUVector::~SuperLUVector()
{
  _F_
  this->free();
}

void SuperLUVector::alloc(int n)
{
  _F_
  this->free();
  size = n;
  v = new slu_scalar[n];
  zero();
}

void SuperLUVector::zero()
{
  _F_
  memset(v, 0, size * sizeof(slu_scalar));
}

void SuperLUVector::free()
{
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void SuperLUVector::set(int idx, scalar y)
{
  _F_
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] = y;
#else
  if (idx >= 0) {
    v[idx].r = y.real();
    v[idx].i = y.imag();
  }
#endif
}

void SuperLUVector::add(int idx, scalar y)
{
  _F_
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] += y;
#else
  if (idx >= 0) {
    v[idx].r += y.real();
    v[idx].i += y.imag();
  }
#endif
}

void SuperLUVector::add(int n, int *idx, scalar *y)
{
  _F_
  for (int i = 0; i < n; i++)
    if (idx[i] >= 0) {
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
      v[idx[i]] += y[i];
#else
      v[idx[i]].r += y[i].real();
      v[idx[i]].i += y[i].imag();
#endif
    }
}

bool SuperLUVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  switch (fmt) 
  {
    case DF_NATIVE:
    case DF_PLAIN_ASCII:
      for (int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", SUPERLU_SCALAR(v[i]));
      
      return true;
      
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", SUPERLU_SCALAR(v[i]));
      fprintf(file, " ];\n");
      return true;
      
    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DR\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(v, sizeof(scalar), size, file);
      return true;
    }
    
    default:
      return false;
  }
}

// SUPERLU solver ////////////////////////////////////////////////////////////////////////////////////

bool SuperLUSolver::check_status(int info)
{
  _F_
  if (info < 0)
  {
    warning("SuperLU: %d-th argument had an illegal value.", -info);
    return false;
  }
  else if (info == 0)
  {
    // Success.
    return true;
  }
  else if (info > 0 && info <= m->size)
  {
    warning("SuperLU: Factor U is singular, solution could not be computed.");
    return false;
  }
  else if (info == m->size + 1)
  {
    warning("SuperLU: RCOND is less than machine precision "
            "(system matrix is singular to working precision).");
    return true;
  }
  else if (info > m->size + 1)
  {
    warning("SuperLU: Not enough memory.\n Failure when %.3f MB were allocated.",
            (info - m->size)/1e6);
    return false;
  }
  
  return false;
}
  
SuperLUSolver::SuperLUSolver(SuperLUMatrix *m, SuperLUVector *rhs) 
  : LinearSolver(HERMES_FACTORIZE_FROM_SCRATCH), m(m), rhs(rhs), 
      local_Ai(NULL), local_Ap(NULL), local_Ax(NULL), local_rhs(NULL)
{
  _F_
#ifdef WITH_SUPERLU
  R = NULL;
  C = NULL;
  perm_r = NULL;
  perm_c = NULL;
  etree = NULL; 
#ifndef SLU_MT
  *equed = '\0';
#endif  

  // Set the default input options:
#ifdef SLU_MT
  // I am not sure if this will work well on Windows:
  // http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c
  char *nt_var = getenv("OMP_NUM_THREADS");
  if (nt_var)
    options.nprocs          = std::max(1, atoi(nt_var));
  else
    options.nprocs          = 1;
  
  options.fact              = EQUILIBRATE;  // Rescale the matrix if neccessary.
  options.trans             = NOTRANS;      // Not solving the transposed problem.
  options.refact            = NO;           // Factorize from scratch for the first time.
  options.diag_pivot_thresh = 1.0;          // Use partial pivoting during GEM.
  options.usepr             = NO;           // Let SuperLU compute the row permutations.
  options.drop_tol          = 0.0;          // Not yet implemented in SuperLU_MT 2.0.
  options.SymmetricMode     = NO;           // Assume general non-symmetric problem.
    
  // Default options related to the supernodal algorithm.
  options.panel_size        = sp_ienv(1);
  options.relax             = sp_ienv(2);
#else
  /*
  options.Fact = DOFACT;
  options.Equil = YES;
  options.ColPerm = COLAMD;
  options.DiagPivotThresh = 1.0;
  options.Trans = NOTRANS;
  options.IterRefine = NOREFINE;
  options.SymmetricMode = NO;
  options.PivotGrowth = NO;
  options.ConditionNumber = NO;
  options.PrintStat = YES;
  */
  set_default_options(&options);  // This function is only present in the sequential SLU.
#endif  

  options.PrintStat = YES;   // Set to NO to suppress output.
  
  has_A = has_B = inited = false;
#else
  error(SUPERLU_NOT_COMPILED);
#endif
}

SuperLUSolver::~SuperLUSolver()
{
  _F_
  free_factorization_data();
  free_matrix();
  free_rhs();
  
  if (local_Ai)  delete [] local_Ai;
  if (local_Ap)  delete [] local_Ap;
  if (local_Ax)  delete [] local_Ax;
  if (local_rhs) delete [] local_rhs;
}

bool SuperLUSolver::solve()
{
  _F_
#ifdef WITH_SUPERLU
  assert(m != NULL);
  assert(rhs != NULL);
  
  TimePeriod tmr;
  
  // Initialize the statistics variable.
  slu_stat_t stat;
  SLU_INIT_STAT(&stat);
  
  // Prepare data structures serving as input for the solver driver 
  // (according to the chosen factorization reuse strategy).
  void *work = NULL;        // Explicit pointer to the factorization workspace 
                            // (unused, see below).
  int lwork = 0;            // Space for the factorization will be allocated 
                            // internally by system malloc.
  double ferr = 1.0;        // Estimated relative forward error 
                            // (unused unless iterative refinement is performed).
  double berr = 1.0;        // Estimated relative backward error 
                            // (unused unless iterative refinement is performed).
  slu_memusage_t memusage;  // Record the memory usage statistics.
  double rpivot_growth;     // The reciprocal pivot growth factor.
  double rcond;             // The estimate of the reciprocal condition number.                          
#ifdef SLU_MT
  options.work = work;
  options.lwork = lwork;
#endif

  if ( !setup_factorization() )
  {
    warning("LU factorization could not be completed.");
    return false;
  }
  
  // If the previous factorization of A is to be fully reused as an input for the solver driver,
  // keep the (possibly rescaled) matrix from the last factorization, otherwise recreate it 
  // from the master SuperLUMatrix pointed to by this->m (this also applies to the case when 
  // A does not yet exist).
  if (!has_A || factorization_scheme != HERMES_REUSE_FACTORIZATION_COMPLETELY)
  {
    if (A_changed) 
      free_matrix();
    
    if (!has_A)
    {
      // A will be created from the local copy of the value and index arrays, because these
      // may be modified by the solver driver.
      if (local_Ai) delete [] local_Ai;
      local_Ai = new int [m->nnz];
      memcpy(local_Ai, m->Ai, m->nnz * sizeof(int));
      
      if (local_Ap) delete [] local_Ap;
      local_Ap = new int [m->size+1];
      memcpy(local_Ap, m->Ap, (m->size+1) * sizeof(int));
      
      if (local_Ax) delete [] local_Ax;
      local_Ax = new slu_scalar [m->nnz];
      memcpy(local_Ax, m->Ax, m->nnz * sizeof(slu_scalar));
      
      // Create new general (non-symmetric), column-major, non-supernodal, size X size matrix.
      SLU_CREATE_CSC_MATRIX(&A, m->size, m->size, m->nnz, local_Ax, local_Ai, local_Ap, SLU_NC, SLU_DTYPE, SLU_GE);
      
      has_A = true;
    }
  }

  // Recreate the input rhs for the solver driver from a local copy of the new value array.
  free_rhs();
 
  if (local_rhs) delete [] local_rhs;
  local_rhs = new slu_scalar [rhs->size];
  memcpy(local_rhs, rhs->v, rhs->size * sizeof(slu_scalar));
  
  SLU_CREATE_DENSE_MATRIX(&B, rhs->size, 1, local_rhs, rhs->size, SLU_DN, SLU_DTYPE, SLU_GE);
  
  has_B = true;
  
  // Initialize the solution variable.
  SuperMatrix X;
  slu_scalar *x;
  if ( !(x = SLU_SCALAR_MALLOC(m->size)) ) 
    error("Malloc fails for x[].");
  SLU_CREATE_DENSE_MATRIX(&X, m->size, 1, x, m->size, SLU_DN, SLU_DTYPE, SLU_GE);
    
  // Solve the system.
  int info;

#ifdef SLU_MT  
  if (options.refact == NO)
  {
    // Get column permutation vector perm_c[], according to the first argument:
    //  0: natural ordering 
    //  1: minimum degree ordering on structure of A'*A
    //  2: minimum degree ordering on structure of A'+A
    //  3: approximate minimum degree for unsymmetric matrices   
    get_perm_c(1, &A, perm_c);
  }
   
/*
  // Compute reciprocal pivot growth, estimate reciprocal condition number of A, solve,
  // perform iterative refinement of the solution and estimate forward and backward error.
  // Memory usage will be acquired at the end. If A is singular, info will be set to A->ncol+1.
  //
  slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
                        &L, &U, &B, &X, &rpivot_growth, &rcond, &ferr, &berr, 
                        &stat, &memusage, &info );
*/
                        
  // ... OR ...
  
  // Estimate reciprocal condition number of A and solve the system. If A is singular, info
  // will be set to A->ncol+1.
  //
  slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
                        &L, &U, &B, &X, NULL, &rcond, NULL, NULL, 
                        &stat, NULL, &info );

  // ... OR ...

/*  
  // Do not check the regularity of A and just solve the system.
  //
  slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
                        &L, &U, &B, &X, NULL, NULL, NULL, NULL, 
                        &stat, NULL, &info );                        
*/
#else
  SLU_SOLVER_DRIVER(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U,
                    work, lwork, &B, &X, &rpivot_growth, &rcond, &ferr, &berr,
                    &memusage, &stat, &info);
#endif
                    
  // A and B may have been multiplied by the scaling vectors R and C on the output of the 
  // solver. If the next call to the solver should reuse factorization only partially,
  // it will need the original unscaled matrix - this will indicate such situation 
  // (rhs is always recreated anew).
#ifdef SLU_MT  
  A_changed = (equed != NOEQUIL);
#else
  A_changed = (*equed != 'N');
#endif  

  bool factorized = check_status(info);
  
  if (factorized) 
  {
    delete [] sln;
    sln = new scalar[m->size];
    
    slu_scalar *sol = (slu_scalar*) ((DNformat*) X.Store)->nzval; 
    
    for (int i = 0; i < rhs->size; i++)
#if !defined(H1D_COMPLEX) && !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)      
      sln[i] = sol[i];
#else
      sln[i] = cplx(sol[i].r, sol[i].i);
#endif
  }
  
  // If required, print statistics.
  if ( options.PrintStat ) SLU_PRINT_STAT(&stat);
  
  // Free temporary local variables.
  StatFree(&stat);
  SUPERLU_FREE (x);
  Destroy_SuperMatrix_Store(&X);
  
  tmr.tick();
  time = tmr.accumulated();
  
  return factorized;
#else
  return false;
#endif
}

bool SuperLUSolver::setup_factorization()
{
  _F_
#ifdef WITH_SUPERLU
  if (has_A && factorization_scheme != HERMES_FACTORIZE_FROM_SCRATCH && A.nrow != m->size)
  {
    warning("You cannot reuse factorization structures for factorizing matrices of different sizes.");
    return false;
  }
  
  // Always factorize from scratch for the first time.
  int eff_fact_scheme;
  if (!inited)
    eff_fact_scheme = HERMES_FACTORIZE_FROM_SCRATCH;
  else
    eff_fact_scheme = factorization_scheme;
  
  // Prepare factorization structures. In case of a particular reuse scheme, comments are given
  // to clarify which arguments will be reused and which will be reset by the dgssvx (zgssvx) routine. 
  // It was determined empirically by running the dlinsolx2 example from SuperLU, setting options.Fact
  // to the appropriate value and reallocating the various structures before the second run of dgssvx,
  // and observing when segfault will happen and when not. It is actually not needed to reallocate 
  // the structures by hand, but comments at various places of the SuperLU 4.0 library contradict 
  // each other and often lead to segfault when the structures are reallocated according to them.
  // It might thus bring some insight into how SuperLU works and how to correctly use it 
  // (the PDF documentation is, unfortunately, even less helpful).
  switch (eff_fact_scheme)
  {
    case HERMES_FACTORIZE_FROM_SCRATCH:
      // This case should generally allow for solving a completely new system, i.e. for a change of 
      // matrix and rhs size - for simplicity, we reallocate the structures every time.
      
      // Clear the structures emanating from previous factorization.
      free_factorization_data();
      
      // Allocate the row/column reordering vectors.
      if ( !(perm_c = intMalloc(m->size)) ) 
        error("Malloc fails for perm_c[].");
      if ( !(perm_r = intMalloc(m->size)) ) 
        error("Malloc fails for perm_r[].");
      
      // Allocate vectors with row/column scaling factors.
      if ( !(R = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) ) 
        error("SUPERLU_MALLOC fails for R[].");
      if ( !(C = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) )
        error("SUPERLU_MALLOC fails for C[].");

#ifdef SLU_MT
      options.fact = EQUILIBRATE;
      options.refact = NO;      
      options.perm_c = perm_c;
      options.perm_r = perm_r;
#else 
      // Allocate additional structures used by the driver routine of sequential SuperLU.
      // Elimination tree is contained in the options structure in SuperLU_MT.
      if ( !(etree = intMalloc(m->size)) )    
        error("Malloc fails for etree[].");

      options.Fact = DOFACT;
#endif      
      A_changed = true;
      break;
    case HERMES_REUSE_MATRIX_REORDERING:
      // needed from previous:      etree, perm_c
      // not needed from previous:  perm_r, R, C, L, U, equed     
#ifdef SLU_MT
      options.fact = EQUILIBRATE;
      options.refact = YES;
#else
      options.Fact = SamePattern;
#endif      
      // L,U matrices may be reused without reallocating.
      // SLU_DESTROY_L(&L);
      // SLU_DESTROY_U(&U);
      break;
    case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
      // needed from previous:      etree, perm_c, perm_r, L, U
      // not needed from previous:  R, C, equed
#ifdef SLU_MT
      // MT version of SLU cannot reuse the equilibration factors (R, C), so
      // this is the same as the previous case.
      options.fact = EQUILIBRATE; 
      options.refact = YES;
#else
      options.Fact = SamePattern_SameRowPerm;
#endif      
      break;
    case HERMES_REUSE_FACTORIZATION_COMPLETELY:
      // needed from previous:      perm_c, perm_r, equed, L, U
      // not needed from previous:  etree, R, C
#ifdef SLU_MT
      options.fact = FACTORED;
      options.refact = YES;
#else      
      options.Fact = FACTORED;
#endif      
      break;
  }
  
  inited = true;
  
  return true;
#else
  return false;
#endif
}

void SuperLUSolver::free_matrix()
{
  _F_
#ifdef WITH_SUPERLU  
  if (has_A)
  {
    Destroy_SuperMatrix_Store(&A);
    has_A = false;
  }
#endif  
}

void SuperLUSolver::free_rhs()
{
  _F_
  #ifdef WITH_SUPERLU  
  if (has_B)
  {
    Destroy_SuperMatrix_Store(&B);
    has_B = false;
  }
  #endif  
}

void SuperLUSolver::free_factorization_data()
{ 
  _F_
#ifdef WITH_SUPERLU
  if (inited)
  {
#ifdef SLU_MT    
    SUPERLU_FREE(options.etree);
    SUPERLU_FREE(options.colcnt_h);
    SUPERLU_FREE(options.part_super_h);
    Destroy_CompCol_Permuted(&AC);
#else
    SUPERLU_FREE (etree);
#endif    
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SLU_DESTROY_L(&L);
    SLU_DESTROY_U(&U);
    inited = false;
  }
#endif
}

#ifdef SLU_MT
// This is a modification of the original p*gssvx routines from the SuperLU_MT library.
//
// The original routines have been changed in view of our applications, i.e. 
//  * only one right hand side is allowed, 
//  * some initial parameter checks have been omitted, 
//  * macros allowing abstraction from the fundamental scalar datatype have been used
//  * some phases of the calculation may be omitted for speed-up (less information about
//    the matrix/solution can then be acquired, however),
//  * deallocation at the end of the routine has been removed (this was neccessary to 
//    enable factorization reuse).
//
// See the correspondingly named attributes of SuperLUSolver class for brief description 
// of most parameters or the library source code for pdgssvx for more details. You may pass
// NULL for
//  * recip_pivot_growth  - reciprocal pivot growth factor will then not be computed;
//                          reip_pivot_growth much less than one may indicate poor 
//                          stability of the factorization;
//  * rcond               - estimate of the reciprocal condition number of matrix A will
//                          then not be computed; this will prevent detection of singularity
//                          of matrix A;
//  * ferr or berr        - iterative refinement of the solution will then not be performed;
//                          this also prevents computation of forward and backward error 
//                          estimates of the computed solution;
//  * memusage            - memory usage during the factorization/solution will not be queried.
//
void slu_mt_solver_driver(slu_options_t *options, SuperMatrix *A, 
                          int *perm_c, int *perm_r, SuperMatrix *AC,
                          equed_t *equed, double *R, double *C,
                          SuperMatrix *L, SuperMatrix *U,
                          SuperMatrix *B, SuperMatrix *X, 
                          double *recip_pivot_growth, double *rcond, 
                          double *ferr, double *berr, 
                          slu_stat_t *stat, slu_memusage_t *memusage,
                          int *info)
{
  /* Profiling variables. */
  double    t0;
  flops_t   flopcnt;
  
  /* State variables. */
  int dofact = (options->fact == DOFACT);
  int equil = (options->fact == EQUILIBRATE);
  int notran = (options->trans == NOTRANS);
  int colequ, rowequ;
  
  /* Right hand side and solution vectors. */
  DNformat *Bstore = (DNformat*) B->Store;
  DNformat *Xstore = (DNformat*) X->Store;
  slu_scalar *Bmat = (slu_scalar*) Bstore->nzval;
  slu_scalar *Xmat = (slu_scalar*) Xstore->nzval;
    
  *info = 0;
  
  /* ------------------------------------------------------------
  Diagonal scaling to equilibrate the matrix.
  ------------------------------------------------------------*/
  if (dofact || equil) 
  {
    *equed = NOEQUIL;
    rowequ = colequ = FALSE;
  } 
  else 
  {
    rowequ = (*equed == ROW) || (*equed == BOTH);
    colequ = (*equed == COL) || (*equed == BOTH);
  }
  
  if ( equil ) 
  {
    t0 = SuperLU_timer_();
    /* Compute row and column scalings to equilibrate the matrix A. */
    int info1;
    double rowcnd, colcnd, amax;
    SLU_GSEQU(A, R, C, &rowcnd, &colcnd, &amax, &info1);
    
    if ( info1 == 0 ) {
      /* Equilibrate matrix A. */
      SLU_LAQGS(A, R, C, rowcnd, colcnd, amax, equed);
      rowequ = (*equed == ROW) || (*equed == BOTH);
      colequ = (*equed == COL) || (*equed == BOTH);
    }
    stat->utime[EQUIL] = SuperLU_timer_() - t0;
  }
  
  /* ------------------------------------------------------------
  Scale the right hand side.
  ------------------------------------------------------------*/
  if ( notran ) 
  {
    if ( rowequ ) 
      for (int i = 0; i < A->nrow; ++i) 
        SLU_MULT(Bmat[i], R[i]);
  } 
  else if ( colequ ) 
  {
    for (int i = 0; i < A->nrow; ++i)
      SLU_MULT(Bmat[i], C[i]);
  }
  
  /* ------------------------------------------------------------
  Perform the LU factorization.
  ------------------------------------------------------------*/
  if ( dofact || equil ) 
  {  
    /* Obtain column etree, the column count (colcnt_h) and supernode
    partition (part_super_h) for the Householder matrix. */
    if (options->refact == NO)
    {
      t0 = SuperLU_timer_();
      SLU_SP_COLORDER(A, perm_c, options, AC);
      stat->utime[ETREE] = SuperLU_timer_() - t0;
    }
     
    /* Compute the LU factorization of A*Pc. */
    t0 = SuperLU_timer_();
    SLU_GSTRF(options, AC, perm_r, L, U, stat, info);
    stat->utime[FACT] = SuperLU_timer_() - t0;
    
    flopcnt = 0;
    for (int i = 0; i < options->nprocs; ++i) flopcnt += stat->procstat[i].fcops;
    stat->ops[FACT] = flopcnt;
    
    if ( options->lwork == -1 ) 
    {
      if (memusage)
        memusage->total_needed = *info - A->ncol;
      return;
    }
  }
  
  if ( *info > 0 ) 
  {
    if ( *info <= A->ncol ) 
    {
      /* Compute the reciprocal pivot growth factor of the leading
        rank-deficient *info columns of A. */
      if (recip_pivot_growth)
      *recip_pivot_growth = SLU_PIVOT_GROWTH(*info, A, perm_c, L, U);
    }
  } 
  else 
  {
    /* ------------------------------------------------------------
      Compute the reciprocal pivot growth factor *recip_pivot_growth.
      ------------------------------------------------------------*/
    if (recip_pivot_growth)
      *recip_pivot_growth = SLU_PIVOT_GROWTH(A->ncol, A, perm_c, L, U);

    /* ------------------------------------------------------------
      Estimate the reciprocal of the condition number of A.
      ------------------------------------------------------------*/
    if (rcond) 
    {
      t0 = SuperLU_timer_();
      
      // Next two lines are a bit complicated, but taken as they appear
      // in the original library function.
      char norm[1];
      *(unsigned char *)norm = (notran) ? '1' : 'I';
      
      double anorm = SLU_LANGS(norm, A);
      SLU_GSCON(norm, L, U, anorm, rcond, info);
      stat->utime[RCOND] = SuperLU_timer_() - t0;
    }  

    /* ------------------------------------------------------------
      Compute the solution matrix X.
      ------------------------------------------------------------*/
    // Save a copy of the right hand side.
    memcpy(Xmat, Bmat, B->nrow * sizeof(slu_scalar)); 
            
    t0 = SuperLU_timer_();
    SLU_GSTRS(options->trans, L, U, perm_r, perm_c, X, stat, info);
    stat->utime[SOLVE] = SuperLU_timer_() - t0;
    stat->ops[SOLVE] = stat->ops[TRISOLVE];
      
    /* ------------------------------------------------------------
      Use iterative refinement to improve the computed solution and
      compute error bounds and backward error estimates for it.
      ------------------------------------------------------------*/
    if (ferr && berr)
    {
      t0 = SuperLU_timer_();
      SLU_GSRFS(options->trans, A, L, U, perm_r, perm_c, *equed,
                R, C, B, X, ferr, berr, stat, info);
      stat->utime[REFINE] = SuperLU_timer_() - t0;
    }

    /* ------------------------------------------------------------
      Transform the solution matrix X to a solution of the original
      system.
      ------------------------------------------------------------*/
    if ( notran ) 
    {
      if ( colequ ) 
        for (int i = 0; i < A->nrow; ++i)
          SLU_MULT(Xmat[i], C[i]);
    } 
    else if ( rowequ ) 
    {
      for (int i = 0; i < A->nrow; ++i)
        SLU_MULT(Xmat[i], R[i]);
    }

    /* Set INFO = A->ncol+1 if the matrix is singular to 
      working precision.*/
    char param[1]; param[0] = 'E';
    if ( rcond && *rcond < SLU_LAMCH_(param) ) *info = A->ncol + 1; 
  }

  if (memusage)
    SLU_QUERY_SPACE(options->nprocs, L, U, options->panel_size, memusage);
}
#endif
