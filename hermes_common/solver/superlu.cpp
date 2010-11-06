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
  _F_
  /*           Ax               Ai                 Ap                      nnz     */
  return (sizeof(scalar) + sizeof(int)) * nnz + sizeof(int)*(size+1) + sizeof(int);
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
  R = C = NULL;
  perm_r = perm_c = etree = NULL;
  *equed = '\0';
  
  /* Set the default input options:
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
  set_default_options(&options);
  options.PrintStat = YES;   // Set to NO to suppress output.
  
  has_A = has_B = inited = false;
#else
  error(SUPERLU_NOT_COMPILED);
#endif
}

SuperLUSolver::~SuperLUSolver()
{
  _F_
  free_factorization_structures();
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
  
  // Initialize SuperLU.
  void *work = NULL;      // Explicit pointer to the factorization workspace 
                          // (unused, see below).
  int lwork = 0;          // Space for the factorization will be allocated 
                          // internally by system malloc.
  double ferr = 1.0;      // Estimated relative forward error 
                          // (unused unless iterative refinement is performed).
  double berr = 1.0;      // Estimated relative backward error 
                          // (unused unless iterative refinement is performed).
  mem_usage_t mem_usage;  // Record the memory usage statistics.
  double rpivot_growth;   // The reciprocal pivot growth factor.
  double rcond;           // The estimate of the reciprocal condition number.

  int info;
    
  // Prepare data structures serving as input for the solver driver 
  // (according to the chosen factorization reuse strategy).
  if ( !prepare_factorization_structures() )
  {
    warning("LU factorization could not be completed.");
    return false;
  }
  
  // If the previous factorization of A is to be fully reused as an input for the solver driver,
  // keep the (possibly rescaled) matrix from the last factorization, otherwise recreate it 
  // from the master SuperLUMatrix pointed to by this->m (this also applies to the case when 
  // A does not yet exist).
  if (factorization_scheme != HERMES_REUSE_FACTORIZATION_COMPLETELY)
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
    ABORT("Malloc fails for x[].");
  SLU_CREATE_DENSE_MATRIX(&X, m->size, 1, x, m->size, SLU_DN, SLU_DTYPE, SLU_GE);
  
  // Initialize the statistics variable.
  SuperLUStat_t stat;
  StatInit(&stat);
  
  // Solve the system.
  SLU_SOLVER_DRIVER(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U,
                    work, lwork, &B, &X, &rpivot_growth, &rcond, &ferr, &berr,
                    &mem_usage, &stat, &info);
  
  // A and B may have been multiplied by the scaling vectors R and C on the output of the 
  // solver. If the next call to the solver should reuse factorization only partially,
  // it will need the original unscaled matrix - this will indicate such situation 
  // (rhs is always recreated anew).
  A_changed = (options.Equil == YES && *equed != 'N');
     
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
  if ( options.PrintStat ) StatPrint(&stat);
  
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

bool SuperLUSolver::prepare_factorization_structures()
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
  switch (factorization_scheme)
  {
    case HERMES_FACTORIZE_FROM_SCRATCH:
      // This case should generally allow for solving a completely new system, i.e. for a change of 
      // matrix and rhs size - for simplicity, we reallocate the structures every time.
      
      // Clear the structures emanating from previous factorization.
      free_factorization_structures();
      
      // Allocate the new structures (internal arrays of L, U are allocated automatically by SuperLU).
      if ( !(etree = intMalloc(m->size)) )  
        ABORT("Malloc fails for etree[].");
      if ( !(perm_c = intMalloc(m->size)) ) 
        ABORT("Malloc fails for perm_c[].");
      if ( !(perm_r = intMalloc(m->size)) ) 
        ABORT("Malloc fails for perm_r[].");
      if ( !(R = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
      if ( !(C = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
            
      options.Fact = DOFACT;
      A_changed = true;
      break;
    case HERMES_REUSE_MATRIX_REORDERING:
      // needed from previous:      etree, perm_c
      // not needed from previous:  perm_r, R, C, L, U, equed     
      options.Fact = SamePattern;
      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
      break;
    case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
      // needed from previous:      etree, perm_c, perm_r, L, U
      // not needed from previous:  R, C, equed
      options.Fact = SamePattern_SameRowPerm;
      break;
    case HERMES_REUSE_FACTORIZATION_COMPLETELY:
      // needed from previous:      perm_c, perm_r, equed, L, U
      // not needed from previous:  etree, R, C
      options.Fact = FACTORED;
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

void SuperLUSolver::free_factorization_structures()
{ 
  _F_
#ifdef WITH_SUPERLU
  if (inited)
  {
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    inited = false;
  }
#endif
}
