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

#include "pardiso.h"

#include "../trace.h"
#include "../error.h"
#include "../utils.h"
#include "../callstack.h"

#ifdef WITH_PARDISO
  #ifdef __cplusplus
    extern "C" {
  #endif
      extern int pardisoinit_(void *, int *, int *, int*, double *, int *);
      extern int pardiso_(void *, int *, int *, int *, int *, int *,
                          scalar *, int *, int *, int *, int *, int *, int *, 
                          scalar *, scalar *, int *, double *);
  #ifdef __cplusplus
    }
  #endif
  
  #define PARDISOINIT pardisoinit_
  #define PARDISO pardiso_
#else
  #define PARDISOINIT
  #define PARDISO
#endif

// Binary search for the location of a particular CSC/CSR matrix entry.
//
// Typically, we search for the index into Ax that corresponds to a given 
// row (CSC) or column (CSR) ('idx') among indices of nonzero values in 
// a particular column (CSC) or row (CSR) ('Ai').
//
static int find_position(int *Ai, int Alen, int idx) {
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

PardisoMatrix::PardisoMatrix() {
  _F_
  size = 0; nnz = 0;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;
}

PardisoMatrix::~PardisoMatrix() {
  _F_
  free();
}

void PardisoMatrix::pre_add_ij(int row, int col) {
  _F_
  SparseMatrix::pre_add_ij(col, row);
}

void PardisoMatrix::alloc() {
  _F_
  assert(pages != NULL);
  assert(size > 0);

  // initialize the arrays Ap and Ai
  Ap = new int[size + 1];
  MEM_CHECK(Ap);
  int aisize = get_num_indices();
  Ai = new int[aisize];
  MEM_CHECK(Ai);

  // sort the indices and remove duplicities, insert into Ai
  int i, pos = 0;
  for (i = 0; i < size; i++) {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;

  delete[] pages;
  pages = NULL;

  nnz = Ap[size];
  
  Ax = new scalar[nnz];
  MEM_CHECK(Ax);
  zero();
}

void PardisoMatrix::free() {
  _F_
  nnz = 0;
  delete [] Ap; Ap = NULL;
  delete [] Ai; Ai = NULL;
  delete [] Ax; Ax = NULL;
}

scalar PardisoMatrix::get(int m, int n)
{
  _F_
  // Find n-th column in the m-th row.
  int mid = find_position(Ai + Ap[m], Ap[m + 1] - Ap[m], n);

  if (mid < 0) // if the entry has not been found
    return 0.0;   
  else 
    return Ax[Ap[m]+mid];
}

void PardisoMatrix::zero() {
  _F_
    memset(Ax, 0, sizeof(scalar) * nnz);
}

void PardisoMatrix::add(int m, int n, scalar v) {
  _F_
  if (v != 0.0 && m >= 0 && n >= 0) // ignore dirichlet DOFs
  {   
    // Find n-th column in the m-th row.
    int pos = find_position(Ai + Ap[m], Ap[m + 1] - Ap[m], n);
    // Make sure we are adding to an existing non-zero entry.
    if (pos < 0) 
      error("Sparse matrix entry not found");
    
    Ax[Ap[m]+pos] += v;
  }
}

void PardisoMatrix::add(int m, int n, scalar **mat, int *rows, int *cols) {
  _F_
  for (int i = 0; i < m; i++)       // rows
    for (int j = 0; j < n; j++)     // cols
      add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool PardisoMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
  _F_
  switch (fmt) 
  {
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, nnz, nnz);
      for (int j = 0; j < size; j++)
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
          fprintf(file, "%d %d " SCALAR_FMT ";\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
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
      hermes_fwrite(Ax, sizeof(scalar), nnz, file);
      return true;
    }

    case DF_PLAIN_ASCII:
      EXIT(HERMES_ERR_NOT_IMPLEMENTED);
      return false;

    default:
      return false;
  }
}

int PardisoMatrix::get_matrix_size() const {
  _F_
  assert(Ap != NULL);
  /*          Ai             Ax                     Ap                    nnz       */    
  return (sizeof(int) + sizeof(scalar)) * nnz + sizeof(int)*(size+1) + sizeof(int);
}

double PardisoMatrix::get_fill_in() const {
  _F_
  return nnz / (double) (size * size);
}


// PardisoVector ///////

PardisoVector::PardisoVector() {
  _F_
  v = NULL;
  size = 0;
}

PardisoVector::~PardisoVector() {
  _F_
  free();
}

void PardisoVector::alloc(int n) {
  _F_
  free();
  v = new scalar[n];
  MEM_CHECK(v);
  size = n;
  zero();
}

void PardisoVector::zero() {
  _F_
  memset(v, 0, size * sizeof(scalar));
}

void PardisoVector::free() {
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void PardisoVector::set(int idx, scalar y) {
  _F_
  if (idx >= 0) v[idx] = y;
}

void PardisoVector::add(int idx, scalar y) {
  _F_
  if (idx >= 0) v[idx] += y;
}

void PardisoVector::extract(scalar *v) const
{
  return;
}

void PardisoVector::add(int n, int *idx, scalar *y) {
  _F_
  for (int i = 0; i < n; i++)
    if (idx[i] >= 0) v[idx[i]] += y[i];
}

bool PardisoVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt) {
  _F_
  switch (fmt) 
  {
    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (int i = 0; i < this->size; i++)
        fprintf(file, SCALAR_FMT ";\n", SCALAR(v[i]));
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

    case DF_PLAIN_ASCII:
      EXIT(HERMES_ERR_NOT_IMPLEMENTED);
      return false;

    default:
      return false;
  }
}

// PARDISO solver //////

PardisoLinearSolver::PardisoLinearSolver(PardisoMatrix *m, PardisoVector *rhs)
  : LinearSolver(), m(m), rhs(rhs)
{
  _F_
#ifdef WITH_PARDISO
#else
  error(PARDISO_NOT_COMPILED);
#endif
}

PardisoLinearSolver::~PardisoLinearSolver() {
  _F_
#ifdef WITH_PARDISO  
  //if (m != NULL) delete m;
  //if (rhs != NULL) delete rhs;
#endif
}

bool PardisoLinearSolver::solve() {
  _F_
#ifdef WITH_PARDISO
  assert(m != NULL);
  assert(rhs != NULL);

  bool res = true;
  int n = m->size;

  try {
    // Numbers of processors, value of OMP_NUM_THREADS
    int num_procs;
    char *var = getenv("OMP_NUM_THREADS");
    if (var != NULL) sscanf(var, "%d", &num_procs);
    else num_procs = 1;

#if defined (H2D_COMPLEX) || defined (H3D_COMPLEX)
    int mtype = 13;		// Complex unsymmetric matrix
#else    
    int mtype = 11;   // Real unsymmetric matrix
#endif

    int nrhs = 1;		// Number of right hand sides
    int nnz = m->Ap[n];	// The number of nonzero elements

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *pt[64];
    // Pardiso control parameters. Consult Pardiso manual for interpretation.
    int iparm[64];
    double dparm[64];
    int maxfct, mnum, phase, err, msglvl;
    // Auxiliary variables.
    scalar ddum; // Double dummy
    int idum; // Integer dummy.

    iparm[2] = num_procs; // Number of SMP threads. Remaining entries in iparm
                          // will be filled by default values by PARDISOINIT.
                          
    int solver = 0;   // Sparse direct solver. Set solver = 1 for multi-recursive
                      // iterative solver.

    maxfct = 1;		// Maximum number of numerical factorizations.
    mnum = 1;		// Which factorization to use.
    msglvl = 0;		// Do not print statistical information
    err = 0;		// Initialize error flag

    // Convert matrix from 0-based C-notation to Fortran 1-based notation.
    for (int i = 0; i < n + 1; i++) m->Ap[i] += 1;
    for (int i = 0; i < nnz; i++) m->Ai[i] += 1;

    TimePeriod tmr;
    
    // Setup Pardiso control parameters.
    PARDISOINIT(pt, &mtype, &solver, iparm, dparm, &err);

    // .. Reordering and Symbolic Factorization. This step also allocates
    // all memory that is necessary for the factorization.
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err, dparm);
    if (err != 0) {
      // ERROR during symbolic factorization: err
      throw ERR_FAILURE;
    }

    // .. Numerical factorization.
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err, dparm);
    if (err != 0) {
      // ERROR during numerical factorization: err
      throw ERR_FAILURE;
    }

    // .. Back substitution and iterative refinement.
    delete [] sln;
    sln = new scalar[m->size];
    MEM_CHECK(sln);
    memset(sln, 0, (m->size) * sizeof(scalar));

    phase = 33;
    iparm[7] = 1; // Max numbers of iterative refinement steps.
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, m->Ax, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, rhs->v, sln, &err, dparm);
    if (err != 0) {
      // ERROR during solution: err
      throw ERR_FAILURE;
    }

    tmr.tick();
    time = tmr.accumulated();

    //  Convert matrix back to 0-based C-notation.
    for (int i = 0; i < n + 1; i++) m->Ap[i] -= 1;
    for (int i = 0; i < nnz; i++) m->Ai[i] -= 1;

    // .. Termination and release of memory.
    phase = -1; // Release internal memory.
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, m->Ap, m->Ai, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &err, dparm);
  }
  catch (int e) {
    error = e;
    res = false;
  }

  return res;
#else
  return false;
#endif
}
