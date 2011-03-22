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

#include "mumps.h"
#include "../trace.h"
#include "../error.h"
#include "../utils.h"
#include "../callstack.h"

#ifndef HERMES_COMMON_COMPLEX
  #define MUMPS         dmumps_c
#else
  #define MUMPS         zmumps_c
#endif

#define USE_COMM_WORLD  -987654

#ifdef WITH_MUMPS

extern "C" {
  extern void MUMPS(MUMPS_STRUCT *mumps_param_ptr);
}

#else
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

MumpsMatrix::MumpsMatrix()
{
  _F_
  nnz = 0;
  size = 0;
  irn = NULL;
  jcn = NULL;
  Ax = NULL;
  Ap = NULL;
  Ai = NULL;
}

MumpsMatrix::~MumpsMatrix()
{
  _F_
  free();
}

void MumpsMatrix::alloc()
{
  _F_
  assert(pages != NULL);

  // initialize the arrays Ap and Ai
  Ap = new unsigned int [size + 1];
  MEM_CHECK(Ap);
  int aisize = get_num_indices();
  Ai = new int [aisize];
  MEM_CHECK(Ai);

  // sort the indices and remove duplicities, insert into Ai
  unsigned int i, pos = 0;
  for (i = 0; i < size; i++) {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;

  delete [] pages;
  pages = NULL;

  nnz = Ap[size];

  Ax = new mumps_scalar[nnz];
  memset(Ax, 0, sizeof(mumps_scalar) * nnz);

  irn = new int[nnz];
  jcn = new int[nnz];
  for (unsigned int i = 0; i < nnz; i++)
  {
    irn[i] = 1;
    jcn[i] = 1;
  }  
}

void MumpsMatrix::free()
{
  _F_
  nnz = 0;
  delete[] Ap; Ap = NULL;
  delete[] Ai; Ai = NULL;
  delete[] Ax; Ax = NULL;
  delete[] irn; irn = NULL;
  delete[] jcn; jcn = NULL;
}

scalar MumpsMatrix::get(unsigned int m, unsigned int n)
{
  _F_
  // Find m-th row in the n-th column.
  int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
  // Return 0 if the entry has not been found.
  if (mid < 0) return 0.0;
  // Otherwise, add offset to the n-th column and return the value.
  if (mid >= 0) mid += Ap[n];
#ifndef HERMES_COMMON_COMPLEX
  return Ax[mid];
#else
  return cplx(Ax[mid].r, Ax[mid].i);
#endif
}

void MumpsMatrix::zero()
{
  _F_
  memset(Ax, 0, sizeof(mumps_scalar) * Ap[size]);
}

void MumpsMatrix::add(unsigned int m, unsigned int n, scalar v)
{
  _F_
  // WARNING: The additional condition v != 0.0 used in (Umfpack)Matrix
  //          produced an error in neutronics-2-group-adapt (although tutorial-07
  //          ran well).
  // Find m-th row in the n-th column.
  int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
  // Make sure we are adding to an existing non-zero entry.
  if (pos < 0) 
    error("Sparse matrix entry not found");
  // Add offset to the n-th column.
  pos += Ap[n];
#ifndef HERMES_COMMON_COMPLEX
  Ax[pos] += v;
#else
  Ax[pos].r += v.real();
  Ax[pos].i += v.imag();
#endif
  irn[pos] = m + 1;  // MUMPS is indexing from 1
  jcn[pos] = n + 1;
}

void MumpsMatrix::add(unsigned int m, unsigned int n, scalar **mat, int *rows, int *cols)
{
  _F_
  for (unsigned int i = 0; i < m; i++)       // rows
    for (unsigned int j = 0; j < n; j++)     // cols
      if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
        add(rows[i], cols[j], mat[i][j]);
}

/// Add a number to each diagonal entry.
void MumpsMatrix::add_to_diagonal(scalar v) 
{
  for (unsigned int i = 0; i < size; i++) {
    add(i, i, v);
  }
};

/// dumping matrix and right-hand side
///
bool MumpsMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  // TODO
  switch (fmt) 
  {
    case DF_NATIVE:
    case DF_PLAIN_ASCII:
      fprintf(file, "%d\n", size);
      fprintf(file, "%d\n", nnz);
      for (unsigned int i = 0; i < nnz; i++)
        fprintf(file, "%d %d " SCALAR_FMT "\n", irn[i], jcn[i], MUMPS_SCALAR(Ax[i]));
      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
      for (unsigned int j = 0; j < size; j++)
        for (unsigned int i = Ap[j]; i < Ap[j + 1]; i++)
#ifndef HERMES_COMMON_COMPLEX          
          fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, MUMPS_SCALAR(Ax[i]));
#else          
          fprintf(file, "%d %d %lf+%lfi\n", Ai[i] + 1, j + 1, MUMPS_SCALAR(Ax[i]));
#endif          
      fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("HERMESX\001", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(&nnz, sizeof(int), 1, file);
      hermes_fwrite(Ap, sizeof(int), size + 1, file);
      hermes_fwrite(Ai, sizeof(int), nnz, file);
      hermes_fwrite(Ax, sizeof(mumps_scalar), nnz, file);
      return true;
    }

    default:
      return false;
  }
}

unsigned int MumpsMatrix::get_matrix_size() const
{
  _F_
  return size;
}

unsigned int MumpsMatrix::get_nnz() const
{
  _F_
  return nnz;
}

/* THIS WAS WRONG
int MumpsMatrix::get_matrix_size() const
{
  _F_
  //           Ax               Ai                 Ap                 
  return (sizeof(scalar) + sizeof(int)) * nnz + sizeof(int)*(size+1)
          + 2 * sizeof(int) * nnz + sizeof(int);
  //          irn, jcn                  nnz                             
}
*/

double MumpsMatrix::get_fill_in() const
{
  _F_
  return Ap[size] / (double) (size * size);
}

void MumpsMatrix::add_matrix(MumpsMatrix* mat){
  _F_
  add_as_block(0,0,mat);
};

void MumpsMatrix::add_to_diagonal_blocks(int num_stages, MumpsMatrix* mat){
  _F_
  int ndof = mat->get_size();
  if (this->get_size() != (unsigned int) num_stages * ndof) 
    error("Incompatible matrix sizes in PetscMatrix::add_to_diagonal_blocks()");

  for (int i = 0; i < num_stages; i++) {
    this->add_as_block(ndof*i, ndof*i, mat);
  }
}

void MumpsMatrix::add_as_block(unsigned int i, unsigned int j, MumpsMatrix* mat){
  _F_
  int idx;
  for (unsigned int col=0;col<mat->get_size();col++){
    for (unsigned int n=mat->Ap[col];n<mat->Ap[col+1];n++){
      idx=find_position(Ai + Ap[col+j], Ap[col + 1 + j] - Ap[col],mat->Ai[n]+i);
      if (idx<0)
        error("Sparse matrix entry not found");
#ifndef HERMES_COMMON_COMPLEX
      Ax[idx]+=mat->Ax[n];
#else      
      Ax[idx].r+=mat->Ax[n].r;
      Ax[idx].i+=mat->Ax[n].i;
#endif
    }
  }
}

  // Applies the matrix to vector_in and saves result to vector_out.
void MumpsMatrix::multiply_with_vector(scalar* vector_in, scalar* vector_out){
  for(unsigned int i=0;i<size;i++){
    vector_out[i]=0;
  }
  scalar a;
  for (unsigned int i=0;i<nnz;i++){
#ifndef HERMES_COMMON_COMPLEX
    a=Ax[i];
#else
    a=cplx(Ax[i].r,Ax[i].i);
#endif
    vector_out[jcn[i]]+=vector_in[irn[i]]*a;
  }
}
  // Multiplies matrix with a scalar.
void MumpsMatrix::multiply_with_scalar(scalar value){
  int n=nnz;
  scalar a;
  for(int i=0;i<n;i++){
#ifndef HERMES_COMMON_COMPLEX
    Ax[i]=Ax[i]*value;
#else
    a=cplx(Ax[i].r,Ax[i].i);
    a=a*value;
    Ax[i].r=a.real();
    Ax[i].i=a.imag();
#endif
  }
}
  // Creates matrix using size, nnz, and the three arrays.
void MumpsMatrix::create(unsigned int size, unsigned int nnz, int* ap, int* ai, scalar* ax){
  this->nnz = nnz;
  this->size = size;
  this->Ap = new unsigned int[size+1]; assert(this->Ap != NULL);
  this->Ai = new int[nnz];    assert(this->Ai != NULL);
  this->Ax = new mumps_scalar[nnz]; assert(this->Ax != NULL);
  irn=new int[nnz];           assert(this->irn !=NULL);     // Row indices.
  jcn=new int[nnz];           assert(this->jcn !=NULL);     // Column indices.

  for (unsigned int i = 0; i < size; i++){
    this->Ap[i] = ap[i];
    for (int j=ap[i];j<ap[i+1];j++) jcn[j]=i;
  }
  this->Ap[size]=ap[size];
  for (unsigned int i = 0; i < nnz; i++) {
#ifndef HERMES_COMMON_COMPLEX
    this->Ax[i] = ax[i]; 
#else
    this->Ax[i].r=ax[i].real();
    this->Ax[i].i=ax[i].imag();
#endif
    this->Ai[i] = ai[i];
    irn[i]=ai[i];
  } 
}
  // Duplicates a matrix (including allocation).
MumpsMatrix* MumpsMatrix::duplicate(){
  MumpsMatrix * nmat=new MumpsMatrix();

  nmat->nnz = nnz;
  nmat->size = size;
  nmat->Ap = new unsigned int[size+1]; assert(nmat->Ap != NULL);
  nmat->Ai = new int[nnz];    assert(nmat->Ai != NULL);
  nmat->Ax = new mumps_scalar[nnz]; assert(nmat->Ax != NULL);
  nmat->irn=new int[nnz];           assert(nmat->irn !=NULL);     // Row indices.
  nmat->jcn=new int[nnz];           assert(nmat->jcn !=NULL);     // Column indices.
  for (unsigned int i = 0;i<nnz;i++){
    nmat->Ai[i]=Ai[i];
    nmat->Ax[i]=Ax[i];
    nmat->irn[i]=irn[i];
    nmat->jcn[i]=jcn[i];
  }
  for (unsigned int i = 0;i<size+1;i++){
    nmat->Ap[i]=Ap[i];
  }
  return nmat;
}

// MumpsVector /////////////////////////////////////////////////////////////////////////////////////

MumpsVector::MumpsVector()
{
  _F_
  v = NULL;
  size = 0;
}

MumpsVector::~MumpsVector()
{
  _F_
  free();
}

void MumpsVector::alloc(unsigned int n)
{
  _F_
  free();
  size = n;
  v = new mumps_scalar[n];
  zero();
}

void MumpsVector::change_sign()
{
  _F_
#ifndef HERMES_COMMON_COMPLEX
  for (unsigned int i = 0; i < size; i++) v[i] *= -1.;
#else
  for (unsigned int i = 0; i < size; i++) {
    v[i].r *= -1.;
    v[i].i *= -1.;
  }
#endif
}

void MumpsVector::zero()
{
  _F_
  memset(v, 0, size * sizeof(mumps_scalar));
}

void MumpsVector::free()
{
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void MumpsVector::set(unsigned int idx, scalar y)
{
  _F_
#ifndef HERMES_COMMON_COMPLEX
  v[idx] = y;
#else
  v[idx].r = y.real();
  v[idx].i = y.imag();
#endif
}

void MumpsVector::add(unsigned int idx, scalar y)
{
  _F_
#ifndef HERMES_COMMON_COMPLEX
  v[idx] += y;
#else
  v[idx].r += y.real();
  v[idx].i += y.imag();
#endif
}

void MumpsVector::add(unsigned int n, unsigned int *idx, scalar *y)
{
  _F_
  for (unsigned int i = 0; i < n; i++) {
#ifndef HERMES_COMMON_COMPLEX
    v[idx[i]] += y[i];
#else
    v[idx[i]].r += y[i].real();
    v[idx[i]].i += y[i].imag();
#endif
  }
}

bool MumpsVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  switch (fmt) 
  {
    case DF_NATIVE:
    case DF_PLAIN_ASCII:
      for (unsigned int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", MUMPS_SCALAR(v[i]));

      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (unsigned int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", MUMPS_SCALAR(v[i]));
      fprintf(file, " ];\n");
      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("HERMESR\001", 1, 8, file);
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

// MUMPS solver ////////////////////////////////////////////////////////////////////////////////////

#ifdef WITH_MUMPS

// Macros allowing to use indices according to the Fortran documentation to index C arrays.
#define ICNTL(I)            icntl[(I)-1]
#define MUMPS_INFO(param,I) (param).infog[(I)-1]
#define INFOG(I)            infog[(I)-1]

// Job definitions according to MUMPS documentation.
#define JOB_INIT                    -1
#define JOB_END                     -2
#define JOB_ANALYZE_FACTORIZE_SOLVE  6
#define JOB_FACTORIZE_SOLVE          5
#define JOB_SOLVE                    3

bool MumpsSolver::check_status()
{
  _F_
  switch (param.INFOG(1)) {
    case 0: return true; // no error
    case -1: warning("Error occured on processor %d", MUMPS_INFO(param, 2)); break;
    // TODO: add the rest according to the MUMPS docs
    default: warning("INFOG(1) = %d", param.INFOG(1)); break;
  }
  return false;
}

bool MumpsSolver::reinit()
{
  _F_
  if (inited)
  {
    // If there is already an instance of MUMPS running, 
    // terminate it.
    param.job = JOB_END;
    MUMPS(&param);
  }
  
  param.job = JOB_INIT;
  param.par = 1; // host also performs calculations
  param.sym = 0; // 0 = unsymmetric
  param.comm_fortran=USE_COMM_WORLD;
  
  MUMPS(&param);
  inited = check_status();
  
  if (inited)
  {
    // No printings.
    param.ICNTL(1) = -1;
    param.ICNTL(2) = -1;
    param.ICNTL(3) = -1;
    param.ICNTL(4) = 0;
    
    param.ICNTL(20) = 0; // centralized dense RHS
    param.ICNTL(21) = 0; // centralized dense solution
    
    // Specify the matrix.
    param.n = m->size;
    param.nz = m->nnz;
    param.irn = m->irn;
    param.jcn = m->jcn;
    param.a = m->Ax;
  }
  
  return inited;
}

#endif

MumpsSolver::MumpsSolver(MumpsMatrix *m, MumpsVector *rhs) :
  LinearSolver(), m(m), rhs(rhs)
{
  _F_
#ifdef WITH_MUMPS
  inited = false;
  
  // Initial values for some fields of the MUMPS_STRUC structure that may be accessed
  // before MUMPS has been initialized.
  param.rhs = NULL;
  param.INFOG(33) = -999; // see the case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING 
                          // in setup_factorization()
#else
  error(MUMPS_NOT_COMPILED);
#endif
}

MumpsSolver::~MumpsSolver()
{
  _F_
#ifdef WITH_MUMPS
  // Terminate the current instance of MUMPS.
  if (inited)
  {
    param.job = JOB_END;
    MUMPS(&param);
  }
  
  if (param.rhs != NULL) delete [] param.rhs;
#endif
}

bool MumpsSolver::solve()
{
  _F_
#ifdef WITH_MUMPS
  bool ret = false;
  assert(m != NULL);
  assert(rhs != NULL);

  TimePeriod tmr;

  // Prepare the MUMPS data structure with input for the solver driver 
  // (according to the chosen factorization reuse strategy), as well as
  // the system matrix.
  if ( !setup_factorization() )
  {
    warning("LU factorization could not be completed.");
    return false;
  }
  
  // Specify the right-hand side (will be replaced by the solution).
  param.rhs = new mumps_scalar[m->size];
  memcpy(param.rhs, rhs->v, m->size * sizeof(mumps_scalar));
  
  // Do the jobs specified in setup_factorization().
  MUMPS(&param);
  
  ret = check_status();

  if (ret) 
  {
    delete [] sln;
    sln = new scalar[m->size];
#ifndef HERMES_COMMON_COMPLEX
    for (unsigned int i = 0; i < rhs->size; i++)
      sln[i] = param.rhs[i];
#else
    for (unsigned int i = 0; i < rhs->size; i++)
      sln[i] = cplx(param.rhs[i].r, param.rhs[i].i);
#endif
  }

  tmr.tick();
  time = tmr.accumulated();

  delete [] param.rhs;
  param.rhs = NULL;

  return ret;
#else
  return false;
#endif
}

bool MumpsSolver::setup_factorization()
{
  _F_
#ifdef WITH_MUMPS
  // When called for the first time, all three phases (analysis, factorization,
  // solution) must be performed. 
  int eff_fact_scheme = factorization_scheme;
  if (!inited)
    if( factorization_scheme == HERMES_REUSE_MATRIX_REORDERING || 
        factorization_scheme == HERMES_REUSE_FACTORIZATION_COMPLETELY )
      eff_fact_scheme = HERMES_FACTORIZE_FROM_SCRATCH;
  
  switch (eff_fact_scheme)
  {
    case HERMES_FACTORIZE_FROM_SCRATCH: 
      // (Re)initialize new instance.
      reinit();
      
      // Let MUMPS decide when and how to compute matrix reordering and scaling.
      param.ICNTL(6) = 7;
      param.ICNTL(8) = 77;
      param.job = JOB_ANALYZE_FACTORIZE_SOLVE;
      
      break;
    case HERMES_REUSE_MATRIX_REORDERING:
      // Let MUMPS reuse results of the symbolic analysis and perform 
      // scaling during each factorization (values 1-8 may be set here, 
      // corresponding to different scaling algorithms during factorization; 
      // see the MUMPS documentation for details).
      param.ICNTL(8) = 7; 
      param.job = JOB_FACTORIZE_SOLVE;
      
      break;
    case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
      // Perform scaling along with reordering during the symbolic analysis phase
      // and then reuse it during subsequent factorizations. New instance of MUMPS
      // has to be created before the analysis phase. 
      if (param.INFOG(33) != -2)
      {
        reinit();
        param.ICNTL(6) = 5;
        param.job = JOB_ANALYZE_FACTORIZE_SOLVE;
        // After analysis is done, INFOG(33) will be set to -2 by MUMPS.
      }
      else
      {
        param.job = JOB_FACTORIZE_SOLVE;
      }
      break;
    case HERMES_REUSE_FACTORIZATION_COMPLETELY:
      param.job = JOB_SOLVE;
      break;
  }
  
  return true;
#else
  return false;
#endif
}
