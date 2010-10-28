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

#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  #define MUMPS			dmumps_c
  #define MUMPS_STRUCT	DMUMPS_STRUC_C
#else
  #define MUMPS			zmumps_c
  #define MUMPS_STRUCT	ZMUMPS_STRUC_C
#endif

#define USE_COMM_WORLD -987654

#ifdef WITH_MUMPS

extern "C" {
  extern void MUMPS(MUMPS_STRUCT *idptr);
}

#else
#endif

static
int find_position(int *Ai, int Alen, int idx) {
  _F_
  if (idx >= 0) {
    register int lo = 0, hi = Alen - 1, mid;

    while (1) {
      mid = (lo + hi) >> 1;

      if (idx < Ai[mid]) hi = mid - 1;
      else if (idx > Ai[mid]) lo = mid + 1;
      else break;

      if (lo > hi) error("Sparse matrix entry not found.");
    }

    return mid;
  }
  return -1;
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
  assert(size != 0);

  // initialize the arrays Ap and Ai
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
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  Ax = new scalar[nnz];
  memset(Ax, 0, sizeof(scalar) * nnz);
#else
  Ax = new ZMUMPS_COMPLEX[nnz];
  memset(Ax, 0, sizeof(ZMUMPS_COMPLEX) * nnz);
#endif

  irn = new int[nnz];
  memset(irn, 0, sizeof(int) * nnz);
  jcn = new int[nnz];
  memset(jcn, 0, sizeof(int) * nnz);
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

scalar MumpsMatrix::get(int m, int n)
{
  _F_
  int mid = Ap[n] + find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  return Ax[mid];
#else
  return cplx(Ax[mid].r, Ax[mid].i);
#endif
}

void MumpsMatrix::zero()
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  memset(Ax, 0, sizeof(scalar) * Ap[size]);
#else
  memset(Ax, 0, sizeof(ZMUMPS_COMPLEX) * Ap[size]);
#endif
}

void MumpsMatrix::add(int m, int n, scalar v)
{
  _F_
  if (m >= 0 && n >= 0) {		// ignore dirichlet DOFs
    int pos = Ap[n] + find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    Ax[pos] += v;
#else
    Ax[pos].r += v.real();
    Ax[pos].i += v.imag();
#endif
    irn[pos] = m + 1;			// MUMPS is indexing from 1
    jcn[pos] = n + 1;
  }
}

void MumpsMatrix::add(int m, int n, scalar **mat, int *rows, int *cols)
{
  _F_
  for (int i = 0; i < m; i++)				// rows
    for (int j = 0; j < n; j++)			// cols
      add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool MumpsMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  // TODO
  switch (fmt) 
  {
    case DF_NATIVE:
      fprintf(file, "%d\n", size);
      fprintf(file, "%d\n", nnz);
      for (int i = 0; i < nnz; i++)
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
        fprintf(file, "%d %d %lf\n", irn[i], jcn[i], Ax[i]);
#else
        fprintf(file, "%d %d (%lf,%lf)\n", irn[i], jcn[i], Ax[i].r, Ax[i].i);
#endif
      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
      for (int j = 0; j < size; j++)
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
          fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, SCALAR(Ax[i]));
      fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DX\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      int nnz = Ap[size];
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

int MumpsMatrix::get_matrix_size() const
{
  _F_
  return (sizeof(scalar) + 2 * sizeof(int)) * nnz;
}

double MumpsMatrix::get_fill_in() const
{
  _F_
  return Ap[size] / (double) (size * size);
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

void MumpsVector::alloc(int n)
{
  _F_
  free();
  size = n;
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  v = new scalar[n];
#else
  v = new ZMUMPS_COMPLEX[n];
#endif
  zero();
}

void MumpsVector::zero()
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  memset(v, 0, size * sizeof(scalar));
#else
  memset(v, 0, size * sizeof(ZMUMPS_COMPLEX));
#endif
}

void MumpsVector::free()
{
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void MumpsVector::set(int idx, scalar y)
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] = y;
#else
  if (idx >= 0) {
    v[idx].r = y.real();
    v[idx].i = y.imag();
  }
#endif
}

void MumpsVector::add(int idx, scalar y)
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] += y;
#else
  if (idx >= 0) {
    v[idx].r += y.real();
    v[idx].i += y.imag();
  }
#endif
}

void MumpsVector::add(int n, int *idx, scalar *y)
{
  _F_
  for (int i = 0; i < n; i++)
    if (idx[i] >= 0) {
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
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
  switch (fmt) {
    case DF_NATIVE:
      for (int i = 0; i < size; i++)
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
        fprintf(file, "%lf\n", v[i]);
#else
        fprintf(file, "(%lf,%lf)\n", v[i].r, v[i].i);
#endif
      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (int i = 0; i < size; i++)
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
        fprintf(file, SCALAR_FMT "\n", SCALAR(v[i]));
#else
      fprintf(file, "(%lf, %lf)\n", v[i].r, v[i].i);
#endif
      fprintf(file, " ];\n");
      return true;

    case DF_HERMES_BIN: {
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

// MUMPS solver ////////////////////////////////////////////////////////////////////////////////////

MumpsSolver::MumpsSolver(MumpsMatrix *m, MumpsVector *rhs) :
  LinearSolver(), m(m), rhs(rhs)
{
  _F_
#ifdef WITH_MUMPS
#else
  error(MUMPS_NOT_COMPILED);
#endif
}

MumpsSolver::~MumpsSolver()
{
  _F_
#ifdef WITH_MUMPS
  //if (m != NULL) delete m;
  //if (rhs != NULL) delete rhs;
#endif
}

#ifdef WITH_MUMPS

// macro s.t. indices match Fortran documentation
#define ICNTL(I)						icntl[(I)-1]
#define MUMPS_INFO(id, I)				id->infog[(I)-1]
#define INFOG(I)						infog[(I)-1]

#define JOB_INIT						-1
#define JOB_END							-2

static bool check_status(MUMPS_STRUCT *id)
{
  _F_
  switch (id->INFOG(1)) {
    case 0: return true; // no error
    case -1: warning("Error occured on processor %d", MUMPS_INFO(id, 2)); break;
    // TODO: add the rest according to the MUMPS docs
    default: warning("INFOG(1) = %d", id->INFOG(1)); break;
  }
  return false;
}

#endif

bool MumpsSolver::solve()
{
  _F_
#ifdef WITH_MUMPS
  bool ret = false;
  assert(m != NULL);
  assert(rhs != NULL);

  TimePeriod tmr;

  MUMPS_STRUCT id;

  // Initialize a MUMPS instance
  id.job = JOB_INIT;
  id.par = 1; // host also performs calculations
  id.sym = 0; // 0 = unsymmetric
  id.comm_fortran=USE_COMM_WORLD;
  
  MUMPS(&id);
  check_status(&id);

  // matrix
  id.n = m->size;
  id.nz = m->nnz;
  id.irn = m->irn;
  id.jcn = m->jcn;
  id.a = m->Ax;

  // right-hand side
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  id.rhs = new double[m->size];
  memcpy(id.rhs, rhs->v, m->size * sizeof(double));
#else
  id.rhs = new ZMUMPS_COMPLEX[m->size];
  memcpy(id.rhs, rhs->v, m->size * sizeof(ZMUMPS_COMPLEX));
#endif

  // No printings
  id.ICNTL(1) = -1;
  id.ICNTL(2) = -1;
  id.ICNTL(3) = -1;
  id.ICNTL(4) = 0;

  id.ICNTL(20) = 0; // centralized dense RHS
  id.ICNTL(21) = 0; // centralized dense solution
  
  id.job = 6; // 6 means Analysis + factorization + solve (FIXME: remove magic constant - see MUMPS docs)
  MUMPS(&id);
  
  ret = check_status(&id);

  if (ret) {
    delete [] sln;
    sln = new scalar[m->size];
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    for (int i = 0; i < rhs->size; i++)
      sln[i] = id.rhs[i];
#else
    for (int i = 0; i < rhs->size; i++)
      sln[i] = cplx(id.rhs[i].r, id.rhs[i].i);
#endif
  }

  // Terminate/free current instance
  id.job = JOB_END;
  MUMPS(&id);

  tmr.tick();
  time = tmr.accumulated();

  delete [] id.rhs;

  return ret;
#else
  return false;
#endif
}
