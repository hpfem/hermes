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

#include "epetra.h"
#include "../error.h"
#include "../callstack.h"

// EpetraMatrix ////////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRA
  // A communicator for Epetra objects (serial version)
  static Epetra_SerialComm seq_comm;
#endif

EpetraMatrix::EpetraMatrix()
{
  _F_
#ifdef HAVE_EPETRA
  this->mat = NULL;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  this->mat_im = NULL;
#endif
  this->grph = NULL;
  this->std_map = NULL;
  this->owner = true;

  this->row_storage = true;
  this->col_storage = false;
#else
  error(EPETRA_NOT_COMPILED);
#endif
}

#ifdef HAVE_EPETRA
EpetraMatrix::EpetraMatrix(Epetra_RowMatrix &op)
{
  _F_
  this->mat = dynamic_cast<Epetra_CrsMatrix *>(&op);
  assert(mat != NULL);
  this->grph = (Epetra_CrsGraph *) &this->mat->Graph();
  this->std_map = (Epetra_BlockMap *) &this->grph->Map();
  this->owner = false;

  this->row_storage = true;
  this->col_storage = false;
}
#endif

EpetraMatrix::~EpetraMatrix()
{
  _F_
#ifdef HAVE_EPETRA
  free();
#endif
}

void EpetraMatrix::prealloc(unsigned int n)
{
  _F_
#ifdef HAVE_EPETRA
  this->size = n;
  // alloc trilinos structs
  std_map = new Epetra_Map(n, 0, seq_comm); MEM_CHECK(std_map);
  grph = new Epetra_CrsGraph(Copy, *std_map, 0); MEM_CHECK(grph);
#endif
}

void EpetraMatrix::pre_add_ij(unsigned int row, unsigned int col)
{
  _F_
#ifdef HAVE_EPETRA
  int col_to_pass = col;
  grph->InsertGlobalIndices(row, 1, &col_to_pass);
#endif
}

void EpetraMatrix::finish()
{
  _F_
#ifdef HAVE_EPETRA
  mat->FillComplete();
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  mat_im->FillComplete();
#endif
#endif
}

void EpetraMatrix::alloc()
{
  _F_
#ifdef HAVE_EPETRA
  grph->FillComplete();
  // create the matrix
  mat = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat);
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  mat_im = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat_im);
#endif
#endif
}

void EpetraMatrix::free()
{
  _F_
#ifdef HAVE_EPETRA
  if (owner) {
    delete mat; mat = NULL;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
    delete mat_im; mat_im = NULL;
#endif
    delete grph; grph = NULL;
    delete std_map; std_map = NULL;
  }
#endif
}

scalar EpetraMatrix::get(unsigned int m, unsigned int n)
{
  _F_
#ifdef HAVE_EPETRA
    int n_entries = mat->NumGlobalEntries(m);
    std::vector<double> vals(n_entries);
    std::vector<int> idxs(n_entries);
    mat->ExtractGlobalRowCopy(m, n_entries, n_entries, &vals[0], &idxs[0]);
    for (int i = 0; i < n_entries; i++) if (idxs[i] == (int)n) return vals[i];
#endif
    return 0.0;
}

int EpetraMatrix::get_num_row_entries(unsigned int row)
{
  _F_
#ifdef HAVE_EPETRA
    return mat->NumGlobalEntries(row);
#else
  return 0;
#endif
}

void EpetraMatrix::extract_row_copy(unsigned int row, unsigned int len, unsigned int &n_entries, double *vals, unsigned int *idxs)
{
  _F_
#ifdef HAVE_EPETRA
  int* idxs_to_pass = new int[len];
  for(unsigned int i = 0; i < len; i++)
    idxs_to_pass[i] = idxs[i];
  int n_entries_to_pass = n_entries;
  mat->ExtractGlobalRowCopy(row, len, n_entries_to_pass, vals, idxs_to_pass);
  delete [] idxs_to_pass;
#endif
}

void EpetraMatrix::zero()
{
  _F_
#ifdef HAVE_EPETRA
  mat->PutScalar(0.0);
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  mat_im->PutScalar(0.0);
#endif
#endif
}

void EpetraMatrix::add(unsigned int m, unsigned int n, scalar v)
{
  _F_
#ifdef HAVE_EPETRA
  if (v != 0.0) {		// ignore zero values
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    int n_to_pass = n;
    int ierr = mat->SumIntoGlobalValues(m, 1, &v, &n_to_pass);
    if (ierr != 0) error("Failed to insert into Epetra matrix");
#else
    double v_r = std::real<double>(v);
    int n_to_pass = n;
    int ierr = mat->SumIntoGlobalValues(m, 1, &v_r, &n_to_pass);
    assert(ierr == 0);
    double v_i = std::imag<double>(v);
    ierr = mat_im->SumIntoGlobalValues(m, 1, &v_i, &n_to_pass);
    assert(ierr == 0);
#endif
  }
#endif
}

/// Add a number to each diagonal entry.
void EpetraMatrix::add_to_diagonal(scalar v) 
{
  for (unsigned int i=0; i<size; i++) {
    add(i, i, v);
  }
};

void EpetraMatrix::add(unsigned int m, unsigned int n, scalar **mat, int *rows, int *cols)
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < m; i++)				// rows
    for (unsigned int j = 0; j < n; j++)			// cols
      if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
        add(rows[i], cols[j], mat[i][j]);
#endif
}

/// dumping matrix and right-hand side
///
bool EpetraMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  return false;
}

unsigned int EpetraMatrix::get_matrix_size() const
{
  _F_
#ifdef HAVE_EPETRA
  return size;
#else
  return -1;
#endif
}

double EpetraMatrix::get_fill_in() const
{
  _F_
#ifdef HAVE_EPETRA
  return mat->NumGlobalNonzeros() / ((double)size*size);
#else
  return -1;
#endif
}

unsigned int EpetraMatrix::get_nnz() const
{
  _F_
#ifdef HAVE_EPETRA
  return mat->NumGlobalNonzeros();
#else
  return -1;
#endif
}

// EpetraVector ////////////////////////////////////////////////////////////////////////////////////

EpetraVector::EpetraVector()
{
  _F_
#ifdef HAVE_EPETRA
  this->std_map = NULL;
  this->vec = NULL;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  this->vec_im = NULL;
#endif
  this->size = 0;
  this->owner = true;
#endif
}

#ifdef HAVE_EPETRA
EpetraVector::EpetraVector(const Epetra_Vector &v)
{
  _F_
  this->vec = (Epetra_Vector *) &v;
  this->std_map = (Epetra_BlockMap *) &v.Map();
  this->size = v.MyLength();
  this->owner = false;
}
#endif

EpetraVector::~EpetraVector()
{
  _F_
#ifdef HAVE_EPETRA
  if (owner) free();
#endif
}

void EpetraVector::alloc(unsigned int n)
{
  _F_
#ifdef HAVE_EPETRA
  free();
  size = n;
  std_map = new Epetra_Map(size, 0, seq_comm); MEM_CHECK(std_map);
  vec = new Epetra_Vector(*std_map); MEM_CHECK(vec);
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  vec_im = new Epetra_Vector(*std_map); MEM_CHECK(vec_im);
#endif
  zero();
#endif
}

void EpetraVector::zero()
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < size; i++) (*vec)[i] = 0.0;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  for (unsigned int i = 0; i < size; i++) (*vec_im)[i] = 0.0;
#endif
#endif
}

void EpetraVector::change_sign()
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < size; i++) (*vec)[i] *= -1.;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  for (unsigned int i = 0; i < size; i++) (*vec_im)[i] *= -1.;
#endif
#endif
}

void EpetraVector::free()
{
  _F_
#ifdef HAVE_EPETRA
  delete std_map; std_map = NULL;
  delete vec; vec = NULL;
#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)
  delete vec_im; vec_im = NULL;
#endif
  size = 0;
#endif
}

void EpetraVector::set(unsigned int idx, scalar y)
{
  _F_
#ifdef HAVE_EPETRA
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  (*vec)[idx] = y;
#else
  (*vec)[idx] = std::real(y);
  (*vec_im)[idx] = std::imag(y);
#endif
#endif
}

void EpetraVector::add(unsigned int idx, scalar y)
{
  _F_
#ifdef HAVE_EPETRA
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  (*vec)[idx] += y;
#else
  (*vec)[idx] += std::real(y);
  (*vec_im)[idx] += std::imag(y);
#endif
#endif
}

void EpetraVector::add(unsigned int n, unsigned int *idx, scalar *y)
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < n; i++)
    add(idx[i], y[i]);
#endif
}

bool EpetraVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  return false;
}
