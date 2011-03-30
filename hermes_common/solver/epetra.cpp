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

// EpetraMatrix<Scalar> ////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_EPETRA
  // A communicator for Epetra objects (serial version)
  static Epetra_SerialComm seq_comm;
#endif

template<typename Scalar>
EpetraMatrix<Scalar>::EpetraMatrix()
{
  _F_
#ifdef HAVE_EPETRA
  this->mat = NULL;
  this->mat_im = NULL;
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
template<typename Scalar>
EpetraMatrix<Scalar>::EpetraMatrix(Epetra_RowMatrix &op)
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

template<typename Scalar>
EpetraMatrix<Scalar>::~EpetraMatrix()
{
  _F_
#ifdef HAVE_EPETRA
  free();
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::prealloc(unsigned int n)
{
  _F_
#ifdef HAVE_EPETRA
  this->size = n;
  // alloc trilinos structs
  std_map = new Epetra_Map(n, 0, seq_comm); MEM_CHECK(std_map);
  grph = new Epetra_CrsGraph(Copy, *std_map, 0); MEM_CHECK(grph);
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
{
  _F_
#ifdef HAVE_EPETRA
  int col_to_pass = col;
  grph->InsertGlobalIndices(row, 1, &col_to_pass);
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::finish()
{
  _F_
#ifdef HAVE_EPETRA
  mat->FillComplete();
#ifdef HERMES_COMMON_COMPLEX
  mat_im->FillComplete();
#endif
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::alloc()
{
  _F_
#ifdef HAVE_EPETRA
  grph->FillComplete();
  // create the matrix
  mat = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat);
#ifdef HERMES_COMMON_COMPLEX
  mat_im = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat_im);
#endif
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::free()
{
  _F_
#ifdef HAVE_EPETRA
  if (owner) {
    delete mat; mat = NULL;
#ifdef HERMES_COMMON_COMPLEX
    delete mat_im; mat_im = NULL;
#endif
    delete grph; grph = NULL;
    delete std_map; std_map = NULL;
  }
#endif
}

template<typename Scalar>
Scalar EpetraMatrix<Scalar>::get(unsigned int m, unsigned int n)
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

template<typename Scalar>
int EpetraMatrix<Scalar>::get_num_row_entries(unsigned int row)
{
  _F_
#ifdef HAVE_EPETRA
    return mat->NumGlobalEntries(row);
#else
  return 0;
#endif
}

template<typename Scalar>
void EpetraMatrix<Scalar>::extract_row_copy(unsigned int row, unsigned int len, unsigned int &n_entries, double *vals, unsigned int *idxs)
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

template<typename Scalar>
void EpetraMatrix<Scalar>::zero()
{
  _F_
#ifdef HAVE_EPETRA
  mat->PutScalar(0.0);
#ifdef HERMES_COMMON_COMPLEX
  mat_im->PutScalar(0.0);
#endif
#endif
}

template<>
void EpetraMatrix<double>::add(unsigned int m, unsigned int n, double v)
{
  _F_
#ifdef HAVE_EPETRA
  if (v != 0.0) {		// ignore zero values
    int n_to_pass = n;
    int ierr = mat->SumIntoGlobalValues(m, 1, &v, &n_to_pass);
    if (ierr != 0) error("Failed to insert into Epetra matrix");
  }
#endif
}

template<>
void EpetraMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
{
  _F_
#ifdef HAVE_EPETRA
  if (v != 0.0) {		// ignore zero values
    double v_r = std::real<double>(v);
    int n_to_pass = n;
    int ierr = mat->SumIntoGlobalValues(m, 1, &v_r, &n_to_pass);
    assert(ierr == 0);
    double v_i = std::imag<double>(v);
    ierr = mat_im->SumIntoGlobalValues(m, 1, &v_i, &n_to_pass);
    assert(ierr == 0);
  }
#endif
}

/// Add a number to each diagonal entry.
template<typename Scalar>
void EpetraMatrix<Scalar>::add_to_diagonal(Scalar v) 
{
  for (unsigned int i=0; i < this->size; i++) {
    add(i, i, v);
  }
};

template<typename Scalar>
void EpetraMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
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
template<typename Scalar>
bool EpetraMatrix<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  return false;
}

template<typename Scalar>
unsigned int EpetraMatrix<Scalar>::get_matrix_size() const
{
  _F_
#ifdef HAVE_EPETRA
  return this->size;
#else
  return -1;
#endif
}

template<typename Scalar>
double EpetraMatrix<Scalar>::get_fill_in() const
{
  _F_
#ifdef HAVE_EPETRA
  return mat->NumGlobalNonzeros() / ((double)size*size);
#else
  return -1;
#endif
}

template<typename Scalar>
unsigned int EpetraMatrix<Scalar>::get_nnz() const
{
  _F_
#ifdef HAVE_EPETRA
  return mat->NumGlobalNonzeros();
#else
  return -1;
#endif
}

// EpetraVector<Scalar> ////////////////////////////////////////////////////////////////////////////////////

template<typename Scalar>
EpetraVector<Scalar>::EpetraVector()
{
  _F_
#ifdef HAVE_EPETRA
  this->std_map = NULL;
  this->vec = NULL;
  this->vec_im = NULL;
  this->size = 0;
  this->owner = true;
#endif
}

#ifdef HAVE_EPETRA
template<typename Scalar>
EpetraVector<Scalar>::EpetraVector(const Epetra_Vector &v)
{
  _F_
  this->vec = (Epetra_Vector *) &v;
  this->std_map = (Epetra_BlockMap *) &v.Map();
  this->size = v.MyLength();
  this->owner = false;
}
#endif

template<typename Scalar>
EpetraVector<Scalar>::~EpetraVector()
{
  _F_
#ifdef HAVE_EPETRA
  if (owner) free();
#endif
}

template<typename Scalar>
void EpetraVector<Scalar>::alloc(unsigned int n)
{
  _F_
#ifdef HAVE_EPETRA
  free();
  this->size = n;
  std_map = new Epetra_Map(this->size, 0, seq_comm); MEM_CHECK(std_map);
  vec = new Epetra_Vector(*std_map); MEM_CHECK(vec);
  vec_im = new Epetra_Vector(*std_map); MEM_CHECK(vec_im);
  zero();
#endif
}

template<typename Scalar>
void EpetraVector<Scalar>::zero()
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < this->size; i++) (*vec)[i] = 0.0;
  for (unsigned int i = 0; i < this->size; i++) (*vec_im)[i] = 0.0;
#endif
}

template<typename Scalar>
void EpetraVector<Scalar>::change_sign()
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < this->size; i++) (*vec)[i] *= -1.;
  for (unsigned int i = 0; i < this->size; i++) (*vec_im)[i] *= -1.;
#endif
}

template<typename Scalar>
void EpetraVector<Scalar>::free()
{
  _F_
#ifdef HAVE_EPETRA
  delete std_map; std_map = NULL;
  delete vec; vec = NULL;
  delete vec_im; vec_im = NULL;
  this->size = 0;
#endif
}

template<>
void EpetraVector<double>::set(unsigned int idx, double y)
{
  _F_
#ifdef HAVE_EPETRA
  (*vec)[idx] = y;
#endif
}

template<>
void EpetraVector<std::complex<double> >::set(unsigned int idx, std::complex<double> y)
{
  _F_
#ifdef HAVE_EPETRA
  (*vec)[idx] = std::real(y);
  (*vec_im)[idx] = std::imag(y);
#endif
}

template<>
void EpetraVector<double>::add(unsigned int idx, double y)
{
  _F_
#ifdef HAVE_EPETRA
  (*vec)[idx] += y;
#endif
}

template<>
void EpetraVector<std::complex<double> >::add(unsigned int idx, std::complex<double> y)
{
  _F_
#ifdef HAVE_EPETRA
  (*vec)[idx] += std::real(y);
  (*vec_im)[idx] += std::imag(y);
#endif
}

template<typename Scalar>
void EpetraVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
{
  _F_
#ifdef HAVE_EPETRA
  for (unsigned int i = 0; i < n; i++)
    add(idx[i], y[i]);
#endif
}

template<typename Scalar>
bool EpetraVector<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  return false;
}
template class HERMES_API EpetraMatrix<double>;
template class HERMES_API EpetraMatrix<std::complex<double> >;
template class HERMES_API EpetraVector<double>;
template class HERMES_API EpetraVector<std::complex<double> >;