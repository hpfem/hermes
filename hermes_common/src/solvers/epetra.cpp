// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file epetra.cpp
\brief EpetraMatrix and EpetraVector storage classes for Amesos, AztecOO, ... .
*/
#include "config.h"
#ifdef HAVE_EPETRA
#include "epetra.h"
#include "error.h"
#include "callstack.h"

using namespace Hermes::Error;

namespace Hermes 
{
  namespace Algebra 
  {
    /// \brief A communicator for Epetra objects (serial version).
    static Epetra_SerialComm seq_comm;

    template<typename Scalar>
    EpetraMatrix<Scalar>::EpetraMatrix()
    {
      _F_;
      this->mat = NULL;
      this->mat_im = NULL;
      this->grph = NULL;
      this->std_map = NULL;
      this->owner = true;

      this->row_storage = true;
      this->col_storage = false;
    }

    template<typename Scalar>
    EpetraMatrix<Scalar>::EpetraMatrix(Epetra_RowMatrix &op)
    {
      _F_;
      this->mat = dynamic_cast<Epetra_CrsMatrix *>(&op);
      assert(mat != NULL);
      this->grph = (Epetra_CrsGraph *) &this->mat->Graph();
      this->std_map = (Epetra_BlockMap *) &this->grph->Map();
      this->owner = false;

      this->row_storage = true;
      this->col_storage = false;
    }

    template<typename Scalar>
    EpetraMatrix<Scalar>::~EpetraMatrix()
    {
      _F_;
      free();
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::prealloc(unsigned int n)
    {
      _F_;
      this->size = n;
      // alloc trilinos structs
      std_map = new Epetra_Map(n, 0, seq_comm); MEM_CHECK(std_map);
      grph = new Epetra_CrsGraph(Copy, *std_map, 0); MEM_CHECK(grph);
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      _F_;
      int col_to_pass = col;
      grph->InsertGlobalIndices(row, 1, &col_to_pass);
    }

    template<>
    void EpetraMatrix<double>::finish()
    {
      _F_;
      mat->FillComplete();
    }

    template<>
    void EpetraMatrix<std::complex<double> >::finish()
    {
      _F_;
      mat->FillComplete();
      mat_im->FillComplete();
    }

    template<>
    void EpetraMatrix<double>::alloc()
    {
      _F_;
      grph->FillComplete();
      // create the matrix
      mat = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat);
    }

    template<>
    void EpetraMatrix<std::complex<double> >::alloc()
    {
      _F_;
      grph->FillComplete();
      // create the matrix
      mat = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat);
      mat_im = new Epetra_CrsMatrix(Copy, *grph); MEM_CHECK(mat_im);
    }

    template<>
    void EpetraMatrix<double>::free()
    {
      _F_;
      if (owner) 
      {
        delete mat; mat = NULL;
        delete grph; grph = NULL;
        delete std_map; std_map = NULL;
      }
    }

    template<>
    void EpetraMatrix<std::complex<double> >::free()
    {
      _F_;
      if (owner) 
      {
        delete mat; mat = NULL;
        delete mat_im; mat_im = NULL;
        delete grph; grph = NULL;
        delete std_map; std_map = NULL;
      }
    }

    template<typename Scalar>
    Scalar EpetraMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      _F_;
      int n_entries = mat->NumGlobalEntries(m);
      Hermes::vector<double> vals(n_entries);
      Hermes::vector<int> idxs(n_entries);
      mat->ExtractGlobalRowCopy(m, n_entries, n_entries, &vals[0], &idxs[0]);
      for (int i = 0; i < n_entries; i++) 
        if (idxs[i] == (int)n) 
          return vals[i];
      return 0.0;
    }

    template<typename Scalar>
    int EpetraMatrix<Scalar>::get_num_row_entries(unsigned int row)
    {
      _F_;
      return mat->NumGlobalEntries(row);
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::extract_row_copy(unsigned int row, unsigned int len, unsigned int &n_entries, double *vals, unsigned int *idxs)
    {
      _F_;
      int* idxs_to_pass = new int[len];
      for(unsigned int i = 0; i < len; i++)
        idxs_to_pass[i] = idxs[i];
      int n_entries_to_pass = n_entries;
      mat->ExtractGlobalRowCopy(row, len, n_entries_to_pass, vals, idxs_to_pass);
      delete [] idxs_to_pass;
    }

    template<>
    void EpetraMatrix<double>::zero()
    {
      _F_;
      mat->PutScalar(0.0);
    }

    template<>
    void EpetraMatrix<std::complex<double> >::zero()
    {
      _F_;
      mat->PutScalar(0.0);
      mat_im->PutScalar(0.0);
    }

    template<>
    void EpetraMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      _F_;
      if (v != 0.0) {		// ignore zero values
        int n_to_pass = n;
        int ierr = mat->SumIntoGlobalValues(m, 1, &v, &n_to_pass);
        if (ierr != 0) error("Failed to insert into Epetra matrix");
      }
    }

    template<>
    void EpetraMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      _F_;
      if (v != 0.0) {		// ignore zero values
        double v_r = std::real<double>(v);
        int n_to_pass = n;
        int ierr = mat->SumIntoGlobalValues(m, 1, &v_r, &n_to_pass);
        if (ierr != 0) error("Failed to insert into Epetra matrix");
        assert(ierr == 0);
        double v_i = std::imag<double>(v);
        ierr = mat_im->SumIntoGlobalValues(m, 1, &v_i, &n_to_pass);
        assert(ierr == 0);
      }
    }

    /// Add a number to each diagonal entry.
    template<typename Scalar>
    void EpetraMatrix<Scalar>::add_to_diagonal(Scalar v) 
    {
      for (unsigned int i=0; i < this->size; i++) 
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void EpetraMatrix<Scalar>::add_to_diagonal_blocks(int num_stages, EpetraMatrix<Scalar>* mat_block)
    {
      _F_;
      int ndof = mat_block->get_size();
      if (this->get_size() != (unsigned int) num_stages * ndof) 
        error("Incompatible matrix sizes in CSCMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++) 
      {
        this->add_as_block(ndof*i, ndof*i, mat_block);
      }
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::add_as_block(unsigned int i, unsigned int j, EpetraMatrix<Scalar>* mat)
    {
      if ((this->get_size() < i+mat->get_size() )||(this->get_size() < j+mat->get_size() )) 
        error("Incompatible matrix sizes in Epetra<Scalar>::add_as_block()");
      unsigned int block_size=mat->get_size();
      for (unsigned int r=0;r<block_size;r++)
      {
        for (unsigned int c=0;c<block_size;c++)
        {
          this->add(i+r,j+c,mat->get(r,c));
        }
      }
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out){
      for (unsigned int i=0;i<this->size;i++) //probably can be optimized by use native vectors
      {
        vector_out[i]=0;
        for (unsigned int j=0;j<this->size;j++)
        {
          vector_out[i]+=vector_in[j]*get(i,j);
        }
      }
    }

    template<typename Scalar>
    void EpetraMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      _F_;
      for (unsigned int i = 0; i < m; i++)				// rows
        for (unsigned int j = 0; j < n; j++)			// cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    template<typename Scalar>
    bool EpetraMatrix<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      return false;
    }

    template<typename Scalar>
    unsigned int EpetraMatrix<Scalar>::get_matrix_size() const
    {
      _F_;
      return this->size;
    }

    template<typename Scalar>
    double EpetraMatrix<Scalar>::get_fill_in() const
    {
      _F_;
      return mat->NumGlobalNonzeros() / ((double)this->size*this->size);
    }

    template<typename Scalar>
    unsigned int EpetraMatrix<Scalar>::get_nnz() const
    {
      _F_;
      return mat->NumGlobalNonzeros();
    }

    template<typename Scalar>
    EpetraVector<Scalar>::EpetraVector()
    {
      _F_;
      this->std_map = NULL;
      this->vec = NULL;
      this->vec_im = NULL;
      this->size = 0;
      this->owner = true;
    }

    template<typename Scalar>
    EpetraVector<Scalar>::EpetraVector(const Epetra_Vector &v)
    {
      _F_;
      this->vec = (Epetra_Vector *) &v;
      this->std_map = (Epetra_BlockMap *) &v.Map();
      this->size = v.MyLength();
      this->owner = false;
    }

    template<typename Scalar>
    EpetraVector<Scalar>::~EpetraVector()
    {
      _F_;
      if (owner) free();
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::alloc(unsigned int n)
    {
      _F_;
      free();
      this->size = n;
      std_map = new Epetra_Map(this->size, 0, seq_comm); MEM_CHECK(std_map);
      vec = new Epetra_Vector(*std_map); MEM_CHECK(vec);
      vec_im = new Epetra_Vector(*std_map); MEM_CHECK(vec_im);
      zero();
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::zero()
    {
      _F_;
      for (unsigned int i = 0; i < this->size; i++) (*vec)[i] = 0.0;
      for (unsigned int i = 0; i < this->size; i++) (*vec_im)[i] = 0.0;
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::change_sign()
    {
      _F_;
      for (unsigned int i = 0; i < this->size; i++) (*vec)[i] *= -1.;
      for (unsigned int i = 0; i < this->size; i++) (*vec_im)[i] *= -1.;
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::free()
    {
      _F_;
      if(this->owner) 
      {
        delete std_map; std_map = NULL;
        delete vec; vec = NULL;
        delete vec_im; vec_im = NULL;
        this->size = 0;
      }
    }

    template<>
    void EpetraVector<double>::set(unsigned int idx, double y)
    {
      _F_;
      (*vec)[idx] = y;
    }

    template<>
    void EpetraVector<std::complex<double> >::set(unsigned int idx, std::complex<double> y)
    {
      _F_;
      (*vec)[idx] = std::real(y);
      (*vec_im)[idx] = std::imag(y);
    }

    template<>
    void EpetraVector<double>::add(unsigned int idx, double y)
    {
      _F_;
      (*vec)[idx] += y;
    }

    template<>
    void EpetraVector<std::complex<double> >::add(unsigned int idx, std::complex<double> y)
    {
      _F_;
      (*vec)[idx] += std::real(y);
      (*vec_im)[idx] += std::imag(y);
    }

    template<typename Scalar>
    void EpetraVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      _F_;
      for (unsigned int i = 0; i < n; i++)
        add(idx[i], y[i]);
    }

    template<typename Scalar>
    bool EpetraVector<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      return false;
    }

    template class HERMES_API EpetraMatrix<double>;
    template class HERMES_API EpetraMatrix<std::complex<double> >;
    template class HERMES_API EpetraVector<double>;
    template class HERMES_API EpetraVector<std::complex<double> >;
  }
}
#endif
