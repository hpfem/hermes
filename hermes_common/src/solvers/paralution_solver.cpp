// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file paralution_solver.cpp
\brief PARALUTION solver interface.
*/
#include "config.h"
#ifdef WITH_PARALUTION
#include "paralution_solver.h"

namespace Hermes
{
  namespace Algebra
  {
    template<typename Scalar>
    ParalutionMatrix<Scalar>::ParalutionMatrix(ParalutionMatrixType type) 
      : paralutionMatrixType(type)
    {
      this->paralutionMatrix = new paralution::LocalMatrix();
    }

    template<typename Scalar>
    ParalutionMatrix<Scalar>::~ParalutionMatrix()
    {
      delete this->paralutionMatrix;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::alloc()
    {
      assert(this->pages != NULL);

      // initialize the arrays Ap and Ai
      int* Ap = new int[this->size + 1];
      int aisize = this->get_num_indices();
      int* Ai = new int[aisize];

      // sort the indices and remove duplicities, insert into Ai
      unsigned int i;
      int pos = 0;
      for (i = 0; i < this->size; i++)
      {
        Ap[i] = pos;
        pos += this->sort_and_store_indices(this->pages[i], Ai + pos, Ai + aisize);
      }
      Ap[i] = pos;

      delete [] this->pages;
      this->pages = NULL;

      nnz = Ap[this->size];

      switch(this->paralutionMatrixType)
      {
      case ParalutionMatrixTypeCSR :
        this->paralutionMatrix->AllocateCSR("Paralution matrix", nnz, size, size);
        break;
      case ParalutionMatrixTypeMCSR :
        this->paralutionMatrix->AllocateMCSR("Paralution matrix", nnz, size, size);
        break;
      case ParalutionMatrixTypeCOO :
        this->paralutionMatrix->AllocateCOO("Paralution matrix", nnz, size, size);
        break;
      case ParalutionMatrixTypeDENSE :
        this->paralutionMatrix->AllocateDENSE("Paralution matrix", nnz, size, size);
        break;
      }

      delete [] Ai;
      delete [] Ap;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::free()
    {
      this->paralutionMatrix->Clear();
    }

    template<typename Scalar>
    Scalar ParalutionMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      Scalar toReturn;
      int row, column;
      paralution::LocalMatrix<Scalar> localMatrix;
      localMatrix.AllocateDENSE("temporary", 1, 1);
      return this->paralutionMatrix->ExtractSubMatrix(m, n, 1, 1, &localMatrix);
      localMatrix.CopyToCSR(&row, &col, &toReturn);
      return toReturn;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::zero()
    {
      this->paralutionMatrix->Zeros();
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar v)
    {
      this->paralutionMatrix->SetDataPtrCOO(&&n, &&m, &&v, "temporary", 1, 1, 1);
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::add_to_diagonal(Scalar v)
    {
      this->paralutionMatrix->AddScalarDiagonal(v);
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      // not implemented yet.
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      this->paralutionMatrix->SetDataPtrCOO(&rows, &cols, &&m, mat, "temporary", 1, 1, 1);
    }
    template<typename Scalar>
    bool ParalutionMatrix<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf")
    {
    }

    template<typename Scalar>
    unsigned int ParalutionMatrix<Scalar>::get_matrix_size() const
    {
    }
    template<typename Scalar>
    unsigned int ParalutionMatrix<Scalar>::get_nnz() const
    {
    }
    template<typename Scalar>
    double ParalutionMatrix<Scalar>::get_fill_in() const
    {
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out)
    {
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector()
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::~ParalutionVector()
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::alloc(unsigned int ndofs)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::free()
    {
    }

    template<typename Scalar>
    Scalar ParalutionVector<Scalar>::get(unsigned int idx)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::extract(Scalar *v) const
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::zero()
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::change_sign()
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::set(unsigned int idx, Scalar y)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add(unsigned int idx, Scalar y)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Scalar* vec)
    {
    }

    template<typename Scalar>
    bool ParalutionVector<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf")
    {
    }
  }
  namespace Solvers
  {
    template <typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::ParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *m, ParalutionVector<Scalar> *rhs)
    {
    }

    template <typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::~ParalutionLinearMatrixSolver()
    {
    }

    template <typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve()
    {
    }

    template <typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
    }

    template <typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
    }

    template <typename Scalar>
    double ParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
    }

    template <typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
    }
  }
}
#endif
