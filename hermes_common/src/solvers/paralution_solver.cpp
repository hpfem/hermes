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
      : CSRMatrix<Scalar>(), paralutionMatrixType(type)
    {
    }

    template<typename Scalar>
    ParalutionMatrix<Scalar>::~ParalutionMatrix()
    {
      this->paralutionMatrix.Clear();
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::free()
    {
      CSRMatrix<Scalar>::free();
      this->paralutionMatrix.Clear();
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::zero()
    {
      CSRMatrix<Scalar>::zero();
      this->paralutionMatrix.Zeros();
    }

    template<typename Scalar>
    paralution::LocalMatrix<Scalar>& ParalutionMatrix<Scalar>::get_paralutionMatrix()
    {
      return this->paralutionMatrix;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::alloc()
    {
      assert(this->pages != NULL);

      // initialize the arrays Ap and Ai
      this->Ap = new int[this->size + 1];
      int aisize = this->get_num_indices();
      this->Ai = new int[aisize];

      // sort the indices and remove duplicities, insert into Ai
      unsigned int i;
      int pos = 0;
      for (i = 0; i < this->size; i++)
      {
        this->Ap[i] = pos;
        pos += this->sort_and_store_indices(this->pages[i], this->Ai + pos, this->Ai + aisize);
      }
      this->Ap[i] = pos;

      delete [] this->pages;
      this->pages = NULL;

      this->nnz = this->Ap[this->size];

      this->Ax = new Scalar[this->nnz];
      memset(this->Ax, 0, sizeof(Scalar) * this->nnz);

      this->paralutionMatrix.SetDataPtrCSR(&this->Ap, &this->Ai, &this->Ax, "paralutionMatrix", this->nnz, this->size, this->size);
    }

    
    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector() : Vector<Scalar>(), v(NULL)
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector(unsigned int size) : Vector<Scalar>(size), v(NULL)
    {
      this->alloc(size);
      this->paralutionVector.SetDataPtr(&this->v, "paralutionVector", this->size);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      v = new Scalar[n];
      this->zero();
      this->paralutionVector.SetDataPtr(&this->v, "vector", this->size);
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::~ParalutionVector()
    {
      this->paralutionVector.Clear();
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::free()
    {
      delete [] v;
      v = NULL;
      this->size = 0;
      this->paralutionVector.Clear();
    }

    template<typename Scalar>
    paralution::LocalVector<Scalar>& ParalutionVector<Scalar>::get_paralutionVector()
    {
      return this->paralutionVector;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::zero()
    {
      memset(v, 0, this->size * sizeof(Scalar));
      this->paralutionVector.Zeros();
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++) v[i] *= -1.;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      v[idx] = y;
    }

    template<>
    void ParalutionVector<double>::add(unsigned int idx, double y)
    {
#pragma omp atomic
      v[idx] += y;
    }


    template<typename Scalar>
    void ParalutionVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        v[idx[i]] += y[i];
    }

    template<typename Scalar>
    Scalar ParalutionVector<Scalar>::get(unsigned int idx)
    {
      return v[idx];
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::extract(Scalar *v) const
    {
      memcpy(v, this->v, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->v[i] += vec->get(i);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->length(); i++)
        this->v[i] += vec[i];
    }

    template<>
    bool ParalutionVector<double>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx1\n%s =[\n", this->size, var_name);
        for (unsigned int i = 0; i < this->size; i++)
        {
          Hermes::Helpers::fprint_num(file, v[i], number_format);
          fprintf(file, "\n");
        }
        fprintf(file, " ];\n");
        return true;

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESR\001", 1, 8, file);
          int ssize = sizeof(double);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(v, sizeof(double), this->size, file);
          return true;
        }

      case DF_PLAIN_ASCII:
        {
          fprintf(file, "\n");
          for (unsigned int i = 0; i < size; i++)
          {
            Hermes::Helpers::fprint_num(file, v[i], number_format);
            fprintf(file, "\n");
          }

          return true;
        }

      default:
        return false;
      }
    }

    template class HERMES_API ParalutionMatrix<double>;
    template class HERMES_API ParalutionVector<double>;
  }
  namespace Solvers
  {
    template<typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::ParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : IterSolver<Scalar>(), matrix(matrix), rhs(rhs)
    {
    }

    template<typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::~ParalutionLinearMatrixSolver()
    {
    }

    template<typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve()
    {
      paralution::init_paralution();
      paralution::set_omp_threads_paralution(HermesCommonApi.get_integral_param_value(numThreads));
      paralution::info_paralution();

      this->matrix->get_paralutionMatrix().MoveToAccelerator();
      rhs->get_paralutionVector().MoveToAccelerator();

      assert(matrix != NULL);
      assert(rhs != NULL);
      assert(matrix->get_size() == rhs->length());

      this->paralutionSolver.SetOperator(this->matrix->get_paralutionMatrix());
      this->paralutionSolver.SetPreconditioner(this->paralutionPreconditioner);

      this->paralutionSolver.Build();

      paralution::LocalVector<Scalar> x;
      x.Allocate("x", matrix->get_size());
      x.MoveToAccelerator();

      this->paralutionSolver.Solve(rhs->get_paralutionVector(), &x);

      x.LeaveDataPtr(&this->sln);

      paralution::stop_paralution();

      x.Clear();
      return true;
    }

    template<typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return matrix->get_size();
    }

    template<typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
      return 0;
    }

    template<typename Scalar>
    double ParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
      return 0.;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
    }

    template class HERMES_API ParalutionLinearMatrixSolver<double>;
  }
}
#endif
