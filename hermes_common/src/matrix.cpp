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
/*! \file matrix.cpp
\brief Basic matrix classes and operations.
*/
#include "common.h"
#include "matrix.h"
#include "callstack.h"

#include "solvers/linear_matrix_solver.h"
#include "solvers/umfpack_solver.h"
#include "solvers/superlu_solver.h"
#include "solvers/amesos_solver.h"
#include "solvers/petsc_solver.h"
#include "solvers/mumps_solver.h"
#include "solvers/aztecoo_solver.h"
#include "solvers/paralution_solver.h"
#include "qsort.h"
#include "api.h"

namespace Hermes
{
  namespace Algebra
  {
    void DenseMatrixOperations::ludcmp(double **a, int n, int *indx, double *d)
    {
      int i, imax = 0, j, k;
      double big, dum, sum, temp;
      double *vv = new double[n];

      *d = 1.0;
      for (i = 0; i < n; i++)
      {
        big = 0.0;
        for (j = 0; j < n; j++)
          if((temp = fabs(a[i][j])) > big)
            big = temp;
        if(big == 0.0)
        {
          delete [] vv;
          throw Exceptions::Exception("Singular matrix in routine LUDCMP!");
        }
        vv[i] = 1.0 / big;
      }
      for (j = 0; j < n; j++)
      {
        for (i = 0; i < j; i++)
        {
          sum = a[i][j];
          for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
          a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++)
        {
          sum = a[i][j];
          for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
          a[i][j] = sum;
          if((dum = vv[i]*fabs(sum)) >= big)
          {
            big = dum;
            imax = i;
          }
        }
        if(j != imax)
        {
          for (k = 0; k < n; k++)
          {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
          }
          *d = -(*d);
          vv[imax] = vv[j];
        }
        indx[j] = imax;
        if(a[j][j] == 0.0) a[j][j] = 1.0e-20;
        if(j != n-1)
        {
          dum = 1.0 / (a[j][j]);
          for (i = j + 1; i < n; i++) a[i][j] *= dum;
        }
      }
      delete [] vv;
    }

    void DenseMatrixOperations::choldc(double **a, int n, double p[])
    {
      int i, j, k;
      for (i = 0; i < n; i++)
      {
        for (j = i; j < n; j++)
        {
          double sum = a[i][j];
          k = i;
          while (--k >= 0) sum -= a[i][k] * a[j][k];
          if(i == j)
          {
            if(sum <= 0.0)
              throw Exceptions::Exception("CHOLDC failed!");
            else p[i] = sqrt(sum);
          }
          else a[j][i] = sum / p[i];
        }
      }
    }

    template<typename Scalar>
    Matrix<Scalar>::Matrix(unsigned int size)
    {
      this->size = size;
    }

    template<typename Scalar>
    void Matrix<Scalar>::set_row_zero(unsigned int n)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("Matrix<Scalar>::set_row_zero");
    }

    template<typename Scalar>
    void Matrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i < this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void Matrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const
    {
      if(!vector_out_initialized)
        vector_out = new Scalar[this->size];
      for(int i = 0; i < this->size; i++)
      {
        vector_out[i] = Scalar(0.);
        for(int j = 0; j < this->size; j++)
          vector_out[i] += this->get(i, j) * vector_in[j];
      }
    }

    template<typename Scalar>
    void Matrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      throw Hermes::Exceptions::MethodNotOverridenException("multiply_with_Scalar()");
    }

    template<typename Scalar>
    unsigned int Matrix<Scalar>::get_size() const
    {
      return this->size;
    };

    template<typename Scalar>
    void Matrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      for (unsigned int i = 0; i < m; i++)       // rows
        for (unsigned int j = 0; j < n; j++)     // cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix()
    {
      this->size = 0;
      pages = NULL;

      row_storage = false;
      col_storage = false;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
    {
      this->size = size;
      pages = NULL;

      row_storage = false;
      col_storage = false;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::SparseMatrix(const SparseMatrix<Scalar>& mat)
    {
      this->size = mat.get_size();

      if (mat.pages)
      {
        this->pages = new Page *[this->size];
        memset(this->pages, 0, this->size * sizeof(Page *));

        for (int col = 0; col < this->size; col++)
        {
          Page *page = mat.pages[col];
          Page *new_page = new Page;
          new_page->count = page->count;
          memcpy(new_page->idx, page->idx, sizeof(int) * page->count);
          new_page->next = this->pages[col];
          this->pages[col] = new_page; 
        }
      }
      else
      {
        this->pages = NULL;
      }

      row_storage = false;
      col_storage = false;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>::~SparseMatrix()
    {
      if(pages)
      {
        for (unsigned int i = 0; i < this->size; i++)
          if(pages[i])
            delete pages[i];
        delete [] pages;
      }
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::finish()
    {
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_row_entries(unsigned int row) const
    {
      return -1;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::extract_row_copy(unsigned int row, unsigned int len,
      unsigned int &n_entries, double *vals,
      unsigned int *idxs) const
    {
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_col_entries(unsigned int col) const
    {
      return -1;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::extract_col_copy(unsigned int col, unsigned int len,
      unsigned int &n_entries, double *vals,
      unsigned int *idxs) const
    {
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_sparse_matrix(SparseMatrix<Scalar>* mat)
    {
      add_as_block(0, 0, mat);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      int ndof = mat->get_size();
      if(this->get_size() != (unsigned int) num_stages * ndof)
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in SparseMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
        this->add_as_block(ndof*i, ndof*i, mat);
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, SparseMatrix<Scalar>* mat)
    {
      if((this->get_size() < offset_i + mat->get_size() )||(this->get_size() < offset_j + mat->get_size() ))
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in SparseMatrix<Scalar>::add_as_block()");
      unsigned int block_size = mat->get_size();
      for (unsigned int r = 0; r < block_size; r++)
      {
        for (unsigned int c = 0; c < block_size; c++)
        {
          this->add(offset_i + r, offset_j + c, mat->get(r, c));
        }
      }
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* SparseMatrix<Scalar>::duplicate() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("SparseMatrix* duplicate()");
    }

    template<typename Scalar>
    unsigned int SparseMatrix<Scalar>::get_nnz() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("get_nnz()");
      return 0;
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::prealloc(unsigned int n)
    {
      this->size = n;

      pages = new Page *[n];
      memset(pages, 0, n * sizeof(Page *));
    }

    template<typename Scalar>
    void SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      if(pages[col] == NULL || pages[col]->count >= PAGE_SIZE)
      {
        Page *new_page = new Page;
        new_page->count = 0;
        new_page->next = pages[col];
        pages[col] = new_page;
      }
      pages[col]->idx[pages[col]->count++] = row;
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
    {
      // gather all pages in the buffer, deleting them along the way
      int *end = buffer;
      while (page != NULL)
      {
        memcpy(end, page->idx, sizeof(int) * page->count);
        end += page->count;
        Page *tmp = page;
        page = page->next;
        delete tmp;
      }

      // sort the indices and remove duplicities
      qsort_int(buffer, end - buffer);
      int *q = buffer;
      for (int *p = buffer, last = -1; p < end; p++)
        if(*p != last)
          *q++ = last = *p;

      return q - buffer;
    }

    template<typename Scalar>
    int SparseMatrix<Scalar>::get_num_indices()
    {
      int total = 0;
      for (unsigned int i = 0; i < this->size; i++)
        for (Page *page = pages[i]; page != NULL; page = page->next)
          total += page->count;

      return total;
    }


    template<typename Scalar>
    Vector<Scalar>::Vector() : size(0)
    {
    }

    template<typename Scalar>
    Vector<Scalar>::Vector(unsigned int size) : size(size)
    {
    }

    template<typename Scalar>
    void Vector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++) this->set(i, vec->get(i));
    }

    template<typename Scalar>
    void Vector<Scalar>::set_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++) this->set(i, vec[i]);
    }

    template<typename Scalar>
    void Vector<Scalar>::add_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++) this->add(i, vec->get(i));
    }

    template<typename Scalar>
    void Vector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++) this->add(i, vec[i]);
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::export_to_file(char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      if(!v)
        throw Exceptions::MethodNotOverridenException("Vector<Scalar>::export_to_file");

      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real\n");
          else
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate complex\n");

          fprintf(file, "%d 1 %d\n", this->size, this->size);

          for (unsigned int j = 0; j < this->size; j++)
          {
            Hermes::Helpers::fprint_coordinate_num(file, j + 1, 1, v[j], number_format);
            fprintf(file, "\n");
          }

          fclose(file);
        }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
        {
#ifdef WITH_MATIO
          size_t dims[2];
          dims[0] = this->size;
          dims[1] = 1;

          mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);
          matvar_t *matvar;

          // For complex.
          double* v_re = NULL;
          double* v_im = NULL;

          void* data;
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
          {
            data = v;
            matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA);
          }
          else
          {
            v_re = new double[this->size];
            v_im = new double[this->size];
            struct mat_complex_split_t z = {v_re, v_im};

            for(int i = 0; i < this->size; i++)
            {
              v_re[i] = ((std::complex<double>)(this->v[i])).real();
              v_im[i] = ((std::complex<double>)(this->v[i])).imag();
              data = &z;
            }
            matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
          } 

          if (matvar)
          {
            Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
          }

          if(v_re)
            delete [] v_re;
          if(v_im)
            delete [] v_im;
          Mat_Close(mat);

          if(!matvar)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#else
          throw Exceptions::Exception("MATIO not included.");
#endif
        }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          for (unsigned int i = 0; i < this->size; i++)
          {
            Hermes::Helpers::fprint_num(file, v[i], number_format);
            fprintf(file, "\n");
          }
          fclose(file);
        }
      }
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::import_from_file(char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          std::vector<Scalar> data;
          std::ifstream input (filename);
          if(input.bad())
            throw Exceptions::IOException(Exceptions::IOException::Read, filename);
          std::string lineData;

          while(getline(input, lineData))
          {
            Scalar d;
            std::stringstream lineStream(lineData);
            lineStream >> d;
            data.push_back(d);
          }

          this->alloc(data.size());
          memcpy(this->v, &data[0], sizeof(Scalar)*data.size());
        }
        break;
      case EXPORT_FORMAT_MATLAB_MATIO:
#ifdef WITH_MATIO
        mat_t    *matfp;
        matvar_t *matvar;

        matfp = Mat_Open(filename,MAT_ACC_RDONLY);

        if (!matfp )
        {
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
          return;
        }

        matvar = Mat_VarRead(matfp, var_name);
        if (matvar)
        {
          this->alloc(matvar->dims[0]);
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
            memcpy(this->v, matvar->data, sizeof(Scalar)*this->size);
          else
          {
            std::complex<double>* complex_data = new std::complex<double>[this->size];
            double* real_array = (double*)((mat_complex_split_t*)matvar->data)->Re;
            double* imag_array = (double*)((mat_complex_split_t*)matvar->data)->Im;
            for(int i = 0; i < this->size; i++)
              complex_data[i] = std::complex<double>(real_array[i], imag_array[i]);
            memcpy(this->v, complex_data, sizeof(Scalar)*this->size);
            delete [] complex_data;
          }
        }

        Mat_Close(matfp);
        if(!matvar)
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
#else
        throw Exceptions::Exception("MATIO not included.");
#endif
        break;
      case EXPORT_FORMAT_MATRIX_MARKET:
        throw Hermes::Exceptions::MethodNotImplementedException("SimpleVector<Scalar>::import_from_file - Matrix Market");
      }

    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector() : Vector<Scalar>(), v(NULL)
    {
    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector(unsigned int size) : Vector<Scalar>(size), v(NULL)
    {
    }

    template<typename Scalar>
    SimpleVector<Scalar>::~SimpleVector()
    {
      free();
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      this->v = new Scalar[n];
      zero();
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++)
        v[i] *= -1.;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::zero()
    {
      memset(this->v, 0, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::free()
    {
      if (this->v)
        delete [] this->v;
      this->v = NULL;
      this->size = 0;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      this->v[idx] = y;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add(unsigned int idx, Scalar y)
    {
#pragma omp critical (SimpleVector_add)
      this->v[idx] += y;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        this->v[idx[i]] += y[i];
    }

    template<typename Scalar>
    Scalar SimpleVector<Scalar>::get(unsigned int idx) const
    {
      return this->v[idx];
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::extract(Scalar *v) const
    {
      memcpy(v, this->v, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i));
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i]);
    }

    template HERMES_API void Vector<double>::set_vector(Vector<double>* vec);
    template HERMES_API void Vector<double>::set_vector(double* vec);
    template HERMES_API void Vector<double>::add_vector(Vector<double>* vec);
    template HERMES_API void Vector<double>::add_vector(double* vec);

    template HERMES_API void Vector<std::complex<double> >::set_vector(Vector<std::complex<double> >* vec);
    template HERMES_API void Vector<std::complex<double> >::set_vector(std::complex<double> * vec);
    template HERMES_API void Vector<std::complex<double> >::add_vector(Vector<std::complex<double> >* vec);
    template HERMES_API void Vector<std::complex<double> >::add_vector(std::complex<double> * vec);

    template<>
    HERMES_API SparseMatrix<double>* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new CSCMatrix<double>;
        }

      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new CSCMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new CSCMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return NULL;
    }

    template<>
    HERMES_API Vector<double>* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new SimpleVector<double>;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionVector<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return NULL;
    }

    template<>
    HERMES_API SparseMatrix<std::complex<double> >* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new CSCMatrix<std::complex<double> >;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new CSCMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new CSCMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return NULL;
    }

    template<>
    HERMES_API Vector<std::complex<double> >* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new SimpleVector<std::complex<double> >;
        }

      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");

#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return NULL;
    }

    template class Matrix<double>;
    template class Matrix<std::complex<double> >;

    template class SparseMatrix<double>;
    template class SparseMatrix<std::complex<double> >;

    template class Vector<double>;
    template class Vector<std::complex<double> >;

    template class SimpleVector<double>;
    template class SimpleVector<std::complex<double> >;
  }
  namespace Mixins
  {
    template<typename Scalar>
    MatrixRhsOutput<Scalar>::MatrixRhsOutput() : output_matrixOn(false), output_matrixIterations(-1), matrixFilename("Matrix_"),
      matrixVarname("A"), matrixFormat(Hermes::Algebra::EXPORT_FORMAT_PLAIN_ASCII), matrix_number_format("%lf"), output_rhsOn(false), output_rhsIterations(-1),
      RhsFilename("Rhs_"), RhsVarname("b"), RhsFormat(Hermes::Algebra::EXPORT_FORMAT_PLAIN_ASCII), rhs_number_format("%lf")
    {
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix, int iteration)
    {
      if (matrix == NULL)
        return;

      char* fileName = new char[this->matrixFilename.length() + 5];

      if(this->output_matrixOn)
      {
        if(this->only_lastMatrixIteration)
          sprintf(fileName, "%s", this->matrixFilename.c_str());
        else if(this->output_matrixIterations == -1 || this->output_matrixIterations >= iteration)
          sprintf(fileName, "%s%i", this->matrixFilename.c_str(), iteration);
        else
        {
          delete [] fileName;
          return;
        }

        matrix->export_to_file(fileName, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
        delete [] fileName;
      }
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix)
    {
      if (matrix == NULL)
        return;

      if(this->output_matrixOn)
      {
        char* fileName = new char[this->matrixFilename.length() + 5];
        sprintf(fileName, "%s", this->matrixFilename.c_str());
        matrix->export_to_file(fileName, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
        delete [] fileName;
      }
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs, int iteration)
    {
      if (rhs == NULL)
        return;

      char* fileName = new char[this->RhsFilename.length() + 5];

      if(this->output_rhsOn)
      {
        if(this->only_lastRhsIteration)
          sprintf(fileName, "%s", this->RhsFilename.c_str());
        else if(this->output_rhsIterations == -1 || this->output_rhsIterations >= iteration)
          sprintf(fileName, "%s%i", this->RhsFilename.c_str(), iteration);
        else
        {
          delete [] fileName;
          return;
        }

        rhs->export_to_file(fileName, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
        delete [] fileName;
      }
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs)
    {
      if (rhs == NULL)
        return;

      if(this->output_rhsOn)
      {
        char* fileName = new char[this->RhsFilename.length() + 5];
        sprintf(fileName, "%s", this->RhsFilename.c_str());
        rhs->export_to_file(fileName, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
        delete [] fileName;
      }
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::output_matrix(bool only_last_iteration, int firstIterations)
    {
      output_matrixOn = true;
      this->only_lastMatrixIteration = only_last_iteration;
      this->output_matrixIterations = firstIterations;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_matrix_filename(std::string name)
    {
      this->matrixFilename = name;
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_print_zero_matrix_entries(bool to_set)
    {
      this->print_matrix_zero_values = to_set;
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_matrix_varname(std::string name)
    {
      this->matrixVarname = name;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_matrix_export_format(Hermes::Algebra::MatrixExportFormat format)
    {
      this->matrixFormat = format;
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_matrix_number_format(char* number_format)
    {
      this->matrix_number_format = number_format;
    }

    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::output_rhs(bool only_last_iteration, int firstIterations)
    {
      this->output_rhsOn = true;
      this->only_lastRhsIteration = only_last_iteration;
      this->output_rhsIterations = firstIterations;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_rhs_filename(std::string name)
    {
      this->RhsFilename = name;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_rhs_varname(std::string name)
    {
      this->RhsVarname = name;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_rhs_export_format(Hermes::Algebra::MatrixExportFormat format)
    {
      this->RhsFormat = format;
    }
    template<typename Scalar>
    void MatrixRhsOutput<Scalar>::set_rhs_number_format(char* number_format)
    {
      this->rhs_number_format = number_format;
    }

    template HERMES_API class MatrixRhsOutput<double>;
    template HERMES_API class MatrixRhsOutput<std::complex<double> >;
  }
}
