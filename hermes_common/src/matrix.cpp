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
#include "qsort.h"
#include "api.h"

void Hermes::Algebra::DenseMatrixOperations::ludcmp(double **a, int n, int *indx, double *d)
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

void Hermes::Algebra::DenseMatrixOperations::choldc(double **a, int n, double p[])
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
void Hermes::Algebra::Matrix<Scalar>::set_row_zero(unsigned int n)
{
  throw Hermes::Exceptions::MethodNotOverridenException("Matrix<Scalar>::set");
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::SparseMatrix()
{
  this->size = 0;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::SparseMatrix(unsigned int size)
{
  this->size = size;
  pages = NULL;

  row_storage = false;
  col_storage = false;
}

template<typename Scalar>
Hermes::Algebra::SparseMatrix<Scalar>::~SparseMatrix()
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
void Hermes::Algebra::SparseMatrix<Scalar>::prealloc(unsigned int n)
{
  this->size = n;

  pages = new Page *[n];
  memset(pages, 0, n * sizeof(Page *));
}

template<typename Scalar>
void Hermes::Algebra::SparseMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
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
int Hermes::Algebra::SparseMatrix<Scalar>::sort_and_store_indices(Page *page, int *buffer, int *max)
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
  for (int *p = buffer, last = -1; p < end; p++) if(*p != last) *q++= last = *p;

  return q - buffer;
}

template<typename Scalar>
int Hermes::Algebra::SparseMatrix<Scalar>::get_num_indices()
{
  int total = 0;
  for (unsigned int i = 0; i < this->size; i++)
    for (Page *page = pages[i]; page != NULL; page = page->next)
      total += page->count;

  return total;
}

static int find_position(int *Ai, int Alen, int idx)
{
  assert(Ai != NULL);
  assert(Alen > 0);
  assert(idx >= 0);

  register int lo = 0, hi = Alen - 1, mid;

  while (true)
  {
    mid = (lo + hi) >> 1;

    if(idx < Ai[mid]) hi = mid - 1;
    else if(idx > Ai[mid]) lo = mid + 1;
    else break;

    // Sparse matrix entry not found (raise an error when trying to add
    // value to this position, return 0 when obtaining value there).
    if(lo > hi)
    {
      mid = -1;
      break;
    }
  }
  return mid;
}




namespace Hermes
{

  namespace Algebra
  {
    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix()
    {
      this->size = 0; nnz = 0;
      Ap = NULL;
      Ai = NULL;
      Ax = NULL;
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix(unsigned int size)
    {
      this->size = size;
      this->alloc();
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::~CSCMatrix()
    {
      free();
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out)
    {
      int n = this->size;
      for (int j = 0; j<n; j++) vector_out[j] = 0;
      for (int j = 0; j<n; j++)
      {
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
        {
          vector_out[Ai[i]] += vector_in[j]*Ax[i];
        }
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      for (unsigned int i = 0; i < this->nnz; i++) Ax[i] *= value;
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::alloc()
    {
      assert(this->pages != NULL);

      // initialize the arrays Ap and Ai
      Ap = new int[this->size + 1];
      int aisize = this->get_num_indices();
      Ai = new int[aisize];

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

      Ax = new Scalar[nnz];
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::free()
    {
      nnz = 0;
      if(Ap != NULL)
      {
        delete [] Ap;
        Ap = NULL;
      }
      if(Ai != NULL)
      {
        delete [] Ai;
        Ai = NULL;
      }
      if(Ax != NULL)
      {
        delete [] Ax;
        Ax = NULL;
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::set_row_zero(unsigned int n)
    {
      for(int i = 0; i < Ap[n + 1] - Ap[n]; i++)
        Ax[Ap[n] + i] = Scalar(0);
    }

    template<typename Scalar>
    Scalar CSCMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      // Find m-th row in the n-th column.
      int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);

      if(mid < 0) // if the entry has not been found
        return 0.0;
      else
        return Ax[Ap[n] + mid];
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::zero()
    {
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<>
    void CSCMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSCMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp atomic
        Ax[Ap[n] + pos] += v;
      }
    }

    template<>
    void CSCMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSCMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp critical
        Ax[Ap[n] + pos] += v;
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_to_diagonal_blocks(int num_stages, CSCMatrix<Scalar>* mat_block)
    {
      int ndof = mat_block->get_size();
      if(this->get_size() != (unsigned int) num_stages * ndof)
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in CSCMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
      {
        this->add_as_block(ndof*i, ndof*i, mat_block);
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      add_to_diagonal_blocks(num_stages, static_cast<CSCMatrix<Scalar>*>(mat));
    }

    template<typename Scalar>
    unsigned int CSCMatrix<Scalar>::get_nnz() const
    {
      return this->nnz;
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, CSCMatrix<Scalar>* mat)
    {
      Hermes::Solvers::CSCIterator<Scalar> mat_it(mat);
      Hermes::Solvers::CSCIterator<Scalar> this_it(this);

      // Sanity check.
      bool this_not_empty = this_it.init();
      if(!this_not_empty) throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block().");

      // Iterate through the small matrix column by column and add all nonzeros
      // to the large one.
      bool mat_not_finished = mat_it.init();
      if(!mat_not_finished) throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block().");

      int mat_i, mat_j;
      Scalar mat_val;
      while(mat_not_finished)
      {
        mat_it.get_current_position(mat_i, mat_j, mat_val);
        bool found = this_it.move_to_position(mat_i + offset_i, mat_j + offset_j);
        if(!found)
          throw Hermes::Exceptions::Exception("Nonzero matrix entry at %d, %d not found in CSCMatrix<Scalar>::add_as_block().",
          mat_i + offset_i, mat_j + offset_j);
        this_it.add_to_current_position(mat_val);
        mat_not_finished = mat_it.move_ptr();
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_matrix(CSCMatrix<Scalar>* mat)
    {
      assert(this->get_size() == mat->get_size());
      // Create iterators for both matrices.
      Hermes::Solvers::CSCIterator<Scalar> mat_it(mat);
      Hermes::Solvers::CSCIterator<Scalar> this_it(this);
      int mat_i, mat_j;
      Scalar mat_val;
      int this_i, this_j;
      Scalar this_val;

      bool mat_not_finished = mat_it.init();
      bool this_not_finished = this_it.init();
      while(mat_not_finished && this_not_finished)
      {
        mat_it.get_current_position(mat_i, mat_j, mat_val);
        //printf("mat: current position %d %d %g\n", mat_i, mat_j, mat_val);
        this_it.get_current_position(this_i, this_j, this_val);
        //printf("this: current position %d %d %g\n", this_i, this_j, this_val);
        while(mat_i != this_i || mat_j != this_j)
        {
          //printf("SHOULD NOT BE HERE\n");
          this_not_finished = this_it.move_ptr();
          if(!this_not_finished)
          {
            printf("Entry %d %d does not exist in the matrix to which it is contributed.\n", mat_i, mat_j);
            throw Hermes::Exceptions::Exception("Incompatible matrices in add_umfpack_matrix().");
          }
          this_it.get_current_position(this_i, this_j, this_val);
        }
        this_it.add_to_current_position(mat_val);
        mat_not_finished = mat_it.move_ptr();
        this_not_finished = this_it.move_ptr();
        if(mat_not_finished && !this_not_finished)
          throw Hermes::Exceptions::Exception("Incompatible matrices in add_umfpack_matrix().");
      }
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i<this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void CSCMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      for (unsigned int i = 0; i < m; i++)       // rows
        for (unsigned int j = 0; j < n; j++)     // cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    double inline real(double x)
    {
      return x;
    }

    double inline imag(double x)
    {
      return 0;
    }

    double inline real(std::complex<double> x)
    {
      return x.real();
    }

    double inline imag(std::complex<double> x)
    {
      return x.imag();;
    }

    template<>
    bool CSCMatrix<double>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp =[\n",
          this->size, this->size, nnz, nnz);
        for (unsigned int j = 0; j < this->size; j++)
          for (int i = Ap[j]; i < Ap[j + 1]; i++)
          {
            fprintf(file, "%d %d ", Ai[i] + 1, j + 1);
            Hermes::Helpers::fprint_num(file, Ax[i], number_format);
            fprintf(file, "\n");
          }
          fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

          return true;

      case DF_MATRIX_MARKET:
        {
          fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real symmetric\n");
          int nnz_sym = 0;
          for (unsigned int j = 0; j < this->size; j++)
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
              if((int)j <= Ai[i]) nnz_sym++;
          fprintf(file, "%d %d %d\n", this->size, this->size, nnz_sym);
          for (unsigned int j = 0; j < this->size; j++)
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
              // The following line was replaced with the one below, because it gave a warning
                // to cause code abort at runtime.
                  //if(j <= Ai[i]) fprintf(file, "%d %d %24.15e\n", Ai[i] + 1, j + 1, Ax[i]);
                    if((int)j <= Ai[i])
                    {
                      fprintf(file, "%d %d ", Ai[i] + 1, (int)j + 1);
                      Hermes::Helpers::fprint_num(file, Ax[i], number_format);
                      fprintf(file, "\n");
                    }

                    return true;
        }

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESX\001", 1, 8, file);
          int ssize = sizeof(double);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(Ap, sizeof(int), this->size + 1, file);
          hermes_fwrite(Ai, sizeof(int), nnz, file);
          hermes_fwrite(Ax, sizeof(double), nnz, file);
          return true;
        }

      case DF_PLAIN_ASCII:
        exit(1);
        {
          const double zero_cutoff = 1e-10;
          double *ascii_entry_buff = new double[nnz];
          int *ascii_entry_i = new int[nnz];
          int *ascii_entry_j = new int[nnz];
          int k = 0;

          // If real or imaginary part of Scalar entry is below zero_cutoff
          // it's not included in ascii file, and number of non-zeros is reduced by one.
          for (unsigned int j = 0; j < size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              if(real(Ax[i]) > zero_cutoff || imag(Ax[i]) > zero_cutoff)
              {
                ascii_entry_buff[k] = Ax[i];
                ascii_entry_i[k] = Ai[i];
                ascii_entry_j[k] = j;
                k++;
              }
              else
                nnz -= 1;
            }
          }

          fprintf(file, "%d\n", size);
          fprintf(file, "%d\n", nnz);
          for (unsigned int k = 0; k < nnz; k++)
            fprintf(file, "%d %d %f\n", ascii_entry_i[k], ascii_entry_j[k], ascii_entry_buff[k]);

          //Free memory
          delete [] ascii_entry_buff;
          delete [] ascii_entry_i;
          delete [] ascii_entry_j;

          //Clear pointer
          ascii_entry_buff = NULL;
          ascii_entry_i = NULL;
          ascii_entry_j = NULL;

          return true;
        }

      default:
        return false;
      }
    }

    template<>
    bool CSCMatrix<std::complex<double> >::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp =[\n",
          this->size, this->size, nnz, nnz);
        for (unsigned int j = 0; j < this->size; j++)
          for (int i = Ap[j]; i < Ap[j + 1]; i++)
          {
            fprintf(file, "%d %d ", Ai[i] + 1, j + 1);
            Hermes::Helpers::fprint_num(file, Ax[i], number_format);
            fprintf(file, "\n");
          }
          fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

          return true;

      case DF_MATRIX_MARKET:
        {
          fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real symmetric\n");
          int nnz_sym = 0;
          for (unsigned int j = 0; j < this->size; j++)
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
              if((int)j <= Ai[i]) nnz_sym++;
          fprintf(file, "%d %d %d\n", this->size, this->size, nnz_sym);
          for (unsigned int j = 0; j < this->size; j++)
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
              // The following line was replaced with the one below, because it gave a warning
                // to cause code abort at runtime.
                  //if(j <= Ai[i]) fprintf(file, "%d %d %24.15e\n", Ai[i] + 1, j + 1, Ax[i]);
                    if((int)j <= Ai[i])
                    {
                      fprintf(file, "%d %d ", Ai[i] + 1, (int)j + 1);
                      Hermes::Helpers::fprint_num(file, Ax[i], number_format);
                      fprintf(file, "\n");
                    }

                    return true;
        }

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESX\001", 1, 8, file);
          int ssize = sizeof(std::complex<double>);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(Ap, sizeof(int), this->size + 1, file);
          hermes_fwrite(Ai, sizeof(int), nnz, file);
          hermes_fwrite(Ax, sizeof(std::complex<double>), nnz, file);
          return true;
        }

      case DF_PLAIN_ASCII:
        exit(1);
        {
          const double zero_cutoff = 1e-10;
          std::complex<double> *ascii_entry_buff = new std::complex<double>[nnz];
          int *ascii_entry_i = new int[nnz];
          int *ascii_entry_j = new int[nnz];
          int k = 0;

          // If real or imaginary part of Scalar entry is below zero_cutoff
          // it's not included in ascii file, and number of non-zeros is reduced by one.
          for (unsigned int j = 0; j < size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              if(real(Ax[i]) > zero_cutoff || imag(Ax[i]) > zero_cutoff)
              {
                ascii_entry_buff[k] = Ax[i];
                ascii_entry_i[k] = Ai[i];
                ascii_entry_j[k] = j;
                k++;
              }
              else
                nnz -= 1;
            }
          }

          fprintf(file, "%d\n", size);
          fprintf(file, "%d\n", nnz);
          for (unsigned int k = 0; k < nnz; k++)
            fprintf(file, "%d %d %E %E\n", ascii_entry_i[k], ascii_entry_j[k], ascii_entry_buff[k].real(), ascii_entry_buff[k].imag());

          //Free memory
          delete [] ascii_entry_buff;
          delete [] ascii_entry_i;
          delete [] ascii_entry_j;

          //Clear pointer
          ascii_entry_buff = NULL;
          ascii_entry_i = NULL;
          ascii_entry_j = NULL;

          return true;
        }

      default:
        return false;
      }
    }

    template<typename Scalar>
    unsigned int CSCMatrix<Scalar>::get_matrix_size() const
    {
      return this->size;
    }

    template<typename Scalar>
    double CSCMatrix<Scalar>::get_fill_in() const
    {
      return nnz / (double) (this->size * this->size);
    }

    template<typename Scalar>
    void CSCMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
    {
      this->nnz = nnz;
      this->size = size;
      this->Ap = new int[this->size + 1]; assert(this->Ap != NULL);
      this->Ai = new int[nnz];    assert(this->Ai != NULL);
      this->Ax = new Scalar[nnz]; assert(this->Ax != NULL);
      for (unsigned int i = 0; i < this->size + 1; i++) this->Ap[i] = ap[i];
      for (unsigned int i = 0; i < nnz; i++)
      {
        this->Ax[i] = ax[i];
        this->Ai[i] = ai[i];
      }
    }

    template<typename Scalar>
    CSCMatrix<Scalar>* CSCMatrix<Scalar>::duplicate()
    {
      CSCMatrix<Scalar>* new_matrix = new CSCMatrix<Scalar>();
      new_matrix->create(this->get_size(), this->get_nnz(), this->get_Ap(),  this->get_Ai(),  this->get_Ax());
      return new_matrix;
    }

    template<typename Scalar>
    int *CSCMatrix<Scalar>::get_Ap()
    {
      return this->Ap;
    }

    template<typename Scalar>
    int *CSCMatrix<Scalar>::get_Ai()
    {
      return this->Ai;
    }

    template<typename Scalar>
    Scalar *CSCMatrix<Scalar>::get_Ax()
    {
      return this->Ax;
    }
  }

  namespace Solvers
  {
    template<typename Scalar>
    bool CSCIterator<Scalar>::init()
    {
      if(this->size == 0 || this->nnz == 0) return false;
      this->Ap_pos = 0;
      this->Ai_pos = 0;
      return true;
    }

    template<typename Scalar>
    CSCIterator<Scalar>::CSCIterator(CSCMatrix<Scalar>* mat)
    {
      this->size = mat->get_size();
      this->nnz = mat->get_nnz();
      this->Ai = mat->get_Ai();
      this->Ap = mat->get_Ap();
      this->Ax = mat->get_Ax();
      this->Ai_pos = 0;
      this->Ap_pos = 0;
    }

    template<typename Scalar>
    void CSCIterator<Scalar>::get_current_position(int& i, int& j, Scalar& val)
    {
      i = Ai[Ai_pos];
      j = Ap_pos;
      val = Ax[Ai_pos];
    }

    template<typename Scalar>
    bool CSCIterator<Scalar>::move_to_position(int i, int j)
    {
      int ii, jj;
      Scalar val;
      get_current_position(ii, jj, val);
      while (!(ii == i && jj == j))
      {
        if(!this->move_ptr()) return false;
        get_current_position(ii, jj, val);
      }
      return true;
    }

    template<typename Scalar>
    bool CSCIterator<Scalar>::move_ptr()
    {
      if(Ai_pos >= nnz - 1) return false; // It is no longer possible to find next element.
      if(Ai_pos + 1 >= Ap[Ap_pos + 1])
      {
        Ap_pos++;
      }
      Ai_pos++;
      return true;
    }

    template<typename Scalar>
    void CSCIterator<Scalar>::add_to_current_position(Scalar val)
    {
      this->Ax[this->Ai_pos] += val;
    }
  }
}
template<typename Scalar>
SparseMatrix<Scalar>* Hermes::Algebra::create_matrix()
{
  switch (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
  {
  case Hermes::SOLVER_AMESOS:
    {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
      return new EpetraMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_AZTECOO:
    {
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
      return new EpetraMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_MUMPS:
    {
#ifdef WITH_MUMPS
      return new MumpsMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_PETSC:
    {
#ifdef WITH_PETSC
      return new PetscMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_UMFPACK:
    {
#ifdef WITH_UMFPACK
      return new UMFPackMatrix<Scalar>;
#else
      throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_SUPERLU:
    {
#ifdef WITH_SUPERLU
      return new SuperLUMatrix<Scalar>;
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

template<typename Scalar>
Vector<Scalar>* Hermes::Algebra::create_vector()
{
  switch (Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
  {
  case Hermes::SOLVER_AMESOS:
    {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
      return new EpetraVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_AZTECOO:
    {
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
      return new EpetraVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_MUMPS:
    {
#ifdef WITH_MUMPS
      return new MumpsVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_PETSC:
    {
#ifdef WITH_PETSC
      return new PetscVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_UMFPACK:
    {
#ifdef WITH_UMFPACK
      return new UMFPackVector<Scalar>;
#else
      throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
      break;
    }
  case Hermes::SOLVER_SUPERLU:
    {
#ifdef WITH_SUPERLU
      return new SuperLUVector<Scalar>;
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

template class Hermes::Algebra::SparseMatrix<double>;
template class Hermes::Algebra::SparseMatrix<std::complex<double> >;

template class HERMES_API Hermes::Algebra::CSCMatrix<double>;
template class HERMES_API Hermes::Algebra::CSCMatrix<std::complex<double> >;

template HERMES_API Vector<double>* Hermes::Algebra::create_vector();
template HERMES_API SparseMatrix<double>*  Hermes::Algebra::create_matrix();

template HERMES_API Vector<std::complex<double> >* Hermes::Algebra::create_vector();
template HERMES_API SparseMatrix<std::complex<double> >*  Hermes::Algebra::create_matrix();

template class HERMES_API Hermes::Solvers::CSCIterator<double>;
template class HERMES_API Hermes::Solvers::CSCIterator<std::complex<double> >;