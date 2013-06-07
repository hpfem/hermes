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
/*! \file cs_matrix.cpp
\brief Basic cs (Compressed sparse) matrix classes and operations.
*/
#include "cs_matrix.h"

namespace Hermes
{
  namespace Algebra
  {
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

    template<typename Scalar>
    CSMatrix<Scalar>::CSMatrix()
    {
      this->size = 0; nnz = 0;
      Ap = NULL;
      Ai = NULL;
      Ax = NULL;
    }

    template<typename Scalar>
    CSMatrix<Scalar>::CSMatrix(unsigned int size)
    {
      this->size = size;
      this->alloc();
    }

    template<typename Scalar>
    CSMatrix<Scalar>::~CSMatrix()
    {
      free();
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out)
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
    void CSMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      for (unsigned int i = 0; i < this->nnz; i++) Ax[i] *= value;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::alloc()
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
    void CSMatrix<Scalar>::free()
    {
      nnz = 0;
      if(Ap)
      {
        delete [] Ap;
        Ap = NULL;
      }
      if(Ai)
      {
        delete [] Ai;
        Ai = NULL;
      }
      if(Ax)
      {
        delete [] Ax;
        Ax = NULL;
      }
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::set_row_zero(unsigned int n)
    {
      for(int i = 0; i < Ap[n + 1] - Ap[n]; i++)
        Ax[Ap[n] + i] = Scalar(0);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::zero()
    {
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::add_to_diagonal_blocks(int num_stages, CSMatrix<Scalar>* mat_block)
    {
      int ndof = mat_block->get_size();
      if(this->get_size() != (unsigned int) num_stages * ndof)
        throw Hermes::Exceptions::Exception("Incompatible matrix sizes in CSMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
      {
        this->add_as_block(ndof*i, ndof*i, mat_block);
      }
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      add_to_diagonal_blocks(num_stages, static_cast<CSMatrix<Scalar>*>(mat));
    }

    template<typename Scalar>
    unsigned int CSMatrix<Scalar>::get_nnz() const
    {
      return this->nnz;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i<this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void CSMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
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

    template<typename Scalar>
    unsigned int CSMatrix<Scalar>::get_matrix_size() const
    {
      return this->size;
    }

    template<typename Scalar>
    double CSMatrix<Scalar>::get_fill_in() const
    {
      return nnz / (double) (this->size * this->size);
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
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
    int *CSMatrix<Scalar>::get_Ap()
    {
      return this->Ap;
    }

    template<typename Scalar>
    int *CSMatrix<Scalar>::get_Ai()
    {
      return this->Ai;
    }

    template<typename Scalar>
    Scalar *CSMatrix<Scalar>::get_Ax()
    {
      return this->Ax;
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, CSMatrix<Scalar>* mat)
    {
      throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::add_as_block");
    }

    template<typename Scalar>
    void CSMatrix<Scalar>::add_matrix(CSMatrix<Scalar>* mat)
    {
      throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::add_matrix");
    }

    template<typename Scalar>
    CSMatrix<Scalar>* CSMatrix<Scalar>::duplicate()
    {
      throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::duplicate()");
      return NULL;
    }

    template<>
    void CSMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp atomic
        Ax[Ap[n] + pos] += v;
      }
    }

    template<>
    void CSMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      if(v != 0.0)   // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if(pos < 0)
        {
          this->info("CSMatrix<Scalar>::add(): i = %d, j = %d.", m, n);
          throw Hermes::Exceptions::Exception("Sparse matrix entry not found: [%i, %i]", m, n);
        }

#pragma omp critical
        Ax[Ap[n] + pos] += v;
      }
    }

    template<typename Scalar>
    Scalar CSMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      // Find m-th row in the n-th column.
      int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);

      if(mid < 0) // if the entry has not been found
        return 0.0;
      else
        return Ax[Ap[n] + mid];
    }



    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix() : CSMatrix<Scalar>()
    {
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::CSCMatrix(unsigned int size) : CSMatrix<Scalar>(size)
    {
    }

    template<typename Scalar>
    CSCMatrix<Scalar>::~CSCMatrix()
    {
    }

    template<>
    void CSCMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      CSMatrix<double>::add(m, n, v);
    }

    template<>
    void CSCMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      CSMatrix<std::complex<double> >::add(m, n, v);
    }

    template<typename Scalar>
    Scalar CSCMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      return CSMatrix<Scalar>::get(m, n);
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
          const double zero_cutoff = Hermes::epsilon;
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
          const double zero_cutoff = Hermes::epsilon;
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
    void CSCMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, CSMatrix<Scalar>* mat)
    {
      CSCMatrix<Scalar>* mat_cast = dynamic_cast<CSCMatrix<Scalar>*>(mat);
      if(!mat_cast)
        throw Hermes::Exceptions::Exception("Wrong matrix type detected in CSCMatrix<Scalar>::add_as_block().");

      Hermes::Solvers::CSCIterator<Scalar> mat_it(mat_cast);
      Hermes::Solvers::CSCIterator<Scalar> this_it(this);

      // Sanity check.
      bool this_not_empty = this_it.init();
      if(!this_not_empty)
        throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block().");

      // Iterate through the small matrix column by column and add all nonzeros
      // to the large one.
      bool mat_not_finished = mat_it.init();
      if(!mat_not_finished)
        throw Hermes::Exceptions::Exception("Empty matrix detected in CSCMatrix<Scalar>::add_as_block().");

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
    void CSCMatrix<Scalar>::add_matrix(CSMatrix<Scalar>* mat)
    {
      assert(this->get_size() == mat->get_size());

      // Create iterators for both matrices.
      CSCMatrix<Scalar>* mat_cast = dynamic_cast<CSCMatrix<Scalar>*>(mat);
      if(!mat_cast)
        throw Hermes::Exceptions::Exception("Wrong matrix type detected in CSCMatrix<Scalar>::add_as_block().");

      Hermes::Solvers::CSCIterator<Scalar> mat_it(mat_cast);
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
    CSMatrix<Scalar>* CSCMatrix<Scalar>::duplicate()
    {
      CSCMatrix<Scalar>* new_matrix = new CSCMatrix<Scalar>();
      new_matrix->create(this->get_size(), this->get_nnz(), this->get_Ap(),  this->get_Ai(),  this->get_Ax());
      return new_matrix;
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::CSRMatrix() : CSMatrix<Scalar>()
    {
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::CSRMatrix(unsigned int size) : CSMatrix<Scalar>(size)
    {
    }

    template<typename Scalar>
    CSRMatrix<Scalar>::~CSRMatrix()
    {
    }

    template<>
    void CSRMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      CSMatrix<double>::add(n, m, v);
    }

    template<>
    void CSRMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      CSMatrix<std::complex<double> >::add(n, m, v);
    }

    template<typename Scalar>
    Scalar CSRMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      return CSMatrix<Scalar>::get(n, m);
    }

    template<>
    bool CSRMatrix<double>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp =[\n",
          this->size, this->size, nnz, nnz);
        for (unsigned int j = 0; j < this->size; j++)
          for (int i = Ap[j]; i < Ap[j + 1]; i++)
          {
            fprintf(file, "%d %d ", j + 1, Ai[i] + 1);
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
                      fprintf(file, "%d %d ", (int)j + 1, Ai[i] + 1);
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
          const double zero_cutoff = Hermes::epsilon;
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
                ascii_entry_i[k] = j;
                ascii_entry_j[k] = Ai[i];
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

    template<typename Scalar>
    void CSRMatrix<Scalar>::pre_add_ij(unsigned int row, unsigned int col)
    {
      if(this->pages[row] == NULL || this->pages[row]->count >= SparseMatrix<Scalar>::PAGE_SIZE)
      {
        typename SparseMatrix<Scalar>::Page *new_page = new typename SparseMatrix<Scalar>::Page;
        new_page->count = 0;
        new_page->next = this->pages[row];
        this->pages[row] = new_page;
      }
      this->pages[row]->idx[this->pages[row]->count++] = col;
    }

    template<>
    bool CSRMatrix<std::complex<double> >::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      throw Exceptions::MethodNotImplementedException("CSRMatrix<double>::dump");
      return false;
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::add_as_block(unsigned int offset_i, unsigned int offset_j, CSMatrix<Scalar>* mat)
    {
      throw Exceptions::MethodNotImplementedException("CSRMatrix<double>::add_as_block");
    }

    template<typename Scalar>
    void CSRMatrix<Scalar>::add_matrix(CSMatrix<Scalar>* mat)
    {
      throw Exceptions::MethodNotImplementedException("CSRMatrix<double>::add_matrix");
    }

    template<typename Scalar>
    CSMatrix<Scalar>* CSRMatrix<Scalar>::duplicate()
    {
      throw Exceptions::MethodNotImplementedException("CSRMatrix<double>::duplicate");
    }
  }

  namespace Solvers
  {
    template<typename Scalar>
    bool CSIterator<Scalar>::init()
    {
      if(this->size == 0 || this->nnz == 0) return false;
      this->Ap_pos = 0;
      this->Ai_pos = 0;
      return true;
    }

    template<typename Scalar>
    CSIterator<Scalar>::CSIterator(Hermes::Algebra::CSMatrix<Scalar>* mat)
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
    bool CSIterator<Scalar>::move_ptr()
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
    void CSIterator<Scalar>::add_to_current_position(Scalar val)
    {
      this->Ax[this->Ai_pos] += val;
    }

    template<typename Scalar>
    CSCIterator<Scalar>::CSCIterator(Hermes::Algebra::CSCMatrix<Scalar>* mat) : CSIterator<Scalar>(mat)
    {
    }

    template<typename Scalar>
    void CSCIterator<Scalar>::get_current_position(int& i, int& j, Scalar& val)
    {
      i = this->Ai[this->Ai_pos];
      j = this->Ap_pos;
      val = this->Ax[this->Ai_pos];
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
    CSRIterator<Scalar>::CSRIterator(Hermes::Algebra::CSRMatrix<Scalar>* mat) : CSIterator<Scalar>(mat)
    {
    }

    template<typename Scalar>
    void CSRIterator<Scalar>::get_current_position(int& i, int& j, Scalar& val)
    {
      i = this->Ai[this->Ai_pos];
      j = this->Ap_pos;
      val = this->Ax[this->Ai_pos];
    }

    template<typename Scalar>
    bool CSRIterator<Scalar>::move_to_position(int i, int j)
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
  }
}

template class HERMES_API Hermes::Algebra::CSMatrix<double>;
template class HERMES_API Hermes::Algebra::CSMatrix<std::complex<double> >;

template class HERMES_API Hermes::Algebra::CSCMatrix<double>;
template class HERMES_API Hermes::Algebra::CSCMatrix<std::complex<double> >;

template class HERMES_API Hermes::Algebra::CSRMatrix<double>;
template class HERMES_API Hermes::Algebra::CSRMatrix<std::complex<double> >;

template class HERMES_API Hermes::Solvers::CSIterator<double>;
template class HERMES_API Hermes::Solvers::CSIterator<std::complex<double> >;

template class HERMES_API Hermes::Solvers::CSCIterator<double>;
template class HERMES_API Hermes::Solvers::CSCIterator<std::complex<double> >;

template class HERMES_API Hermes::Solvers::CSRIterator<double>;
template class HERMES_API Hermes::Solvers::CSRIterator<std::complex<double> >;