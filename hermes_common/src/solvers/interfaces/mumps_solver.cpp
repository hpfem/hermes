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
/*! \file mumps_solver.cpp
\brief MUMPS solver interface.
*/
#include "config.h"
#ifdef WITH_MUMPS
#include "mumps_solver.h"
#include "callstack.h"
#include "util/memory_handling.h"

namespace Hermes
{
  namespace Algebra
  {
    inline double mumps_to_Scalar(double x)
    {
      return x;
    }

    inline std::complex<double> mumps_to_Scalar(ZMUMPS_COMPLEX x)
    {
      return std::complex<double>(x.r, x.i);
    }

    double inline real(double x)
    {
      return x;
    }

    double inline imag(double x)
    {
      return 0;
    }

    double inline real(ZMUMPS_COMPLEX x)
    {
      return mumps_to_Scalar(x).real();
    }

    double inline imag(ZMUMPS_COMPLEX x)
    {
      return mumps_to_Scalar(x).imag();
    }

    extern "C"
    {
#ifndef _WINDOWS
      extern void dmumps_c(DMUMPS_STRUC_C *mumps_param_ptr);
      extern void zmumps_c(ZMUMPS_STRUC_C *mumps_param_ptr);
#endif
    }

#define USE_COMM_WORLD  -987654

    inline ZMUMPS_COMPLEX& operator +=(ZMUMPS_COMPLEX &a, double b)
    {
      a.i += b;
      return a;
    }

    inline ZMUMPS_COMPLEX& operator +=(ZMUMPS_COMPLEX &a, std::complex<double> b)
    {
      a.r += b.real();
      a.i += b.imag();
      return a;
    }

    inline ZMUMPS_COMPLEX& operator +=(ZMUMPS_COMPLEX &a, ZMUMPS_COMPLEX b)
    {
      a.r += b.r;
      a.i += b.i;
      return a;
    }

    inline ZMUMPS_COMPLEX& operator *=(ZMUMPS_COMPLEX &a, std::complex<double> b)
    {
      std::complex<double> a_c = std::complex<double>(a.r, a.i);
      std::complex<double> result = a_c * b;
      a.r = result.real();
      a.i = result.imag();
      return a;
    }

    inline void mumps_assign_Scalar(ZMUMPS_COMPLEX & a, std::complex<double> b)
    {
      a.r = b.real();
      a.i = b.imag();
    }

    inline void mumps_assign_Scalar(double & a, double b)
    {
      a = b;
    }

    template<typename Scalar>
    MumpsMatrix<Scalar>::MumpsMatrix() : CSCMatrix<Scalar>(), irn(nullptr), jcn(nullptr), Ax(nullptr)
    {
    }

    template<typename Scalar>
    MumpsMatrix<Scalar>::~MumpsMatrix()
    {
      this->free();
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::alloc_data()
    {
      this->Ax = calloc_with_check<MumpsMatrix<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(this->nnz, this);

      irn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
      jcn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
      for (unsigned int i = 0; i < this->nnz; i++)
      {
        irn[i] = 1;
        jcn[i] = 1;
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::free()
    {
      CSCMatrix<Scalar>::free();
      free_with_check(Ax);
      free_with_check(irn);
      free_with_check(jcn);
    }

    template<typename Scalar>
    Scalar MumpsMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      // Find m-th row in the n-th column.
      int mid = CSMatrix<Scalar>::find_position(this->Ai + this->Ap[n], this->Ap[n + 1] - this->Ap[n], m);
      // Return 0 if the entry has not been found.
      if (mid < 0)
        return 0.0;
      else
        mid += this->Ap[n];

      return mumps_to_Scalar(Ax[mid]);
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::zero()
    {
      memset(this->Ax, 0, sizeof(typename mumps_type<Scalar>::mumps_Scalar) * this->Ap[this->size]);
    }

    template<>
    void MumpsMatrix<double>::add(unsigned int m, unsigned int n, double v)
    {
      // Find m-th row in the n-th column.
      int pos = CSMatrix<double>::find_position(this->Ai + this->Ap[n], this->Ap[n + 1] - this->Ap[n], m);
      // Make sure we are adding to an existing non-zero entry.
      if (pos < 0)
        throw Hermes::Exceptions::Exception("Sparse matrix entry not found");
      // Add offset to the n-th column.
      pos += this->Ap[n];
#pragma omp atomic
      Ax[pos] += v;
      irn[pos] = m + 1;  // MUMPS is indexing from 1
      jcn[pos] = n + 1;
    }

    template<>
    void MumpsMatrix<std::complex<double> >::add(unsigned int m, unsigned int n, std::complex<double> v)
    {
      // Find m-th row in the n-th column.
      int pos = CSMatrix<std::complex<double> >::find_position(this->Ai + this->Ap[n], this->Ap[n + 1] - this->Ap[n], m);
      // Make sure we are adding to an existing non-zero entry.
      if (pos < 0)
        throw Hermes::Exceptions::Exception("Sparse matrix entry not found");
      // Add offset to the n-th column.
      pos += this->Ap[n];
#pragma omp critical (MumpsMatrix_add)
      Ax[pos] += v;
      irn[pos] = m + 1;  // MUMPS is indexing from 1
      jcn[pos] = n + 1;
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      bool invert_storage = false;
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
      {
                                        FILE* file = fopen(filename, "w");
                                        if (!file)
                                          throw Exceptions::IOException(Exceptions::IOException::Write, filename);
                                        if (Hermes::Helpers::TypeIsReal<Scalar>::value)
                                          fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
                                        else
                                          fprintf(file, "%%%%MatrixMarket matrix coordinate complex general\n");

                                        fprintf(file, "%d %d %d\n", this->size, this->size, this->nnz);

                                        if (invert_storage)
                                          this->switch_orientation();
                                        for (unsigned int j = 0; j < this->size; j++)
                                        {
                                          for (int i = this->Ap[j]; i < this->Ap[j + 1]; i++)
                                          {
                                            Hermes::Helpers::fprint_coordinate_num(file, this->Ai[i] + 1, j + 1, mumps_to_Scalar(this->Ax[i]), number_format);
                                            fprintf(file, "\n");
                                          }
                                        }
                                        if (invert_storage)
                                          this->switch_orientation();

                                        fclose(file);
      }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
      {
#ifdef WITH_MATIO
                                       mat_sparse_t sparse;
                                       sparse.nzmax = this->nnz;
                                       if (invert_storage)
                                         this->switch_orientation();

                                       sparse.nir = this->nnz;
                                       sparse.ir = this->Ai;
                                       sparse.njc = this->size;
                                       sparse.jc = (int *)this->Ap;
                                       sparse.ndata = this->nnz;

                                       size_t dims[2];
                                       dims[0] = this->size;
                                       dims[1] = this->size;

                                       mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);

                                       matvar_t *matvar;

                                       // For complex. No allocation here.
                                       double* Ax_re = nullptr;
                                       double* Ax_im = nullptr;

                                       // For real.
                                       if (Hermes::Helpers::TypeIsReal<Scalar>::value)
                                       {
                                         sparse.data = this->Ax;
                                         matvar = Mat_VarCreate(var_name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA);
                                       }
                                       else
                                       {
                                         // For complex.
                                         Ax_re = malloc_with_check<CSMatrix<Scalar>, double>(this->nnz, this);
                                         Ax_im = malloc_with_check<CSMatrix<Scalar>, double>(this->nnz, this);
                                         struct mat_complex_split_t z = { Ax_re, Ax_im };

                                         for (int i = 0; i < this->nnz; i++)
                                         {
                                           Ax_re[i] = ((std::complex<double>)(mumps_to_Scalar(this->Ax[i]))).real();
                                           Ax_im[i] = ((std::complex<double>)(mumps_to_Scalar(this->Ax[i]))).imag();
                                           sparse.data = &z;
                                         }
                                         matvar = Mat_VarCreate(var_name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
                                       }

                                       if (matvar)
                                       {
                                         Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
                                         Mat_VarFree(matvar);
                                       }
                                       if (invert_storage)
                                         this->switch_orientation();
                                        free_with_check(Ax_re);
                                        free_with_check(Ax_im);
                                       Mat_Close(mat);

                                       if (!matvar)
                                         throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#endif
      }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
      {
                                      FILE* file = fopen(filename, "w");
                                      if (!file)
                                        throw Exceptions::IOException(Exceptions::IOException::Write, filename);

                                      if (invert_storage)
                                        this->switch_orientation();
                                      for (unsigned int j = 0; j < this->size; j++)
                                      {
                                        for (int i = this->Ap[j]; i < this->Ap[j + 1]; i++)
                                        {
                                          Helpers::fprint_coordinate_num(file, this->Ai[i], j, mumps_to_Scalar(this->Ax[i]), number_format);
                                          fprintf(file, "\n");
                                        }
                                      }
                                      if (invert_storage)
                                        this->switch_orientation();

                                      fclose(file);
      }

#ifdef WITH_BSON
      case EXPORT_FORMAT_BSON:
      {
                               // Init bson
                               bson bw;
                               bson_init(&bw);

                               // Matrix size.
                               bson_append_int(&bw, "size", this->size);
                               // Nonzeros.
                               bson_append_int(&bw, "nnz", this->nnz);

                               if (invert_storage)
                                 this->switch_orientation();

                               bson_append_start_array(&bw, "Ap");
                               for (unsigned int i = 0; i < this->size; i++)
                                 bson_append_int(&bw, "p", this->Ap[i]);
                               bson_append_finish_array(&bw);

                               bson_append_start_array(&bw, "Ai");
                               for (unsigned int i = 0; i < this->nnz; i++)
                                 bson_append_int(&bw, "i", this->Ai[i]);
                               bson_append_finish_array(&bw);

                               bson_append_start_array(&bw, "Ax");
                               for (unsigned int i = 0; i < this->nnz; i++)
                                 bson_append_double(&bw, "x", real(this->Ax[i]));
                               bson_append_finish_array(&bw);

                               if (!Hermes::Helpers::TypeIsReal<Scalar>::value)
                               {
                                 bson_append_start_array(&bw, "Ax-imag");
                                 for (unsigned int i = 0; i < this->nnz; i++)
                                   bson_append_double(&bw, "x-i", imag(this->Ax[i]));
                                 bson_append_finish_array(&bw);
                               }
                               bson_append_finish_array(&bw);

                               if (invert_storage)
                                 this->switch_orientation();

                               // Done.
                               bson_finish(&bw);

                               // Write to disk.
                               FILE *fpw;
                               fpw = fopen(filename, "wb");
                               const char *dataw = (const char *)bson_data(&bw);
                               fwrite(dataw, bson_size(&bw), 1, fpw);
                               fclose(fpw);

                               bson_destroy(&bw);
      }
#endif
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      bool invert_storage = false;
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::import_from_file - Matrix Market");
        break;
      case EXPORT_FORMAT_MATLAB_MATIO:
      {
#ifdef WITH_MATIO
                                       mat_t    *matfp;
                                       matvar_t *matvar;

                                       matfp = Mat_Open(filename, MAT_ACC_RDONLY);

                                       if (!matfp)
                                         throw Exceptions::IOException(Exceptions::IOException::Read, filename);

                                       matvar = Mat_VarRead(matfp, var_name);

                                       if (matvar)
                                       {
                                         mat_sparse_t *sparse = (mat_sparse_t *)matvar->data;

                                         this->nnz = sparse->nir;
                                         this->size = sparse->njc;

                                         this->Ap = malloc_with_check<MumpsMatrix<Scalar>, int>(this->size + 1, this);
                                         this->Ai = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
                                         this->Ax = malloc_with_check<MumpsMatrix<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(this->nnz, this);
                                         this->irn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
                                         this->jcn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);

                                         void* data = nullptr;
                                         if (Hermes::Helpers::TypeIsReal<Scalar>::value)
                                           data = sparse->data;
                                         else
                                         {
                                           std::complex<double>* complex_data = malloc_with_check<CSMatrix<Scalar>, std::complex<double> >(this->nnz, this);
                                           double* real_array = (double*)((mat_complex_split_t*)sparse->data)->Re;
                                           double* imag_array = (double*)((mat_complex_split_t*)sparse->data)->Im;
                                           for (int i = 0; i < this->nnz; i++)
                                             complex_data[i] = std::complex<double>(real_array[i], imag_array[i]);
                                           data = (void*)complex_data;
                                         }
                                         memcpy(this->Ax, data, this->nnz * sizeof(Scalar));
                                         if (!Hermes::Helpers::TypeIsReal<Scalar>::value)
                                           free_with_check(data);
                                         memcpy(this->Ap, sparse->jc, this->size * sizeof(int));
                                         this->Ap[this->size] = this->nnz;
                                         memcpy(this->Ai, sparse->ir, this->nnz * sizeof(int));
                                         memcpy(this->irn, this->Ai, this->nnz * sizeof(int));

                                         for (unsigned int i = 0; i < this->size; i++)
                                         {
                                           for (int j = this->Ap[i]; j < this->Ap[i + 1]; j++)
                                             jcn[j] = i;
                                         }

                                         if (invert_storage)
                                           this->switch_orientation();
                                       }

                                       Mat_Close(matfp);

                                       if (!matvar)
                                         throw Exceptions::IOException(Exceptions::IOException::Read, filename);
#else
                                       throw Exceptions::Exception("MATIO not included.");
#endif
      }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        throw Exceptions::MethodNotOverridenException("CSMatrix<Scalar>::import_from_file - Simple format");
        break;
#ifdef WITH_BSON
      case EXPORT_FORMAT_BSON:
      {
                               FILE *fpr;
                               fpr = fopen(filename, "rb");

                               // file size:
                               fseek(fpr, 0, SEEK_END);
                               int size = ftell(fpr);
                               rewind(fpr);

                               // allocate memory to contain the whole file:
                               char *datar = malloc_with_check<char>(size);
                               fread(datar, size, 1, fpr);
                               fclose(fpr);

                               bson br;
                               bson_init_finished_data(&br, datar, 0);

                               bson_iterator it;
                               bson sub;
                               bson_find(&it, &br, "size");
                               this->size = bson_iterator_int(&it);
                               bson_find(&it, &br, "nnz");
                               this->nnz = bson_iterator_int(&it);

                               this->Ap = malloc_with_check<MumpsMatrix<Scalar>, int>(this->size + 1, this);
                               this->Ai = malloc_with_check<MumpsMatrix<Scalar>, int>(nnz, this);
                               this->Ax = malloc_with_check<MumpsMatrix<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(this->nnz, this);
                               this->irn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
                               this->jcn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);

                               // coeffs
                               bson_iterator it_coeffs;
                               bson_find(&it_coeffs, &br, "Ap");
                               bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                               bson_iterator_init(&it, &sub);
                               int index_coeff = 0;
                               while (bson_iterator_next(&it))
                                 this->Ap[index_coeff++] = bson_iterator_int(&it);
                               this->Ap[this->size] = this->nnz;

                               bson_find(&it_coeffs, &br, "Ai");
                               bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                               bson_iterator_init(&it, &sub);
                               index_coeff = 0;
                               while (bson_iterator_next(&it))
                                 this->Ai[index_coeff++] = bson_iterator_int(&it);

                               bson_find(&it_coeffs, &br, "Ax");
                               bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                               bson_iterator_init(&it, &sub);
                               index_coeff = 0;
                               while (bson_iterator_next(&it))
                                 mumps_assign_Scalar(this->Ax[index_coeff++], bson_iterator_double(&it));

                               if (!Hermes::Helpers::TypeIsReal<Scalar>::value)
                               {
                                 bson_find(&it_coeffs, &br, "Ax-imag");
                                 bson_iterator_subobject_init(&it_coeffs, &sub, 0);
                                 bson_iterator_init(&it, &sub);
                                 index_coeff = 0;
                                 while (bson_iterator_next(&it))
                                   this->Ax[index_coeff++] += bson_iterator_double(&it);
                               }

                               memcpy(this->irn, this->Ai, this->nnz * sizeof(int));

                               for (unsigned int i = 0; i < this->size; i++)
                               {
                                 for (int j = this->Ap[i]; j < this->Ap[i + 1]; j++)
                                   jcn[j] = i;
                               }

                               if (invert_storage)
                                 this->switch_orientation();

                               bson_destroy(&br);
                               free_with_check(datar);
      }
#endif
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const
    {
      if (!vector_out_initialized)
        vector_out = malloc_with_check<Scalar>(this->size);

      memset(vector_out, 0, sizeof(Scalar)* this->size);

      Scalar a;
      for (unsigned int i = 0; i < this->nnz; i++)
      {
        a = mumps_to_Scalar(Ax[i]);
        vector_out[jcn[i] - 1] += vector_in[irn[i] - 1] * a;
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      for (int i = 0; i < this->nnz; i++)
      {
        Ax[i] *= value;
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::create(unsigned int size, unsigned int nnz_, int* ap, int* ai, Scalar* ax)
    {
      this->nnz = nnz_;
      this->size = size;
      this->Ap = malloc_with_check<MumpsMatrix<Scalar>, int>(this->size + 1, this);
      this->Ai = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
      this->Ax = malloc_with_check<MumpsMatrix<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(this->nnz, this);
      irn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);
      jcn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, this);

      for (unsigned int i = 0; i < this->size; i++)
      {
        this->Ap[i] = ap[i];
        for (int j = ap[i]; j < ap[i + 1]; j++)
          jcn[j] = i;
      }
      this->Ap[this->size] = ap[this->size];
      for (unsigned int i = 0; i < this->nnz; i++)
      {
        mumps_assign_Scalar(this->Ax[i], ax[i]);
        this->Ai[i] = ai[i];
        irn[i] = ai[i];
      }
    }

    template<typename Scalar>
    CSMatrix<Scalar>* MumpsMatrix<Scalar>::duplicate() const
    {
      MumpsMatrix<Scalar> * nmat = new MumpsMatrix<Scalar>();

      nmat->nnz = this->nnz;
      nmat->size = this->size;
      nmat->Ap = malloc_with_check<MumpsMatrix<Scalar>, int>(this->size + 1, nmat);
      nmat->Ai = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, nmat);
      nmat->Ax = malloc_with_check<MumpsMatrix<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(this->nnz, nmat);
      nmat->irn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, nmat);
      nmat->jcn = malloc_with_check<MumpsMatrix<Scalar>, int>(this->nnz, nmat);
      for (unsigned int i = 0; i < this->nnz; i++)
      {
        nmat->Ai[i] = this->Ai[i];
        nmat->Ax[i] = Ax[i];
        nmat->irn[i] = irn[i];
        nmat->jcn[i] = jcn[i];
      }
      for (unsigned int i = 0; i < this->size + 1; i++)
      {
        nmat->Ap[i] = this->Ap[i];
      }
      return nmat;
    }

    template class HERMES_API MumpsMatrix<double>;
    template class HERMES_API MumpsMatrix<std::complex<double> >;
  }

  namespace Solvers
  {
    /// Macros allowing to use indices according to the Fortran documentation to index C arrays.
#define ICNTL(I)            icntl[(I)-1]
#define MUMPS_INFO(param, I) (param).infog[(I)-1]
#define INFOG(I)            infog[(I)-1]

    /// Job definitions according to MUMPS documentation.
#define JOB_INIT                    -1
#define JOB_END                     -2
#define JOB_ANALYZE_FACTORIZE_SOLVE  6
#define JOB_FACTORIZE_SOLVE          5
#define JOB_SOLVE                    3

    template<typename Scalar>
    MumpsSolver<Scalar>::MumpsSolver(MumpsMatrix<Scalar> *m, SimpleVector<Scalar> *rhs) :
      DirectSolver<Scalar>(m, rhs), m(m), rhs(rhs), icntl_14(init_icntl_14)
    {
        inited = false;

        // Initial values for some fields of the MUMPS_STRUC structure that may be accessed
        // before MUMPS has been initialized.
        param.rhs = nullptr;
        param.INFOG(33) = -999; // see the case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING
        // in setup_factorization()
      }

    template<typename Scalar>
    MumpsSolver<Scalar>::~MumpsSolver()
    {
      free();
    }

    template<typename Scalar>
    void MumpsSolver<Scalar>::free()
    {
      // Terminate the current instance of MUMPS.
      if (inited)
      {
        param.job = JOB_END;
        mumps_c(&param);
      }

      if (param.rhs != nullptr)
        free_with_check(param.rhs);
    }

    template<>
    void MumpsSolver<double>::mumps_c(mumps_type<double>::mumps_struct * param)
    {
      dmumps_c(param);
    }

    template<>
    void MumpsSolver<std::complex<double> >::mumps_c(mumps_type<std::complex<double> >::mumps_struct * param)
    {
      zmumps_c(param);
    }

    template<typename Scalar>
    bool MumpsSolver<Scalar>::check_status()
    {
      switch (param.INFOG(1))
      {
      case 0: return true; // no error
      case -1: throw Hermes::Exceptions::LinearMatrixSolverException("Error occured on processor %d", MUMPS_INFO(param, 2)); break;
      case -2: throw Hermes::Exceptions::LinearMatrixSolverException("Number of nonzeros (NNZ) is out of range."); break;
      case -3: throw Hermes::Exceptions::LinearMatrixSolverException("MUMPS called with an invalid option for JOB."); break;
      case -5: throw Hermes::Exceptions::LinearMatrixSolverException("Problem of REAL or COMPLEX workspace allocation of size %i during analysis.", param.INFOG(2)); break;
      case -6: throw Hermes::Exceptions::LinearMatrixSolverException("Matrix is singular in structure."); break;
      case -7: throw Hermes::Exceptions::LinearMatrixSolverException("Problem of INTEGER workspace allocation of size %i during analysis.", param.INFOG(2)); break;
      case -8:
      case -9: return false;
      case -10: throw Hermes::Exceptions::LinearMatrixSolverException("Numerically singular matrix."); break;
      default: Hermes::Exceptions::LinearMatrixSolverException("Non-detailed exception in MUMPS: INFOG(1) = %d", param.INFOG(1)); break;
      }
      return false;
    }

    template<typename Scalar>
    bool MumpsSolver<Scalar>::reinit()
    {
      if (inited)
      {
        // If there is already an instance of MUMPS running,
        // terminate it.
        param.job = JOB_END;
        mumps_c(&param);
      }

      param.job = JOB_INIT;
      param.par = 1; // host also performs calculations
      param.sym = 0; // 0 = unsymmetric
      param.comm_fortran = USE_COMM_WORLD;

      mumps_c(&param);
      inited = check_status();

      if (inited)
      {
        // No printings.
        param.ICNTL(1) = -1;
        param.ICNTL(2) = -1;
        param.ICNTL(3) = -1;
        param.ICNTL(4) = 0;

        param.ICNTL(5) = 0;  // =/ both centralized assembled matrix
        param.ICNTL(18) = 0; // =\ both centralized assembled matrix
        param.ICNTL(20) = 0; // centralized dense RHS
        param.ICNTL(21) = 0; // centralized dense solution

        // Fixing the memory problems - this parameter specifies the maximum
        // extra fill-in.
        // Extract from the docs (4.10, page 27):
        // ICNTL(14) is accessed by the host both during the analysis and the factorization phases. It corresponds
        // to the percentage increase in the estimated working space. When significant extra fill-in is caused
        // by numerical pivoting, increasing ICNTL(14) may help. Except in special cases, the default value
        // is 20 (which corresponds to a 20 % increase).
        param.ICNTL(14) = 100 * this->icntl_14;

        // Specify the matrix.
        param.n = m->size;
        param.nz = m->nnz;
        param.irn = m->irn;
        param.jcn = m->jcn;
        param.a = m->Ax;
      }

      return inited;
    }

    template<typename Scalar>
    int MumpsSolver<Scalar>::get_matrix_size()
    {
      return m->size;
    }

    template<typename Scalar>
    void MumpsSolver<Scalar>::solve()
    {
      assert(m != nullptr);
      assert(rhs != nullptr);

      this->tick();

      // Prepare the MUMPS data structure with input for the solver driver
      // (according to the chosen factorization reuse strategy), as well as
      // the system matrix.
      if (!setup_factorization())
        throw Hermes::Exceptions::LinearMatrixSolverException("LU factorization could not be completed.");

      // Specify the right-hand side (will be replaced by the solution).
      param.rhs = malloc_with_check<MumpsSolver<Scalar>, typename mumps_type<Scalar>::mumps_Scalar>(m->size, this);
      memcpy(param.rhs, rhs->v, m->size * sizeof(typename mumps_type<Scalar>::mumps_Scalar));

      // Do the jobs specified in setup_factorization().
      mumps_c(&param);

      // Throws appropriate exception.
      if (check_status())
      {
        free_with_check(this->sln);
        this->sln = malloc_with_check<MumpsSolver<Scalar>, Scalar>(m->size, this);
        for (unsigned int i = 0; i < rhs->get_size(); i++)
          this->sln[i] = mumps_to_Scalar(param.rhs[i]);
      }
      else
      {
        free_with_check(param.rhs);

        icntl_14 *= 2;
        if (icntl_14 > max_icntl_14)
          throw Hermes::Exceptions::LinearMatrixSolverException("MUMPS memory overflow - potentially singular matrix");
        else
        {
          this->reinit();
          this->solve();
        }
        return;
        /* From the MUMPS docs - these two cases are forwarded from check_status().
        case –8: throw Hermes::Exceptions::LinearMatrixSolverException("Main internal integer workarray IS too small for factorization. This may happen, for example, if
        numerical pivoting leads to significantly more fill-in than was predicted by the analysis. The user
        should increase the value of ICNTL(14) before calling the factorization again (JOB=2).
        case –9: Main internal real/complex workarray S too small. If INFO(2) is positive, then the number of entries
        that are missing in S at the moment when the error is raised is available in INFO(2). If INFO(2) is
        negative, then its absolute value should be multiplied by 1 million. If an error –9 occurs, the user
        should increase the value of ICNTL(14) before calling the factorization (JOB=2) again, except if
        ICNTL(23) is provided, in which case ICNTL(23) should be increased.
        */
      }

      this->tick();
      this->time = this->accumulated();

      free_with_check(param.rhs);
      param.rhs = nullptr;
    }

    template<typename Scalar>
    bool MumpsSolver<Scalar>::setup_factorization()
    {
      // When called for the first time, all three phases (analysis, factorization,
      // solution) must be performed.
      int eff_fact_scheme = this->reuse_scheme;
      if (!inited)
      if (this->reuse_scheme == HERMES_REUSE_MATRIX_REORDERING || this->reuse_scheme == HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY)
        eff_fact_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH;

      switch (eff_fact_scheme)
      {
      case HERMES_CREATE_STRUCTURE_FROM_SCRATCH:
        // (Re)initialize new_ instance.
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
        // and then reuse it during subsequent factorizations. new_ instance of MUMPS
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
      case HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY:
        param.job = JOB_SOLVE;
        break;
      }

      return true;
    }

    template class HERMES_API MumpsSolver<double>;
    template class HERMES_API MumpsSolver<std::complex<double> >;
  }
}
#endif
