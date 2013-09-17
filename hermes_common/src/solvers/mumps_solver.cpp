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

namespace Hermes
{
  namespace Algebra
  {
    extern "C"
    {
#ifndef _WINDOWS
      extern void dmumps_c(DMUMPS_STRUC_C *mumps_param_ptr);
      extern void zmumps_c(ZMUMPS_STRUC_C *mumps_param_ptr);
#endif
    }

#define USE_COMM_WORLD  -987654

    /// Binary search for the location of a particular CSC/CSR matrix entry.
    ///
    /// Typically, we search for the index into Ax that corresponds to a given
    /// row (CSC) or column (CSR) ('idx') among indices of nonzero values in
    /// a particular column (CSC) or row (CSR) ('Ai').
    ///
    static int find_position(int *Ai, int Alen, int idx)
    {
      assert(idx >= 0);

      register int lo = 0, hi = Alen - 1, mid;

      while (1)
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
    MumpsMatrix<Scalar>::MumpsMatrix()
    {
      nnz = 0;
      this->size = 0;
      irn = NULL;
      jcn = NULL;
      Ax = NULL;
      Ap = NULL;
      Ai = NULL;
    }

    template<typename Scalar>
    MumpsMatrix<Scalar>::~MumpsMatrix()
    {
      free();
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::alloc()
    {
      assert(this->pages != NULL);

      // initialize the arrays Ap and Ai
      Ap = new int[this->size + 1];
      int aisize = this->get_num_indices();
      Ai = new int[aisize];

      // sort the indices and remove duplicities, insert into Ai
      unsigned int i, pos = 0;
      for (i = 0; i < this->size; i++)
      {
        Ap[i] = pos;
        pos += this->sort_and_store_indices(this->pages[i], Ai + pos, Ai + aisize);
      }
      Ap[i] = pos;

      delete [] this->pages;
      this->pages = NULL;

      nnz = Ap[this->size];

      Ax = new typename mumps_type<Scalar>::mumps_Scalar[nnz];
      memset(Ax, 0, sizeof(Scalar) * nnz);

      irn = new int[nnz];
      jcn = new int[nnz];
      for (unsigned int i = 0; i < nnz; i++)
      {
        irn[i] = 1;
        jcn[i] = 1;
      }
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::free()
    {
      nnz = 0;
      delete [] Ap; Ap = NULL;
      delete [] Ai; Ai = NULL;
      delete [] Ax; Ax = NULL;
      delete [] irn; irn = NULL;
      delete [] jcn; jcn = NULL;
    }

    inline double mumps_to_Scalar(double x)
    {
      return x;
    }

    inline std::complex<double> mumps_to_Scalar(ZMUMPS_COMPLEX x)
    {
      return std::complex<double>(x.r, x.i);
    }

    template<typename Scalar>
    Scalar MumpsMatrix<Scalar>::get(unsigned int m, unsigned int n) const
    {
      // Find m-th row in the n-th column.
      int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
      // Return 0 if the entry has not been found.
      if(mid < 0) return 0.0;
      // Otherwise, add offset to the n-th column and return the value.
      if(mid >= 0) mid += Ap[n];
      return mumps_to_Scalar(Ax[mid]);
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::zero()
    {
      memset(Ax, 0, sizeof(Scalar) * Ap[this->size]);
    }

    inline ZMUMPS_COMPLEX& operator +=(ZMUMPS_COMPLEX &a, std::complex<double> b)
    {
      a.r +=b.real();
      a.i +=b.imag();
      return a;
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar v)
    {
      //          produced an error in neutronics-2-group-adapt (although tutorial-07
      //          ran well).
      // Find m-th row in the n-th column.
      int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
      // Make sure we are adding to an existing non-zero entry.
      if(pos < 0)
        throw Hermes::Exceptions::Exception("Sparse matrix entry not found");
      // Add offset to the n-th column.
      pos += Ap[n];
#pragma omp critical (MumpsMatrix_add)
      Ax[pos] += v;
      irn[pos] = m + 1;  // MUMPS is indexing from 1
      jcn[pos] = n + 1;
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      for (unsigned int i = 0; i < m; i++)       // rows
        for (unsigned int j = 0; j < n; j++)     // cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real\n");
          fprintf(file, "%d %d %d\n", this->size, this->size, this->nnz);

          for (unsigned int j = 0; j < this->size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              Helpers::fprint_coordinate_num(file, irn[i]-1, jcn[i]-1, mumps_to_Scalar(Ax[i]), number_format);
              fprintf(file, "\n");
            }
          }

          fclose(file);
        }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
        {
#ifdef WITH_MATIO
          mat_sparse_t sparse;
          sparse.nzmax = this->nnz;

          // For complex.
          double* Ax_re = NULL;
          double* Ax_im = NULL;

          sparse.nir = this->nnz;
          sparse.ir = Ai;
          sparse.njc = this->size + 1;
          sparse.jc = (int *) Ap;
          sparse.ndata = this->nnz;

          size_t dims[2];
          dims[0] = this->size;
          dims[1] = this->size;

          mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);

          matvar_t *matvar;

          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
          {
            sparse.data = Ax;
            matvar = Mat_VarCreate("matrix", MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA);
          }
          else
          {
            Ax_re = new double[this->nnz];
            Ax_im = new double[this->nnz];
            struct mat_complex_split_t z = {Ax_re, Ax_im};

            for(int i = 0; i < this->nnz; i++)
            {
              Ax_re[i] = ((std::complex<double>)(mumps_to_Scalar(this->Ax[i]))).real();
              Ax_im[i] = ((std::complex<double>)(mumps_to_Scalar(this->Ax[i]))).imag();
              sparse.data = &z;
            }
            matvar = Mat_VarCreate("matrix", MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
          }

          if (matvar)
          {
            Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
          }

          if(Ax_re)
            delete [] Ax_re;
          if(Ax_im)
            delete [] Ax_im;
          Mat_Close(mat);

          if(!matvar)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#endif
        }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          for (unsigned int j = 0; j < this->size; j++)
          {
            for (int i = Ap[j]; i < Ap[j + 1]; i++)
            {
              Helpers::fprint_coordinate_num(file, irn[i]-1, jcn[i]-1, mumps_to_Scalar(Ax[i]), number_format);
              fprintf(file, "\n");
            }
          }

          fclose(file);
        }
      }
    }

    template<typename Scalar>
    unsigned int MumpsMatrix<Scalar>::get_nnz() const
    {
      return nnz;
    }

    template<typename Scalar>
    double MumpsMatrix<Scalar>::get_fill_in() const
    {
      return Ap[this->size] / (double) (this->size * this->size);
    }

    inline ZMUMPS_COMPLEX& operator +=(ZMUMPS_COMPLEX &a, ZMUMPS_COMPLEX b)
    {
      a.r +=b.r;
      a.i +=b.i;
      return a;
    }

    template<typename Scalar>
    void MumpsMatrix<Scalar>::add_as_block(unsigned int i, unsigned int j, MumpsMatrix<Scalar>* mat)
    {
      int idx;
      for (unsigned int col = 0; col<mat->get_size(); col++)
      {
        for (unsigned int n = mat->Ap[col]; n < mat->Ap[col + 1]; n++)
        {
          idx = find_position(Ai + Ap[col + j], Ap[col + 1 + j] - Ap[col + j], mat->Ai[n] + i);
          if(idx < 0)
            throw Hermes::Exceptions::Exception("Sparse matrix entry not found");
          idx += Ap[col + j];
          Ax[idx] += mat->Ax[n];
        }
      }
    }

    // Applies the matrix to vector_in and saves result to vector_out.
    template<typename Scalar>
    void MumpsMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized) const
    {
      if(!vector_out_initialized)
        vector_out = new Scalar[this->size];
      for(unsigned int i = 0; i < this->size; i++)
        vector_out[i] = Scalar(0.);
      Scalar a;
      for (unsigned int i = 0;i<nnz;i++)
      {
        a = mumps_to_Scalar(Ax[i]);
        vector_out[jcn[i] - 1] += vector_in[irn[i] - 1] * a;
      }
    }

    template<>
    void MumpsMatrix<double>::multiply_with_Scalar(double value)
    {
      int n = nnz;
      for(int i = 0;i<n;i++)
      {
        Ax[i] = Ax[i] * value;
      }
    }

    template<>
    void MumpsMatrix<std::complex<double> >::multiply_with_Scalar(std::complex<double> value)
    {
      int n = nnz;
      std::complex<double> a;
      for(int i = 0;i<n;i++)
      {
        a = std::complex<double>(Ax[i].r, Ax[i].i);
        a = a*value;
        Ax[i].r = a.real();
        Ax[i].i = a.imag();
      }
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
    void MumpsMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
    {
      this->nnz = nnz;
      this->size = size;
      this->Ap = new int[this->size + 1]; assert(this->Ap != NULL);
      this->Ai = new int[nnz];    assert(this->Ai != NULL);
      this->Ax = new typename mumps_type<Scalar>::mumps_Scalar[nnz]; assert(this->Ax != NULL);
      irn = new int[nnz];           assert(this->irn !=NULL);     // Row indices.
      jcn = new int[nnz];           assert(this->jcn !=NULL);     // Column indices.

      for (unsigned int i = 0; i < this->size; i++)
      {
        this->Ap[i] = ap[i];
        for (int j = ap[i];j<ap[i + 1];j++) jcn[j] = i;
      }
      this->Ap[this->size] = ap[this->size];
      for (unsigned int i = 0; i < nnz; i++)
      {
        mumps_assign_Scalar(this->Ax[i], ax[i]);
        this->Ai[i] = ai[i];
        irn[i] = ai[i];
      }
    }
    // Duplicates a matrix (including allocation).
    template<typename Scalar>
    SparseMatrix<Scalar>* MumpsMatrix<Scalar>::duplicate() const
    {
      MumpsMatrix<Scalar> * nmat = new MumpsMatrix<Scalar>();

      nmat->nnz = nnz;
      nmat->size = this->size;
      nmat->Ap = new int[this->size + 1]; assert(nmat->Ap != NULL);
      nmat->Ai = new int[nnz];    assert(nmat->Ai != NULL);
      nmat->Ax = new typename mumps_type<Scalar>::mumps_Scalar[nnz]; assert(nmat->Ax != NULL);
      nmat->irn = new int[nnz];           assert(nmat->irn !=NULL);     // Row indices.
      nmat->jcn = new int[nnz];           assert(nmat->jcn !=NULL);     // Column indices.
      for (unsigned int i = 0;i<nnz;i++)
      {
        nmat->Ai[i] = Ai[i];
        nmat->Ax[i] = Ax[i];
        nmat->irn[i] = irn[i];
        nmat->jcn[i] = jcn[i];
      }
      for (unsigned int i = 0;i<this->size + 1;i++)
      {
        nmat->Ap[i] = Ap[i];
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
      case -1: this->warn("Error occured on processor %d", MUMPS_INFO(param, 2)); break;
        /// \todo add the rest according to the MUMPS docs
      default: this->warn("INFOG(1) = %d", param.INFOG(1)); break;
      }
      return false;
    }

    template<typename Scalar>
    bool MumpsSolver<Scalar>::reinit()
    {
      if(inited)
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

      if(inited)
      {
        // No printings.
        param.ICNTL(1) = -1;
        param.ICNTL(2) = -1;
        param.ICNTL(3) = -1;
        param.ICNTL(4) = 0;

        param.ICNTL(20) = 0; // centralized dense RHS
        param.ICNTL(21) = 0; // centralized dense solution

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
    MumpsSolver<Scalar>::MumpsSolver(MumpsMatrix<Scalar> *m, SimpleVector<Scalar> *rhs) :
      DirectSolver<Scalar>(), m(m), rhs(rhs)
    {
      inited = false;

      // Initial values for some fields of the MUMPS_STRUC structure that may be accessed
      // before MUMPS has been initialized.
      param.rhs = NULL;
      param.INFOG(33) = -999; // see the case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING
      // in setup_factorization()
    }

    template<typename Scalar>
    MumpsSolver<Scalar>::~MumpsSolver()
    {
      // Terminate the current instance of MUMPS.
      if(inited)
      {
        param.job = JOB_END;
        mumps_c(&param);
      }

      if(param.rhs != NULL)
        delete [] param.rhs;
    }

    template<typename Scalar>
    int MumpsSolver<Scalar>::get_matrix_size()
    {
      return m->size;
    }

    template<typename Scalar>
    void MumpsSolver<Scalar>::solve()
    {
      bool ret = false;
      assert(m != NULL);
      assert(rhs != NULL);

      this->tick();

      // Prepare the MUMPS data structure with input for the solver driver
      // (according to the chosen factorization reuse strategy), as well as
      // the system matrix.
      if( !setup_factorization() )
        throw Hermes::Exceptions::LinearMatrixSolverException("LU factorization could not be completed.");

      // Specify the right-hand side (will be replaced by the solution).
      param.rhs = new typename mumps_type<Scalar>::mumps_Scalar[m->size];
      memcpy(param.rhs, rhs->v, m->size * sizeof(Scalar));

      // Do the jobs specified in setup_factorization().
      mumps_c(&param);

      ret = check_status();

      if(ret)
      {
        delete [] this->sln;
        this->sln = new Scalar[m->size];
        for (unsigned int i = 0; i < rhs->get_size(); i++)
          this->sln[i] = mumps_to_Scalar(param.rhs[i]);
      }
      else
        throw Hermes::Exceptions::LinearMatrixSolverException("MUMPS failed.");

      this->tick();
      this->time = this->accumulated();

      delete [] param.rhs;
      param.rhs = NULL;
    }

    template<typename Scalar>
    bool MumpsSolver<Scalar>::setup_factorization()
    {
      // When called for the first time, all three phases (analysis, factorization,
      // solution) must be performed.
      int eff_fact_scheme = this->reuse_scheme;
      if(!inited)
        if( this->reuse_scheme == HERMES_REUSE_MATRIX_REORDERING ||
          this->reuse_scheme == HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY )
          eff_fact_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH;

      switch (eff_fact_scheme)
      {
      case HERMES_CREATE_STRUCTURE_FROM_SCRATCH:
        // (Re)initialize new instance.
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
        // and then reuse it during subsequent factorizations. New instance of MUMPS
        // has to be created before the analysis phase.
        if(param.INFOG(33) != -2)
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
