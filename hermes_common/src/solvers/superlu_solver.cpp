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
/*! \file superlu_solver.cpp
\brief SuperLU solver interface.
*/
#include "config.h"
#ifdef WITH_SUPERLU
#include "superlu_solver.h"
#include "trace.h"
#include "error.h"
#include "callstack.h"

using namespace Hermes::Error;

namespace Hermes
{
  namespace Solvers
  {
#ifdef SLU_MT
    template <>
    void SuperLu<double>::sequ(SuperMatrix *A, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info)
    {
      dsequ (A, r, c, rowcnd, colcnd, amax, info);
    }

    template <>
    void SuperLu<double>::laqgs (SuperMatrix *A, float *r, float *c, float rowcnd, float colcnd, float amax, char *equed)
    {
      dlaqgs (A, r, c, rowcnd, colcnd, amax, equed);
    }

    template <>
    int SuperLu<double>::gstrf (superlu_options_t *options, int m, int n, double anorm, LUstruct_t *LUstruct, gridinfo_t *grid, SuperLUStat_t *stat, int *info)
    {
      return dgstrf (options, m, n, anorm, LUstruct, grid, stat, info);
    }

    template <>
    float SuperLu<double>::pivotGrowth (int ncols, SuperMatrix *A, int *perm_c, SuperMatrix *L, SuperMatrix *U)
    {
      return dPivotGrowth (ncols, A, perm_c, L, U);
    }

    template <>
    float SuperLu<double>::langs (char *norm, SuperMatrix *A)
    {
      return dlangs (norm, A);
    }
    template <>
    void  SuperLu<double>::gscon (char *norm, SuperMatrix *L, SuperMatrix *U, float anorm, float *rcond, SuperLUStat_t *stat, int *info)
    {
      dgscon (norm, L, U, anorm, rcond, stat, info);
    }

    template <>
    void  SuperLu<double>::gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r, SuperMatrix *B, SuperLUStat_t *stat, int *info)
    {
      dgstrs (trans, L, U, perm_c, perm_r, B, stat, info);
    }

    template <>
    double SuperLu<double>::lamch_ (char *cmach)
    {
      return dlamch_ (cmach);
    }

    template <>
    int SuperLu<double>::querySpace (SuperMatrix *a, SuperMatrix *b, mem_usage_t *mu)
    {
      return dquerySpace (a, b, mu);
    }
#endif
  }
  namespace Algebra
  {
    /// Binary search for the location of a particular CSC/CSR matrix entry.
    ///
    /// Typically, we search for the index into Ax that corresponds to a given
    /// row (CSC) or column (CSR) ('idx') among indices of nonzero values in
    /// a particular column (CSC) or row (CSR) ('Ai').
    ///
    static int find_position(int *Ai, int Alen, int idx)
    {
      _F_;
      assert (idx >= 0);

      register int lo = 0, hi = Alen - 1, mid;

      while (1)
      {
        mid = (lo + hi) >> 1;

        if (idx < Ai[mid]) hi = mid - 1;
        else if (idx > Ai[mid]) lo = mid + 1;
        else break;

        // Sparse matrix entry not found (raise an error when trying to add
        // value to this position, return 0 when obtaining value there).
        if (lo > hi) mid = -1;
      }
      return mid;
    }

    template<typename Scalar>
    SuperLUMatrix<Scalar>::SuperLUMatrix()
    {
      _F_;
      this->size = 0; nnz = 0;
      Ax = NULL;
      Ap = NULL;
      Ai = NULL;
    }

    template<typename Scalar>
    SuperLUMatrix<Scalar>::~SuperLUMatrix()
    {
      _F_;
      this->free();
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::alloc()
    {
      _F_;
      assert(this->pages != NULL);

      // Initialize the arrays Ap and Ai.
      Ap = new unsigned int [this->size + 1];
      MEM_CHECK(Ap);
      int aisize = this->get_num_indices();
      Ai = new int [aisize];
      MEM_CHECK(Ai);

      // sort the indices and remove duplicities, insert into Ai
      unsigned int i, pos = 0;
      for (i = 0; i < this->size; i++)
      {
        Ap[i] = pos;
        pos += sort_and_store_indices(this->pages[i], Ai + pos, Ai + aisize);
      }
      Ap[i] = pos;

      delete [] this->pages;
      this->pages = NULL;

      nnz = Ap[this->size];

      Ax = new Scalar [nnz];
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::free()
    {
      _F_;
      nnz = 0;
      delete [] Ap; Ap = NULL;
      delete [] Ai; Ai = NULL;
      delete [] Ax; Ax = NULL;
    }

    template<typename Scalar>
    Scalar SuperLUMatrix<Scalar>::get(unsigned int m, unsigned int n)
    {
      _F_;
      // Find m-th row in the n-th column.
      int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
      // Return 0 if the entry has not been found.
      if (mid < 0) return 0.0;
      // Otherwise, add offset to the n-th column and return the value.
      if (mid >= 0) mid += Ap[n];
      return Ax[mid];
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::zero()
    {
      _F_;
      memset(Ax, 0, sizeof(Scalar) * nnz);
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar v)
    {
      _F_;
      if (v != 0.0) // ignore zero values.
      {
        // Find m-th row in the n-th column.
        int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
        // Make sure we are adding to an existing non-zero entry.
        if (pos < 0)
          error("Sparse matrix entry not found");
        // Add offset to the n-th column.
        pos += Ap[n];
        Ax[pos] += v;
      }
    }

    /// Add a number to each diagonal entry.

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i<this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      _F_;
      for (unsigned int i = 0; i < m; i++)       // rows
        for (unsigned int j = 0; j < n; j++)     // cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    /// Save matrix and right-hand side to a file.
    ///
    template<typename Scalar>
    bool SuperLUMatrix<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      // TODO
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", this->size, this->size, Ap[this->size], Ap[this->size]);
        for (unsigned int j = 0; j < this->size; j++)
          for (unsigned int i = Ap[j]; i < Ap[j + 1]; i++)
          {
            fprintf(file, "%d %d " , Ai[i] + 1, j + 1);
            Hermes::Helpers::fprint_num(file, Ax[i]);
            fprintf(file, "\n");
          }
          fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

          return true;

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESX\001", 1, 8, file);
          int ssize = sizeof(Scalar);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(&nnz, sizeof(int), 1, file);
          hermes_fwrite(Ap, sizeof(int), this->size + 1, file);
          hermes_fwrite(Ai, sizeof(int), nnz, file);
          hermes_fwrite(Ax, sizeof(Scalar), nnz, file);
          return true;
        }

      default:
        return false;
      }
    }

    template<typename Scalar>
    unsigned int SuperLUMatrix<Scalar>::get_matrix_size() const
    {
      return this->size;
    }

    template<typename Scalar>
    unsigned int SuperLUMatrix<Scalar>::get_nnz() const
    {
      return nnz;
    }

    template<typename Scalar>
    double SuperLUMatrix<Scalar>::get_fill_in() const
    {
      _F_;
      return nnz / (double) (this->size * this->size);
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add_matrix(SuperLUMatrix<Scalar>* mat)
    {
      _F_;
      add_as_block(0, 0, mat);
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add_to_diagonal_blocks(int num_stages, SuperLUMatrix<Scalar>* mat)
    {
      _F_;
      int ndof = mat->get_size();
      if (this->get_size() != (unsigned int) num_stages * ndof)
        error("Incompatible matrix sizes in PetscMatrix<Scalar>::add_to_diagonal_blocks()");

      for (int i = 0; i < num_stages; i++)
      {
        this->add_as_block(ndof*i, ndof*i, mat);
      }
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      add_to_diagonal_blocks(num_stages, dynamic_cast<SuperLUMatrix*>(mat));
    }

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::add_as_block(unsigned int i, unsigned int j, SuperLUMatrix<Scalar>* mat)
    {
      _F_;
      int idx;
      for (unsigned int col = 0;col<mat->get_size();col++)
      {
        for (unsigned int n = mat->Ap[col];n<mat->Ap[col + 1];n++)
        {
          idx = find_position(Ai + Ap[col + j], Ap[col + 1 + j] - Ap[col + j], mat->Ai[n] + i);
          if (idx<0)
            error("Sparse matrix entry not found");
          idx += Ap[col + j];
          Ax[idx] +=mat->Ax[n];
        }
      }
    }

    // Applies the matrix to vector_in and saves result to vector_out.

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out)
    {
      _F_;
      for(unsigned int i = 0;i<this->size;i++)
      {
        vector_out[i] = 0;
      }
      for (unsigned int c = 0;c<this->size;c++)
      {
        for (unsigned int i = Ap[c];i<Ap[c + 1];i++)
        {
          vector_out[c] +=vector_in[Ai[i]]*Ax[i];
        }
      }
    }

    // Multiplies matrix with a Scalar.
    template<typename Scalar>
    void SuperLUMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      _F_;
      int n = nnz;
      for(int i = 0;i<n;i++)
      {
        Ax[i] = Ax[i]*value;
      }
    }
    // Creates matrix using size, nnz, and the three arrays.

    template<typename Scalar>
    void SuperLUMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
    {
      _F_;
      this->nnz = nnz;
      this->size = size;
      this->Ap = new unsigned int[this->size + 1]; assert(this->Ap != NULL);
      this->Ai = new int[nnz];    assert(this->Ai != NULL);
      this->Ax = new Scalar[nnz]; assert(this->Ax != NULL);

      for (unsigned int i = 0; i < this->size + 1; i++)
      {
        this->Ap[i] = ap[i];
      }
      for (unsigned int i = 0; i < nnz; i++)
      {
        this->Ax[i] = ax[i];
        this->Ai[i] = ai[i];
      }
    }
    // Duplicates a matrix (including allocation).

    template<typename Scalar>
    SuperLUMatrix<Scalar>* SuperLUMatrix<Scalar>::duplicate()
    {
      _F_;
      SuperLUMatrix<Scalar> * nmat = new SuperLUMatrix<Scalar>();

      nmat->nnz = nnz;
      nmat->size = this->size;
      nmat->Ap = new unsigned int[this->size + 1]; assert(nmat->Ap != NULL);
      nmat->Ai = new int[nnz];    assert(nmat->Ai != NULL);
      nmat->Ax = new Scalar[nnz]; assert(nmat->Ax != NULL);
      for (unsigned int i = 0;i<nnz;i++)
      {
        nmat->Ai[i] = Ai[i];
        nmat->Ax[i] = Ax[i];
      }
      for (unsigned int i = 0;i<this->size + 1;i++)
      {
        nmat->Ap[i] = Ap[i];
      }
      return nmat;
    }

    // SuperLUVector<Scalar> /////////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    SuperLUVector<Scalar>::SuperLUVector()
    {
      _F_;
      v = NULL;
      this->size = 0;
    }

    template<typename Scalar>
    SuperLUVector<Scalar>::~SuperLUVector()
    {
      _F_;
      this->free();
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::alloc(unsigned int n)
    {
      _F_;
      this->free();
      this->size = n;
      v = new Scalar[n];
      zero();
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::zero()
    {
      _F_;
      memset(v, 0, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::change_sign()
    {
      _F_;
      for (unsigned int i = 0; i < this->size; i++) v[i] *= -1.;
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::free()
    {
      _F_;
      delete [] v;
      v = NULL;
      this->size = 0;
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      _F_;
      v[idx] = y;
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::add(unsigned int idx, Scalar y)
    {
      _F_;
      v[idx] += y;
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      _F_;
      for (unsigned int i = 0; i < n; i++)
      {
        v[idx[i]] += y[i];
      }
    }

    template<typename Scalar>
    Scalar SuperLUVector<Scalar>::get(unsigned int idx)
    { 
      return v[idx]; 
    }
    
    template<typename Scalar>
    void SuperLUVector<Scalar>::extract(Scalar *v) const 
    { 
      memcpy(v, this->v, this->size * sizeof(Scalar)); 
    }
    
    template<typename Scalar>
    void SuperLUVector<Scalar>::add_vector(Vector<Scalar>* vec) 
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec->get(i));
    }

    template<typename Scalar>
    void SuperLUVector<Scalar>::add_vector(Scalar* vec) 
    {
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec[i]);
    }

    template<typename Scalar>
    bool SuperLUVector<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      switch (fmt)
      {
      case DF_NATIVE:
      case DF_PLAIN_ASCII:
        for (unsigned int i = 0; i < this->size; i++)
        {
          Hermes::Helpers::fprint_num(file, v[i]);
          fprintf(file, "\n");
        }

        return true;

      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx1\n%s = [\n", this->size, var_name);
        for (unsigned int i = 0; i < this->size; i++)
        {
          Hermes::Helpers::fprint_num(file, v[i]);
          fprintf(file, "\n");
        }
        fprintf(file, " ];\n");
        return true;

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESR\001", 1, 8, file);
          int ssize = sizeof(Scalar);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(v, sizeof(Scalar), this->size, file);
          return true;
        }

      default:
        return false;
      }
    }

    template class HERMES_API SuperLUMatrix<double>;
    template class HERMES_API SuperLUMatrix<std::complex<double> >;
    template class HERMES_API SuperLUVector<double>;
    template class HERMES_API SuperLUVector<std::complex<double> >;
  }
  namespace Solvers
  {
    template<typename Scalar>
    bool SuperLUSolver<Scalar>::check_status(unsigned int info)
    {
      _F_;
      if (info == 0)
      {
        // Success.
        return true;
      }
      else if (info <= m->size)
      {
        warning("SuperLU: Factor U is singular, solution could not be computed.");
        return false;
      }
      else if (info == m->size + 1)
      {
        warning("SuperLU: RCOND is less than machine precision "
          "(system matrix is singular to working precision).");
        return true;
      }
      else if (info > m->size + 1)
      {
        warning("SuperLU: Not enough memory.\n Failure when %.3f MB were allocated.",
          (info - m->size)/1e6);
        return false;
      }

      return false;
    }


    template<typename Scalar>
    SuperLUSolver<Scalar>::SuperLUSolver(SuperLUMatrix<Scalar> *m, SuperLUVector<Scalar> *rhs)
      : DirectSolver<Scalar>(HERMES_FACTORIZE_FROM_SCRATCH), m(m), rhs(rhs), local_Ai(NULL), local_Ap(NULL)
      , local_Ax(NULL), local_rhs(NULL)
    {
      _F_;
      R = NULL;
      C = NULL;
      perm_r = NULL;
      perm_c = NULL;
      etree = NULL;
#ifndef SLU_MT
      *equed = '\0';
#endif

      // Set the default input options:
#ifdef SLU_MT
      // I am not sure if this will work well on Windows:
      // http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c
      char *nt_var = getenv("OMP_NUM_THREADS");
      if (nt_var)
        options.nprocs          = std::max(1, atoi(nt_var));
      else
        options.nprocs          = 1;

      options.fact              = EQUILIBRATE;  // Rescale the matrix if neccessary.
      options.trans             = NOTRANS;      // Not solving the transposed problem.
      options.refact            = NO;           // Factorize from scratch for the first time.
      options.diag_pivot_thresh = 1.0;          // Use partial pivoting during GEM.
      options.usepr             = NO;           // Let SuperLU compute the row permutations.
      options.drop_tol          = 0.0;          // Not yet implemented in SuperLU_MT 2.0.
      options.SymmetricMode     = NO;           // Assume general non-symmetric problem.

      // Default options related to the supernodal algorithm.
      options.panel_size        = sp_ienv(1);
      options.relax             = sp_ienv(2);
#else
      /*
      options.Fact = DOFACT;
      options.Equil = YES;
      options.ColPerm = COLAMD;
      options.DiagPivotThresh = 1.0;
      options.Trans = NOTRANS;
      options.IterRefine = NOREFINE;
      options.SymmetricMode = NO;
      options.PivotGrowth = NO;
      options.ConditionNumber = NO;
      options.PrintStat = YES;
      */
      set_default_options(&options);  // This function is only present in the sequential SLU.
#endif

      options.PrintStat = YES;   // Set to NO to suppress output.

      has_A = has_B = inited = false;
    }

    inline SuperLuType<std::complex<double> >::Scalar to_superlu(SuperLuType<std::complex<double> >::Scalar &a, std::complex<double>b)
    {
      a.r = b.real();
      a.i = b.imag();
      return a;
    }

    inline SuperLuType<double>::Scalar to_superlu(SuperLuType<double>::Scalar &a, double b)
    {
      a = b;
      return a;
    }

    template<typename Scalar>
    SuperLUSolver<Scalar>::~SuperLUSolver()
    {
      _F_;
      free_factorization_data();
      free_matrix();
      free_rhs();

      if (local_Ai)  delete [] local_Ai;
      if (local_Ap)  delete [] local_Ap;
      if (local_Ax)  delete [] local_Ax;
      if (local_rhs) delete [] local_rhs;
    }

    template<typename Scalar>
    int SuperLUSolver<Scalar>::get_matrix_size()
    {
      return m->size;
    }

    template<typename Scalar>
    bool SuperLUSolver<Scalar>::solve()
    {
      _F_;
      assert(m != NULL);
      assert(rhs != NULL);

      Hermes::TimePeriod tmr;

      // Initialize the statistics variable.
      slu_stat_t stat;
      SLU_INIT_STAT(&stat);

      // Prepare data structures serving as input for the solver driver
      // (according to the chosen factorization reuse strategy).
      void *work = NULL;        // Explicit pointer to the factorization workspace
      // (unused, see below).
      int lwork = 0;            // Space for the factorization will be allocated
      // internally by system malloc.
      double ferr = 1.0;        // Estimated relative forward error
      // (unused unless iterative refinement is performed).
      double berr = 1.0;        // Estimated relative backward error
      // (unused unless iterative refinement is performed).
      slu_memusage_t memusage;  // Record the memory usage statistics.
      double rpivot_growth;     // The reciprocal pivot growth factor.
      double rcond;             // The estimate of the reciprocal condition number.
#ifdef SLU_MT
      options.work = work;
      options.lwork = lwork;
#endif

      if ( !setup_factorization() )
      {
        warning("LU factorization could not be completed.");
        return false;
      }

      // If the previous factorization of A is to be fully reused as an input for the solver driver,
      // keep the (possibly rescaled) matrix from the last factorization, otherwise recreate it
      // from the master SuperLUMatrix<Scalar> pointed to by this->m (this also applies to the case when
      // A does not yet exist).
      if (!has_A || this->factorization_scheme != HERMES_REUSE_FACTORIZATION_COMPLETELY)
      {
        if (A_changed)
          free_matrix();

        if (!has_A)
        {
          // A will be created from the local copy of the value and index arrays, because these
          // may be modified by the solver driver.
          if (local_Ai) delete [] local_Ai;
          local_Ai = new int [m->nnz];
          memcpy(local_Ai, m->Ai, m->nnz * sizeof(int));

          if (local_Ap) delete [] local_Ap;
          local_Ap = new int [m->size + 1];
          memcpy(local_Ap, m->Ap, (m->size + 1) * sizeof(int));

          if (local_Ax) delete [] local_Ax;
          local_Ax = new typename SuperLuType<Scalar>::Scalar[m->nnz];
          for (unsigned int i = 0;i<m->nnz;i++)
            to_superlu(local_Ax[i], m->Ax[i]);

          // Create new general (non-symmetric), column-major, non-supernodal, size X size matrix.
          create_csc_matrix(&A, m->size, m->size, m->nnz, local_Ax, local_Ai, local_Ap, SLU_NC, SLU_DTYPE, SLU_GE);

          has_A = true;
        }
      }

      // Recreate the input rhs for the solver driver from a local copy of the new value array.
      free_rhs();

      if (local_rhs) delete [] local_rhs;
      local_rhs = new typename SuperLuType<Scalar>::Scalar[rhs->size];
      for (unsigned int i = 0;i<rhs->size;i++)
        to_superlu(local_rhs[i], rhs->v[i]);

      create_dense_matrix(&B, rhs->size, 1, local_rhs, rhs->size, SLU_DN, SLU_DTYPE, SLU_GE);

      has_B = true;

      // Initialize the solution variable.
      SuperMatrix X;
      typename SuperLuType<Scalar>::Scalar *x;
      if ( !(x = new typename SuperLuType<Scalar>::Scalar[m->size]) )
        error("Malloc fails for x[].");
      create_dense_matrix(&X, m->size, 1, x, m->size, SLU_DN, SLU_DTYPE, SLU_GE);

      // Solve the system.
      int info;

#ifdef SLU_MT
      if (options.refact == NO)
      {
        // Get column permutation vector perm_c[], according to the first argument:
        //  0: natural ordering
        //  1: minimum degree ordering on structure of A'*A
        //  2: minimum degree ordering on structure of A' + A
        //  3: approximate minimum degree for unsymmetric matrices
        get_perm_c(1, &A, perm_c);
      }

      /*
      // Compute reciprocal pivot growth, estimate reciprocal condition number of A, solve,
      // perform iterative refinement of the solution and estimate forward and backward error.
      // Memory usage will be acquired at the end. If A is singular, info will be set to A->ncol + 1.
      //
      slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
      &L, &U, &B, &X, &rpivot_growth, &rcond, &ferr, &berr,
      &stat, &memusage, &info );
      */

      // ... OR ...

      // Estimate reciprocal condition number of A and solve the system. If A is singular, info
      // will be set to A->ncol + 1.
      //
      slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
        &L, &U, &B, &X, NULL, &rcond, NULL, NULL,
        &stat, NULL, &info );

      // ... OR ...

      /*
      // Do not check the regularity of A and just solve the system.
      //
      slu_mt_solver_driver( &options, &A, perm_c, perm_r, &AC, &equed, R, C,
      &L, &U, &B, &X, NULL, NULL, NULL, NULL,
      &stat, NULL, &info );
      */
#else
      solver_driver(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U,
        work, lwork, &B, &X, &rpivot_growth, &rcond, &ferr, &berr,
        &memusage, &stat, &info);
#endif

      // A and B may have been multiplied by the scaling vectors R and C on the output of the
      // solver. If the next call to the solver should reuse factorization only partially,
      // it will need the original unscaled matrix - this will indicate such situation
      // (rhs is always recreated anew).
#ifdef SLU_MT
      A_changed = (equed != NOEQUIL);
#else
      A_changed = (*equed != 'N');
#endif

      bool factorized = check_status(info);

      if (factorized)
      {
        delete [] this->sln;
        this->sln = new Scalar[m->size];

        Scalar *sol = (Scalar*) ((DNformat*) X.Store)->nzval;

        for (unsigned int i = 0; i < rhs->size; i++)
          this->sln[i] = sol[i];
      }

      // If required, print statistics.
      if ( options.PrintStat ) SLU_PRINT_STAT(&stat);

      // Free temporary local variables.
      StatFree(&stat);
      //SUPERLU_FREE (x);
      delete x;
      Destroy_SuperMatrix_Store(&X);

      tmr.tick();
      this->time = tmr.accumulated();

      return factorized;
    }

    template<typename Scalar>
    bool SuperLUSolver<Scalar>::setup_factorization()
    {
      _F_;
      unsigned int A_size = A.nrow < 0 ? 0 : A.nrow;
      if (has_A && this->factorization_scheme != HERMES_FACTORIZE_FROM_SCRATCH && A_size != m->size)
      {
        warning("You cannot reuse factorization structures for factorizing matrices of different sizes.");
        return false;
      }

      // Always factorize from scratch for the first time.
      int eff_fact_scheme;
      if (!inited)
        eff_fact_scheme = HERMES_FACTORIZE_FROM_SCRATCH;
      else
        eff_fact_scheme = this->factorization_scheme;

      // Prepare factorization structures. In case of a particular reuse scheme, comments are given
      // to clarify which arguments will be reused and which will be reset by the dgssvx (zgssvx) routine.
      // It was determined empirically by running the dlinsolx2 example from SuperLU, setting options.Fact
      // to the appropriate value and reallocating the various structures before the second run of dgssvx,
      // and observing when segfault will happen and when not. It is actually not needed to reallocate
      // the structures by hand, but comments at various places of the SuperLU 4.0 library contradict
      // each other and often lead to segfault when the structures are reallocated according to them.
      // It might thus bring some insight into how SuperLU works and how to correctly use it
      // (the PDF documentation is, unfortunately, even less helpful).
      switch (eff_fact_scheme)
      {
      case HERMES_FACTORIZE_FROM_SCRATCH:
        // This case should generally allow for solving a completely new system, i.e. for a change of
        // matrix and rhs size - for simplicity, we reallocate the structures every time.

        // Clear the structures emanating from previous factorization.
        free_factorization_data();

        // Allocate the row/column reordering vectors.
        if ( !(perm_c = intMalloc(m->size)) )
          error("Malloc fails for perm_c[].");
        if ( !(perm_r = intMalloc(m->size)) )
          error("Malloc fails for perm_r[].");

        // Allocate vectors with row/column scaling factors.
        if ( !(R = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) )
          error("SUPERLU_MALLOC fails for R[].");
        if ( !(C = (double *) SUPERLU_MALLOC(m->size * sizeof(double))) )
          error("SUPERLU_MALLOC fails for C[].");

#ifdef SLU_MT
        options.fact = EQUILIBRATE;
        options.refact = NO;
        options.perm_c = perm_c;
        options.perm_r = perm_r;
#else
        // Allocate additional structures used by the driver routine of sequential SuperLU.
        // Elimination tree is contained in the options structure in SuperLU_MT.
        if ( !(etree = intMalloc(m->size)) )
          error("Malloc fails for etree[].");

        options.Fact = DOFACT;
#endif
        A_changed = true;
        break;
      case HERMES_REUSE_MATRIX_REORDERING:
        // needed from previous:      etree, perm_c
        // not needed from previous:  perm_r, R, C, L, U, equed
#ifdef SLU_MT
        options.fact = EQUILIBRATE;
        options.refact = YES;
#else
        options.Fact = SamePattern;
#endif
        // L, U matrices may be reused without reallocating.
        // SLU_DESTROY_L(&L);
        // SLU_DESTROY_U(&U);
        break;
      case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
        // needed from previous:      etree, perm_c, perm_r, L, U
        // not needed from previous:  R, C, equed
#ifdef SLU_MT
        // MT version of SLU cannot reuse the equilibration factors (R, C), so
        // this is the same as the previous case.
        options.fact = EQUILIBRATE;
        options.refact = YES;
#else
        options.Fact = SamePattern_SameRowPerm;
#endif
        break;
      case HERMES_REUSE_FACTORIZATION_COMPLETELY:
        // needed from previous:      perm_c, perm_r, equed, L, U
        // not needed from previous:  etree, R, C
#ifdef SLU_MT
        options.fact = FACTORED;
        options.refact = YES;
#else
        options.Fact = FACTORED;
#endif
        break;
      }

      inited = true;

      return true;
    }

    template<typename Scalar>
    void SuperLUSolver<Scalar>::free_matrix()
    {
      _F_;
      if (has_A)
      {
        Destroy_SuperMatrix_Store(&A);
        has_A = false;
      }
    }

    template<typename Scalar>
    void SuperLUSolver<Scalar>::free_rhs()
    {
      _F_;
      if (has_B)
      {
        Destroy_SuperMatrix_Store(&B);
        has_B = false;
      }
    }

    template<typename Scalar>
    void SuperLUSolver<Scalar>::free_factorization_data()
    {
      _F_;
      if (inited)
      {
#ifdef SLU_MT
        SUPERLU_FREE(options.etree);
        SUPERLU_FREE(options.colcnt_h);
        SUPERLU_FREE(options.part_super_h);
        Destroy_CompCol_Permuted(&AC);
#else
        SUPERLU_FREE (etree);
#endif
        SUPERLU_FREE (perm_c);
        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (R);
        SUPERLU_FREE (C);
        SLU_DESTROY_L(&L);
        SLU_DESTROY_U(&U);
        inited = false;
      }
    }

#ifdef SLU_MT
    // This is a modification of the original p*gssvx routines from the SuperLU_MT library.
    //
    // The original routines have been changed in view of our applications, i.e.
    //  * only one right hand side is allowed,
    //  * some initial parameter checks have been omitted,
    //  * macros allowing abstraction from the fundamental Scalar datatype have been used
    //  * some phases of the calculation may be omitted for speed-up (less information about
    //    the matrix/solution can then be acquired, however),
    //  * deallocation at the end of the routine has been removed (this was neccessary to
    //    enable factorization reuse).
    //
    // See the correspondingly named attributes of SuperLUSolver class for brief description
    // of most parameters or the library source code for pdgssvx for more details. You may pass
    // NULL for
    //  * recip_pivot_growth  - reciprocal pivot growth factor will then not be computed;
    //                          reip_pivot_growth much less than one may indicate poor
    //                          stability of the factorization;
    //  * rcond               - estimate of the reciprocal condition number of matrix A will
    //                          then not be computed; this will prevent detection of singularity
    //                          of matrix A;
    //  * ferr or berr        - iterative refinement of the solution will then not be performed;
    //                          this also prevents computation of forward and backward error
    //                          estimates of the computed solution;
    //  * memusage            - memory usage during the factorization/solution will not be queried.
    //
    void slu_mt_solver_driver(slu_options_t *options, SuperMatrix *A,
      int *perm_c, int *perm_r, SuperMatrix *AC,
      equed_t *equed, double *R, double *C,
      SuperMatrix *L, SuperMatrix *U,
      SuperMatrix *B, SuperMatrix *X,
      double *recip_pivot_growth, double *rcond,
      double *ferr, double *berr,
      slu_stat_t *stat, slu_memusage_t *memusage,
      int *info)
    {
      /* Profiling variables. */
      double    t0;
      flops_t   flopcnt;

      /* State variables. */
      int dofact = (options->fact == DOFACT);
      int equil = (options->fact == EQUILIBRATE);
      int notran = (options->trans == NOTRANS);
      int colequ, rowequ;

      /* Right hand side and solution vectors. */
      DNformat *Bstore = (DNformat*) B->Store;
      DNformat *Xstore = (DNformat*) X->Store;
      Scalar *Bmat = (Scalar*) Bstore->nzval;
      Scalar *Xmat = (Scalar*) Xstore->nzval;

      *info = 0;

      /* ------------------------------------------------------------
      Diagonal scaling to equilibrate the matrix.
      ------------------------------------------------------------*/
      if (dofact || equil)
      {
        *equed = NOEQUIL;
        rowequ = colequ = FALSE;
      }
      else
      {
        rowequ = (*equed == ROW) || (*equed == BOTH);
        colequ = (*equed == COL) || (*equed == BOTH);
      }

      if ( equil )
      {
        t0 = SuperLU_timer_();
        /* Compute row and column scalings to equilibrate the matrix A. */
        int info1;
        double rowcnd, colcnd, amax;
        SLU_GSEQU(A, R, C, &rowcnd, &colcnd, &amax, &info1);

        if ( info1 == 0 )
        {
          /* Equilibrate matrix A. */
          SLU_LAQGS(A, R, C, rowcnd, colcnd, amax, equed);
          rowequ = (*equed == ROW) || (*equed == BOTH);
          colequ = (*equed == COL) || (*equed == BOTH);
        }
        stat->utime[EQUIL] = SuperLU_timer_() - t0;
      }

      /* ------------------------------------------------------------
      Scale the right hand side.
      ------------------------------------------------------------*/
      if ( notran )
      {
        if ( rowequ )
          for (int i = 0; i < A->nrow; ++i)
            SLU_MULT(Bmat[i], R[i]);
      }
      else if ( colequ )
      {
        for (int i = 0; i < A->nrow; ++i)
          SLU_MULT(Bmat[i], C[i]);
      }

      /* ------------------------------------------------------------
      Perform the LU factorization.
      ------------------------------------------------------------*/
      if ( dofact || equil )
      {
        /* Obtain column etree, the column count (colcnt_h) and supernode
        partition (part_super_h) for the Householder matrix. */
        if (options->refact == NO)
        {
          t0 = SuperLU_timer_();
          SLU_SP_COLORDER(A, perm_c, options, AC);
          stat->utime[ETREE] = SuperLU_timer_() - t0;
        }

        /* Compute the LU factorization of A*Pc. */
        t0 = SuperLU_timer_();
        SLU_GSTRF(options, AC, perm_r, L, U, stat, info);
        stat->utime[FACT] = SuperLU_timer_() - t0;

        flopcnt = 0;
        for (int i = 0; i < options->nprocs; ++i) flopcnt += stat->procstat[i].fcops;
        stat->ops[FACT] = flopcnt;

        if ( options->lwork == -1 )
        {
          if (memusage)
            memusage->total_needed = *info - A->ncol;
          return;
        }
      }

      if ( *info > 0 )
      {
        if ( *info <= A->ncol )
        {
          /* Compute the reciprocal pivot growth factor of the leading
          rank-deficient *info columns of A. */
          if (recip_pivot_growth)
            *recip_pivot_growth = SLU_PIVOT_GROWTH(*info, A, perm_c, L, U);
        }
      }
      else
      {
        /* ------------------------------------------------------------
        Compute the reciprocal pivot growth factor *recip_pivot_growth.
        ------------------------------------------------------------*/
        if (recip_pivot_growth)
          *recip_pivot_growth = SLU_PIVOT_GROWTH(A->ncol, A, perm_c, L, U);

        /* ------------------------------------------------------------
        Estimate the reciprocal of the condition number of A.
        ------------------------------------------------------------*/
        if (rcond)
        {
          t0 = SuperLU_timer_();

          // Next two lines are a bit complicated, but taken as they appear
          // in the original library function.
          char norm[1];
          *(unsigned char *)norm = (notran) ? '1' : 'I';

          double anorm = SLU_LANGS(norm, A);
          SLU_GSCON(norm, L, U, anorm, rcond, info);
          stat->utime[RCOND] = SuperLU_timer_() - t0;
        }

        /* ------------------------------------------------------------
        Compute the solution matrix X.
        ------------------------------------------------------------*/
        // Save a copy of the right hand side.
        memcpy(Xmat, Bmat, B->nrow * sizeof(Scalar));

        t0 = SuperLU_timer_();
        SLU_GSTRS(options->trans, L, U, perm_r, perm_c, X, stat, info);
        stat->utime[SOLVE] = SuperLU_timer_() - t0;
        stat->ops[SOLVE] = stat->ops[TRISOLVE];

        /* ------------------------------------------------------------
        Use iterative refinement to improve the computed solution and
        compute error bounds and backward error estimates for it.
        ------------------------------------------------------------*/
        if (ferr && berr)
        {
          t0 = SuperLU_timer_();
          SLU_GSRFS(options->trans, A, L, U, perm_r, perm_c, *equed,
            R, C, B, X, ferr, berr, stat, info);
          stat->utime[REFINE] = SuperLU_timer_() - t0;
        }

        /* ------------------------------------------------------------
        Transform the solution matrix X to a solution of the original
        system.
        ------------------------------------------------------------*/
        if ( notran )
        {
          if ( colequ )
            for (int i = 0; i < A->nrow; ++i)
              SLU_MULT(Xmat[i], C[i]);
        }
        else if ( rowequ )
        {
          for (int i = 0; i < A->nrow; ++i)
            SLU_MULT(Xmat[i], R[i]);
        }

        /* Set INFO = A->ncol + 1 if the matrix is singular to
        working precision.*/
        char param[1]; param[0] = 'E';
        if ( rcond && *rcond < SLU_LAMCH_(param) ) *info = A->ncol + 1;
      }

      if (memusage)
        SLU_QUERY_SPACE(options->nprocs, L, U, options->panel_size, memusage);
    }
#endif

    template class HERMES_API SuperLUSolver<double>;
    template class HERMES_API SuperLUSolver<std::complex<double> >;
  }
}
#endif
