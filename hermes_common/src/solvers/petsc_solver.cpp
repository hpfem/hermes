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
/*! \file petsc_solver.cpp
\brief PETSc solver interface.
*/
#include "config.h"
#ifdef WITH_PETSC
#include "petsc_solver.h"
#include "trace.h"
#include "error.h"
#include "callstack.h"

/// \todo Check #ifdef WITH_MPI and use the parallel methods from PETSc accordingly.

namespace Hermes
{
  namespace Algebra
  {
    static int num_petsc_objects = 0;

    Hermes::Helpers::CommandLineArgs cmd_line_args;

    inline void vec_get_value(Vec x, PetscInt ni, const PetscInt ix[], std::complex<double> y[])
    {
      VecGetValues(x, ni, ix, y);
    }

    void vec_get_value(Vec x, PetscInt ni, const PetscInt ix[], double y[])
    {
      PetscScalar *py = new PetscScalar[ni];
      VecGetValues(x, ni, ix, py);
      for (int i = 0;i<ni;i++)y[i] = py[i].real();
      delete [] py;
    }

    int remove_petsc_object()
    {
      _F_;
      PetscTruth petsc_initialized, petsc_finalized;
      int ierr = PetscFinalized(&petsc_finalized); CHKERRQ(ierr);
      ierr = PetscInitialized(&petsc_initialized); CHKERRQ(ierr);
      if (petsc_finalized == PETSC_TRUE || petsc_initialized == PETSC_FALSE)
        // This should never happen here.
        return -1;

      if (--num_petsc_objects == 0)
      {
        int ierr = PetscFinalize();
        CHKERRQ(ierr);
        info("PETSc finalized. No more PETSc usage allowed until application restart.");
      }
    }

    int add_petsc_object()
    {
      _F_;
      int ierr;
      PetscTruth petsc_initialized, petsc_finalized;
      ierr = PetscFinalized(&petsc_finalized); CHKERRQ(ierr);

      if (petsc_finalized == PETSC_TRUE)
        error("PETSc cannot be used once it has been finalized. You must restart the application.");

      ierr = PetscInitialized(&petsc_initialized); CHKERRQ(ierr);

      if (petsc_initialized != PETSC_TRUE)
      {
        bool have_args = cmd_line_args.check();

        if (have_args)
          ierr = PetscInitialize(&cmd_line_args.get_argc(),
          &cmd_line_args.get_argv(),
          PETSC_NULL, PETSC_NULL);
        else
#ifdef WITH_MPI
          cmd_line_args.missing_error();
#else
          ierr = PetscInitializeNoArguments();
#endif

        CHKERRQ(ierr);
      }

      num_petsc_objects++;
    }

    template<typename Scalar>
    PetscMatrix<Scalar>::PetscMatrix()
    {
      _F_;
      inited = false;
      add_petsc_object();
    }

    template<typename Scalar>
    PetscMatrix<Scalar>::~PetscMatrix()
    {
      _F_;
      free();
      remove_petsc_object();
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::alloc()
    {
      _F_;
      assert(this->pages != NULL);

      // calc nnz
      int *nnz_array = new int[this->size];
      MEM_CHECK(nnz_array);

      // fill in nnz_array
      int aisize = this->get_num_indices();
      int *ai = new int[aisize];
      MEM_CHECK(ai);

      // sort the indices and remove duplicities, insert into ai
      int pos = 0;
      for (unsigned int i = 0; i < this->size; i++)
      {
        nnz_array[i] = sort_and_store_indices(this->pages[i], ai + pos, ai + aisize);
        pos += nnz_array[i];
      }
      // stote the number of nonzeros
      nnz = pos;
      delete [] this->pages; this->pages = NULL;
      delete [] ai;

      //
      MatCreateSeqAIJ(PETSC_COMM_SELF, this->size, this->size, 0, nnz_array, &matrix);
      //	MatSetOption(matrix, MAT_ROW_ORIENTED);
      //	MatSetOption(matrix, MAT_ROWS_SORTED);

      delete [] nnz_array;

      inited = true;
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::free()
    {
      _F_;
      if (inited) MatDestroy(matrix);
      inited = false;
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::finish()
    {
      _F_;
      MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    }

    template<>
    double PetscMatrix<double>::get(unsigned int m, unsigned int n)
    {
      _F_;
      double v = 0.0;
      PetscScalar pv;
      MatGetValues(matrix, 1, (PetscInt*) &m, 1, (PetscInt*) &n, &pv);
      v = pv.real();
      return v;
    }

    template<>
    std::complex<double> PetscMatrix<std::complex<double> >::get(unsigned int m, unsigned int n)
    {
      _F_;
      std::complex<double> v = 0.0;
      MatGetValues(matrix, 1, (PetscInt*) &m, 1, (PetscInt*) &n, &v);
      return v;
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::zero()
    {
      _F_;
      MatZeroEntries(matrix);
    }

    inline PetscScalar to_petsc(double x)
    {
      return std::complex<double>(x, 0);
    }

    inline PetscScalar to_petsc(std::complex<double> x)
    {
      return x;
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar v)
    {
      _F_;
      if (v != 0.0)
      {		// ignore zero values.
        MatSetValue(matrix, (PetscInt) m, (PetscInt) n, to_petsc(v), ADD_VALUES);
      }
    }

    /// Add a number to each diagonal entry.

    template<typename Scalar>
    void PetscMatrix<Scalar>::add_to_diagonal(Scalar v)
    {
      for (unsigned int i = 0; i<this->size; i++)
      {
        add(i, i, v);
      }
    };

    template<typename Scalar>
    void PetscMatrix<Scalar>::add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols)
    {
      _F_;
      /// \todo pass in just the block of the matrix without HERMES_DIRICHLET_DOFs (so that can use MatSetValues directly without checking
      // row and cols for -1)
      for (unsigned int i = 0; i < m; i++)				// rows
        for (unsigned int j = 0; j < n; j++)			// cols
          if(rows[i] >= 0 && cols[j] >= 0) // not Dir. dofs.
            add(rows[i], cols[j], mat[i][j]);
    }

    template<typename Scalar>
    bool PetscMatrix<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      switch (fmt)
      {
      case DF_MATLAB_SPARSE: //only to stdout
        PetscViewer  viewer = PETSC_VIEWER_STDOUT_SELF;
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView(matrix, viewer);
        return true;
      }
      return false;
    }

    template<typename Scalar>
    unsigned int PetscMatrix<Scalar>::get_matrix_size() const
    {
      _F_;
      return this->size;
    }

    template<typename Scalar>
    unsigned int PetscMatrix<Scalar>::get_nnz() const
    {
      _F_;
      return nnz;
    }

    template<typename Scalar>
    double PetscMatrix<Scalar>::get_fill_in() const
    {
      _F_;
      return (double) nnz / ((double)this->size*this->size);
    }


    template<typename Scalar>
    void PetscMatrix<Scalar>::multiply_with_vector(Scalar* vector_in, Scalar* vector_out)
    {
      for (unsigned int i = 0;i<this->size;i++)
      {
        vector_out[i] = 0;
        for (unsigned int j = 0;j<this->size;j++)
        {
          vector_out[i] +=vector_in[j]*get(i, j);
        }
      }
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::add_matrix(PetscMatrix<Scalar>* mat)
    {
      MatAXPY(matrix, 1, mat->matrix, DIFFERENT_NONZERO_PATTERN);    //matrix = 1*mat + matrix (matrix and mat have different nonzero structure)
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::add_to_diagonal_blocks(int num_stages, PetscMatrix<Scalar>* mat)
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
    void PetscMatrix<Scalar>::add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat)
    {
      add_to_diagonal_blocks(num_stages, dynamic_cast<PetscMatrix<Scalar>*>(mat));
    }

    template<typename Scalar>
    void PetscMatrix<Scalar>::add_as_block(unsigned int i, unsigned int j, PetscMatrix<Scalar>* mat)
    {
      _F_;
      if ((this->get_size() < i + mat->get_size() )||(this->get_size() < j + mat->get_size() ))
        error("Incompatible matrix sizes in PetscMatrix<Scalar>::add_as_block()");
      unsigned int block_size = mat->get_size();
      for (unsigned int r = 0;r<block_size;r++)
      {
        for (unsigned int c = 0;c<block_size;c++)
        {
          this->add(i + r, j + c, mat->get(r, c));
        }
      }
    }

    // Multiplies matrix with a Scalar.

    template<typename Scalar>
    void PetscMatrix<Scalar>::multiply_with_Scalar(Scalar value)
    {
      _F_;
      MatScale(matrix, to_petsc(value));
    }
    // Creates matrix in PETSC format using size, nnz, and the three arrays.

    template<typename Scalar>
    void PetscMatrix<Scalar>::create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax)
    {
      _F_;
      this->size = size;
      this->nnz = nnz;
      PetscScalar* pax = new PetscScalar[nnz];
      for (unsigned i = 0;i<nnz;i++)
        pax[i] = to_petsc(ax[i]);
      MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size, size, ap, ai, pax, &matrix);
      delete pax;
    }

    template<typename Scalar>
    PetscMatrix<Scalar>* PetscMatrix<Scalar>::duplicate()
    {
      _F_;
      PetscMatrix<Scalar>*ptscmatrix = new PetscMatrix<Scalar>();
      MatDuplicate(matrix, MAT_COPY_VALUES, &(ptscmatrix->matrix));
      ptscmatrix->size = this->size;
      ptscmatrix->nnz = nnz;
      return ptscmatrix;
    };

    template<typename Scalar>
    PetscVector<Scalar>::PetscVector()
    {
      _F_;
      inited = false;
      add_petsc_object();
    }

    template<typename Scalar>
    PetscVector<Scalar>::~PetscVector()
    {
      _F_;
      free();
      remove_petsc_object();
    }

    template<typename Scalar>
    void PetscVector<Scalar>::alloc(unsigned int n)
    {
      _F_;
      free();
      this->size = n;
      VecCreateSeq(PETSC_COMM_SELF, this->size, &vec);
      inited = true;
    }

    template<typename Scalar>
    void PetscVector<Scalar>::free()
    {
      _F_;
      if (inited) VecDestroy(vec);
      inited = false;
    }

    template<typename Scalar>
    void PetscVector<Scalar>::finish()
    {
      _F_;
      VecAssemblyBegin(vec);
      VecAssemblyEnd(vec);
    }

    template<>
    double PetscVector<double>::get(unsigned int idx)
    {
      _F_;
      double y = 0;
      PetscScalar py;
      VecGetValues(vec, 1, (PetscInt*) &idx, &py);
      y = py.real();
      return y;
    }

    template<>
    std::complex<double> PetscVector<std::complex<double> >::get(unsigned int idx)
    {
      _F_;
      std::complex<double> y = 0;
      VecGetValues(vec, 1, (PetscInt*) &idx, &y);
      return y;
    }

    template<typename Scalar>
    void PetscVector<Scalar>::extract(Scalar *v) const
    {
      _F_;
      int *idx = new int [this->size];
      for (unsigned int i = 0; i < this->size; i++) idx[i] = i;
      vec_get_value(vec, this->size, idx, v);
      delete [] idx;
    }

    template<typename Scalar>
    void PetscVector<Scalar>::zero()
    {
      _F_;
      VecZeroEntries(vec);
    }

    template<typename Scalar>
    void PetscVector<Scalar>::change_sign()
    {
      _F_;
      PetscScalar* y = new PetscScalar [this->size];
      int *idx = new int [this->size];
      for (unsigned int i = 0; i < this->size; i++) idx[i] = i;
      VecGetValues(vec, this->size, idx, y);
      for (unsigned int i = 0; i < this->size; i++) y[i] *= -1.;
      VecSetValues(vec, this->size, idx, y, INSERT_VALUES);
      delete [] y;
      delete [] idx;
    }

    template<typename Scalar>
    void PetscVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      _F_;
      VecSetValue(vec, idx, to_petsc(y), INSERT_VALUES);
    }

    template<typename Scalar>
    void PetscVector<Scalar>::add(unsigned int idx, Scalar y)
    {
      _F_;
      VecSetValue(vec, idx, to_petsc(y), ADD_VALUES);
    }

    template<typename Scalar>
    void PetscVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      _F_;
      PetscScalar py;
      for (unsigned int i = 0; i < n; i++)
      {
        VecSetValue(vec, idx[i], to_petsc(y[i]), ADD_VALUES);
      }
    }

    template<typename Scalar>
    void PetscVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec->get(i));
    }

    template<typename Scalar>
    void PetscVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec[i]);
    }

    template<typename Scalar>
    bool PetscVector<Scalar>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
    {
      _F_;
      switch (fmt)
      {
      case DF_MATLAB_SPARSE: //only to stdout
        PetscViewer  viewer = PETSC_VIEWER_STDOUT_SELF;
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        VecView(vec, viewer);
        return true;
      }
      return false;
    }

    template class HERMES_API PetscMatrix<double>;
    template class HERMES_API PetscMatrix<std::complex<double> >;
    template class HERMES_API PetscVector<double>;
    template class HERMES_API PetscVector<std::complex<double> >;
  }
  namespace Solvers
  {
    template<typename Scalar>
    PetscLinearMatrixSolver<Scalar>::PetscLinearMatrixSolver(PetscMatrix<Scalar> *mat, PetscVector<Scalar> *rhs)
      : DirectSolver<Scalar>(), m(mat), rhs(rhs)
    {
      _F_;
      add_petsc_object();
    }

    template<typename Scalar>
    PetscLinearMatrixSolver<Scalar>::~PetscLinearMatrixSolver()
    {
      _F_;
      remove_petsc_object();
    }

    template<typename Scalar>
    bool PetscLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return m->size();
    }

    template<typename Scalar>
    bool PetscLinearMatrixSolver<Scalar>::solve()
    {
      _F_;
      assert(m != NULL);
      assert(rhs != NULL);

      PetscErrorCode ec;
      KSP ksp;
      Vec x;

      Hermes::TimePeriod tmr;

      KSPCreate(PETSC_COMM_WORLD, &ksp);

      KSPSetOperators(ksp, m->matrix, m->matrix, DIFFERENT_NONZERO_PATTERN);
      KSPSetFromOptions(ksp);
      VecDuplicate(rhs->vec, &x);

      ec = KSPSolve(ksp, rhs->vec, x);
      if (ec) return false;

      tmr.tick();
      this->time = tmr.accumulated();

      // allocate memory for solution vector
      delete [] this->sln;
      this->sln = new Scalar [m->size];
      MEM_CHECK(this->sln);
      memset(this->sln, 0, m->size * sizeof(Scalar));

      // index map vector (basic serial code uses the map sln[i] = x[i] for all dofs.
      int *idx = new int [m->size];
      MEM_CHECK(idx);
      for (unsigned int i = 0; i < m->size; i++) idx[i] = i;

      // copy solution to the output solution vector
      vec_get_value(x, m->size, idx, this->sln);
      delete [] idx;

      KSPDestroy(ksp);
      VecDestroy(x);

      return true;
    }

    template class HERMES_API PetscLinearMatrixSolver<double>;
    template class HERMES_API PetscLinearMatrixSolver<std::complex<double> >;
  }
}
#endif
