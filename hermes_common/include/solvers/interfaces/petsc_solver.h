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
/*! \file petsc_solver.h
\brief PETSc solver interface.
*/
#ifndef __HERMES_COMMON_PETSC_SOLVER_H_
#define __HERMES_COMMON_PETSC_SOLVER_H_

#include "config.h"

#include "algebra/matrix.h"
#include "solvers/linear_matrix_solver.h"

#ifdef WITH_PETSC

#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>


namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class PetscLinearMatrixSolver;
  }
}

namespace Hermes
{
  namespace Algebra
  {
    /// \brief Wrapper of PETSc matrix, to store matrices used with PETSc in its native format.
    template <typename Scalar>
    class PetscMatrix : public SparseMatrix<Scalar>
    {
    public:
      PetscMatrix();
      virtual ~PetscMatrix();

      virtual void alloc();
      virtual void free();
      virtual void finish();
      virtual Scalar get(unsigned int m, unsigned int n) const;
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      using Matrix<Scalar>::export_to_file;
      virtual void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;

      /// Add matrix - General.
      /// See add_petsc_matrix
      /// @param[in] mat matrix to be added.
      virtual void add_sparse_matrix(SparseMatrix<Scalar>* mat);

      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);
      /// Creates matrix in PETSC format using size, nnz, and the three arrays.
      /// @param[in] size size of matrix (num of rows and columns)
      /// @param[in] nnz number of nonzero values
      /// @param[in] ap row indices
      /// @param[in] ai column indices
      /// @param[in] ax values
      ///  @todo same input parameters acts differen as in superlu
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);
      // Duplicates a matrix (including allocation).
      PetscMatrix* duplicate() const;
    protected:
      /// Add matrix - PETSc.
      /// @param[in] mat matrix to be added
      virtual void add_petsc_matrix(PetscMatrix* mat);

      /// Petsc matrix data structure.
      Mat matrix;
      /// Number of nonzero values.
      unsigned int nnz;
      /// Is matrix inited (allocated)?
      bool inited;

      friend class Solvers::PetscLinearMatrixSolver<Scalar>;
    };

    /// Wrapper of PETSc vector, to store vectors used with PETSc in its native format.
    ///
    template <typename Scalar>
    class PetscVector : public Vector<Scalar>
    {
    public:
      PetscVector();
      virtual ~PetscVector();

      virtual void alloc(unsigned int ndofs);
      virtual void free();
      /// Finish manipulation with vector.
      virtual void finish();
      virtual Scalar get(unsigned int idx) const;
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual Vector<Scalar>* change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual Vector<Scalar>* add_vector(Vector<Scalar>* vec);
      virtual Vector<Scalar>* add_vector(Scalar* vec);
      using Vector<Scalar>::export_to_file;
      virtual void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");

    protected:
      /// Petsc vector data structure.
      Vec vec;
      /// Is vector initiated (allocated)?
      bool inited;

      friend class Solvers::PetscLinearMatrixSolver<Scalar>;
    };
  }
  namespace Solvers
  {
    /// Encapsulation of PETSc linear solver.
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API PetscLinearMatrixSolver : public DirectSolver<Scalar>
    {
    public:
      PetscLinearMatrixSolver(PetscMatrix<Scalar> *mat, PetscVector<Scalar> *rhs);
      virtual ~PetscLinearMatrixSolver();

      virtual void solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      PetscMatrix<Scalar> *m;
      /// Right hand side vector.
      PetscVector<Scalar> *rhs;
    };
  }
}
#endif
#endif
