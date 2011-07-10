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

#include "../matrix.h"
#include "linear_solver.h"

#ifdef WITH_PETSC
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

namespace Hermes {
  namespace Solvers {
    template <typename Scalar> class PetscLinearSolver;
  }
}

namespace Hermes {
  namespace Algebra {
    /// \brief Wrapper of PETSc matrix, to store matrices used with PETSc in its native format.
    template <typename Scalar>
    class PetscMatrix : public SparseMatrix<Scalar> {
    public:
      PetscMatrix();
      virtual ~PetscMatrix();

      virtual void alloc();
      virtual void free();
      virtual void finish();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;
      virtual void add_matrix(PetscMatrix* mat);
      virtual void add_to_diagonal_blocks(int num_stages, PetscMatrix* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat){
        add_to_diagonal_blocks(num_stages,dynamic_cast<PetscMatrix<Scalar>*>(mat));
      }
      virtual void add_as_block(unsigned int i, unsigned int j, PetscMatrix* mat);

      // Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);

      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);
      // Creates matrix in PETSC format using size, nnz, and the three arrays.
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);
      // Duplicates a matrix (including allocation).
      PetscMatrix* duplicate();
    protected:
      Mat matrix;
      unsigned int nnz;
      bool inited;

      friend class Solvers::PetscLinearSolver<Scalar>;
    };

    /// Wrapper of PETSc vector, to store vectors used with PETSc in its native format
    ///
    template <typename Scalar>
    class PetscVector : public Vector<Scalar> {
    public:
      PetscVector();
      virtual ~PetscVector();

      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual void finish();
      virtual Scalar get(unsigned int idx);
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual void change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual void add_vector(Vector<Scalar>* vec) {
        assert(this->length() == vec->length());
        for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec->get(i));
      };
      virtual void add_vector(Scalar* vec) {
        for (unsigned int i = 0; i < this->length(); i++) this->add(i, vec[i]);
      };
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

    protected:
      Vec vec;
      bool inited;

      friend class Solvers::PetscLinearSolver<Scalar>;
    };
  }
  namespace Solvers
  {
    /// Encapsulation of PETSc linear solver
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API PetscLinearSolver : public DirectSolver<Scalar> {
    public:
      PetscLinearSolver(PetscMatrix<Scalar> *mat, PetscVector<Scalar> *rhs);
      virtual ~PetscLinearSolver();

      virtual bool solve();

    protected:
      PetscMatrix<Scalar> *m;
      PetscVector<Scalar> *rhs;
    };
  }
}
#endif
#endif
