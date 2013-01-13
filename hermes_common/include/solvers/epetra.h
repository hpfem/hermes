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
/*! \file epetra.h
\brief EpetraMatrix and EpetraVector storage classes for Amesos, AztecOO, ... .
*/
#ifndef __HERMES_COMMON_SOLVER_EPETRA_H_
#define __HERMES_COMMON_SOLVER_EPETRA_H_
#include "config.h"
#ifdef HAVE_EPETRA
#define EPETRA_NO_64BIT_GLOBAL_INDICES
#include "matrix.h"
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class AmesosSolver;
    template <typename Scalar> class AztecOOSolver;
    template <typename Scalar> class NewtonSolverNOX;
  }
  namespace Preconditioners
  {
    template <typename Scalar> class IfpackPrecond;
    template <typename Scalar> class MlPrecond;
  }
}

namespace Hermes
{
  namespace Algebra
  {
    template <typename Scalar>
    class HERMES_API EpetraMatrix : public SparseMatrix<Scalar>
    {
    public:
      EpetraMatrix();
      EpetraMatrix(Epetra_RowMatrix &mat);
      virtual ~EpetraMatrix();

      virtual void prealloc(unsigned int n);
      virtual void pre_add_ij(unsigned int row, unsigned int col);
      virtual void finish();

      virtual void alloc();
      virtual void free();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual int get_num_row_entries(unsigned int row);
      virtual void extract_row_copy(unsigned int row, unsigned int len, unsigned int &n_entries, double *vals, unsigned int *idxs);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      virtual void add_to_diagonal_blocks(int num_stages, EpetraMatrix<Scalar>* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      virtual void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      virtual void add_as_block(unsigned int i, unsigned int j, EpetraMatrix<Scalar>* mat);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;

    protected:
      Epetra_BlockMap *std_map;
      Epetra_CrsGraph *grph;
      Epetra_CrsMatrix *mat;
      /// \brief Imaginary part of the matrix, mat holds the real part.
      Epetra_CrsMatrix *mat_im;
      bool owner;

      friend class Hermes::Solvers::AmesosSolver<Scalar>;
      friend class Hermes::Solvers::AztecOOSolver<Scalar>;
      friend class Hermes::Solvers::NewtonSolverNOX<Scalar>;
      friend class Hermes::Preconditioners::IfpackPrecond<Scalar>;
      friend class Hermes::Preconditioners::MlPrecond<Scalar>;
    };

    template <typename Scalar>
    class HERMES_API EpetraVector : public Vector<Scalar>
    {
    public:
      EpetraVector();
      EpetraVector(const Epetra_Vector &v);
      virtual ~EpetraVector();

      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual Scalar get(unsigned int idx);
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual void change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual void add_vector(Vector<Scalar>* vec);
      virtual void add_vector(Scalar* vec);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");

    protected:
      Epetra_BlockMap *std_map;
      Epetra_Vector *vec;
      /// \brief Imaginary part of the vector, vec holds the real part.
      Epetra_Vector *vec_im;
      bool owner;

      friend class Hermes::Solvers::AmesosSolver<Scalar>;
      friend class Hermes::Solvers::AztecOOSolver<Scalar>;
      friend class Hermes::Solvers::NewtonSolverNOX<Scalar>;
    };
  }
}
#endif
#endif
