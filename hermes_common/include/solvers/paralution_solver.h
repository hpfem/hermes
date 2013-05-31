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
/*! \file paralution_solver.h
\brief PARALUTION solver interface.
*/
#ifndef __HERMES_COMMON_PARALUTION_SOLVER_H_ 
#define __HERMES_COMMON_PARALUTION_SOLVER_H_
#include "config.h"
#ifdef WITH_PARALUTION
#include "linear_matrix_solver.h"
#include "matrix.h"

#include "paralution.hpp"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Algebra
  {
    using namespace Hermes::Solvers;

    enum ParalutionMatrixType
    {
      ParalutionMatrixTypeCSR,
      ParalutionMatrixTypeBCSR,
      ParalutionMatrixTypeMCSR,
      ParalutionMatrixTypeCOO,
      ParalutionMatrixTypeDIA,
      ParalutionMatrixTypeELL,
      ParalutionMatrixTypeHYB,
      ParalutionMatrixTypeDENSE
    };

    /// \brief General Paralution matrix.
    template <typename Scalar>
    class HERMES_API ParalutionMatrix : public SparseMatrix<Scalar>
    {
    public:
      /// \brief Default constructor.
      ParalutionMatrix(ParalutionMatrixType type = ParalutionMatrixTypeCSR);
      virtual ~ParalutionMatrix();

      virtual void alloc();
      virtual void free();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      

      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE, char* number_format = "%lf");
      
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;

      // Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);

    private:
      paralution::LocalMatrix<Scalar>* paralutionMatrix;
    };

    /// \brief Class representing the vector for UMFPACK.
    template <typename Scalar>
    class HERMES_API ParalutionVector : public Vector<Scalar>
    {
    public:
      ParalutionVector();
      virtual ~ParalutionVector();
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
    
    private:
      paralution::LocalVector<Scalar>* paralutionVector;
    };
  }
  namespace Solvers
  {
    /// \brief Encapsulation of PARALUTION linear solver.
    ///
    /// @ingroup Solvers
    template <typename Scalar>
    class HERMES_API ParalutionLinearMatrixSolver : public IterSolver<Scalar>
    {
    public:
      /// Constructor of UMFPack solver.
      /// @param[in] m pointer to matrix
      /// @param[in] rhs pointer to right hand side vector
      ParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *m, ParalutionVector<Scalar> *rhs);
      virtual ~ParalutionLinearMatrixSolver();

      virtual bool solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      ParalutionMatrix<Scalar> *m;
      /// Right hand side vector.
      ParalutionVector<Scalar> *rhs;

      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs);
    };
  }
}
#endif
#endif