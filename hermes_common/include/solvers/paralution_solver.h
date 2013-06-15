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
#include "cs_matrix.h"

#include "paralution.hpp"
#include "precond.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Algebra
  {
    using namespace Hermes::Solvers;

    /** Various types of matrix storage PARALUTION can do - commented out are those that are not implemented on our side so far */
    enum ParalutionMatrixType
    {
      ParalutionMatrixTypeCSR,
      // ParalutionMatrixTypeBCSR,
      // ParalutionMatrixTypeMCSR,
      // ParalutionMatrixTypeCOO,
      // ParalutionMatrixTypeDIA,
      // ParalutionMatrixTypeELL,
      // ParalutionMatrixTypeHYB,
      // ParalutionMatrixTypeDENSE
    };

    /// \brief General Paralution matrix.
    template <typename Scalar>
    class HERMES_API ParalutionMatrix : public CSRMatrix<Scalar>
    {
    public:
      /// \brief Default constructor.
      ParalutionMatrix(ParalutionMatrixType type = ParalutionMatrixTypeCSR);
      virtual ~ParalutionMatrix();

      virtual void free();
      virtual void zero();
      virtual void alloc();

      paralution::LocalMatrix<Scalar>& get_paralutionMatrix();

    private:
      paralution::LocalMatrix<Scalar> paralutionMatrix;
      ParalutionMatrixType paralutionMatrixType;

      // Friends.
      template<typename T> friend SparseMatrix<T>*  create_matrix();
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /// \brief Class representing the vector for UMFPACK.
    template <typename Scalar>
    class HERMES_API ParalutionVector : public Vector<Scalar>
    {
    public:
      /// Default constructor.
      ParalutionVector();
      /// Constructor of vector with specific size.
      /// @param[in] size size of vector
      ParalutionVector(unsigned int size);
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

      paralution::LocalVector<Scalar>& get_paralutionVector();

    private:
      Scalar *v;
      paralution::LocalVector<Scalar> paralutionVector;
    };
  }

  namespace Preconditioners
  {
    /// \brief A PARALUTION preconditioner.
    ///
    /// @ingroup preconds
    template <typename Scalar>
    class ParalutionPrecond : public Hermes::Preconditioners::Precond<Scalar>
    {
    public:
      /// The preconditioner type.
      enum ParalutionPreconditionerType
      {
        Jacobi,
        MultiColoredSGS,
        ILU,
        MultiColoredILU,
        IC,
        AIChebyshev
      };

      /// Constructor.
      /// \param[in] paralutionPrecondType The preconditioner type to create.
      ParalutionPrecond(ParalutionPreconditionerType paralutionPrecondType);
      
      paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>& get_paralutionPreconditioner();
    private:
      // Paralution preconditioner
      paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* paralutionPreconditioner;
    };
  }

  namespace Solvers
  {
    /// Utility class for PARALUTION initialization.
    /// Methods called from HermesCommonApi.
    class ParalutionInitialization
    {
    public:
      static void init_paralution()
      {
        paralution::init_paralution();
        paralution::set_omp_threads_paralution(HermesCommonApi.get_integral_param_value(numThreads));
        paralution::info_paralution();
      }

      static void deinit_paralution()
      {
        paralution::stop_paralution();
      }
    };

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

      /// The solver type.
      /// Default: CG
      enum ParalutionSolverType
      {
        CG,
        GMRES,
        BiCGStab,
        AMG
      };

      /// Set current solver type.
      /// This destroys the current solver (NOT the matrix, and rhs).
      void set_solver_type(ParalutionSolverType paralutionSolverType);

      virtual bool solve(Scalar* initial_guess);
      virtual bool solve();

      /// Set the convergence tolerance.
      /// @param[in] tolerance - the tolerance to set
      /// @param[in] toleranceType - the tolerance to set
      virtual void set_tolerance(double tolerance, typename IterSolver<Scalar>::ToleranceType toleranceType);

      /// Set maximum number of iterations to perform.
      /// @param[in] iters - number of iterations
      /// Default:100.
      virtual void set_max_iters(int iters);

      /// Get number of iterations.
      virtual int get_num_iters();

      /// Get the residual value.
      virtual double get_residual();

      virtual void set_precond(Precond<Scalar> *pc);

      /// Sets the verboseness.
      virtual void set_verbose_output(bool to_set);

      /// Utility.
      virtual int get_matrix_size();

    private:
      /// Preconditioner.
      Preconditioners::ParalutionPrecond<Scalar> *preconditioner;

      /// Matrix to solve.
      ParalutionMatrix<Scalar> *matrix;
      /// Right hand side vector.
      ParalutionVector<Scalar> *rhs;

      // Paralution solver type.
      ParalutionSolverType paralutionSolverType;

      // Linear Solver.
      paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* create_paralutionSolver();
      
      // Store num_iters.
      int num_iters;

      // Store final_residual.
      double final_residual;

      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs);
    };
  }
}
#endif
#endif