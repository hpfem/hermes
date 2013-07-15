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
      virtual Scalar get(unsigned int idx) const;
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
        Jacobi = 0,
        MultiColoredSGS = 1,
        ILU = 2,
        MultiColoredILU = 3,
        IC = 4,
        AIChebyshev = 5
      };

      /// Constructor.
      /// \param[in] paralutionPrecondType The preconditioner type to create.
      ParalutionPrecond(ParalutionPreconditionerType paralutionPrecondType);
      
      paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>& get_paralutionPreconditioner();
      static paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* return_paralutionPreconditioner(ParalutionPreconditionerType type);
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

    /// \brief Encapsulation of PARALUTION iterative linear solver.
    ///
    /// @ingroup Solvers
    template <typename Scalar>
    class HERMES_API IterativeParalutionLinearMatrixSolver : public IterSolver<Scalar>
    {
    public:
      /// Constructor of UMFPack solver.
      /// @param[in] m pointer to matrix
      /// @param[in] rhs pointer to right hand side vector
      IterativeParalutionLinearMatrixSolver();
      IterativeParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *m, ParalutionVector<Scalar> *rhs);
      virtual ~IterativeParalutionLinearMatrixSolver();

      /// The solver type.
      /// Default: CG
      enum ParalutionSolverType
      {
        CG = 0,
        GMRES = 1,
        BiCGStab = 2
      };

      /// Set current solver type.
      /// This destroys the current solver (NOT the matrix, and rhs).
      void set_solver_type(ParalutionSolverType paralutionSolverType);

      virtual bool solve(Scalar* initial_guess);
      virtual bool solve();

      /// Get number of iterations.
      virtual int get_num_iters();

      /// Get the residual value.
      virtual double get_residual();

      virtual void set_precond(Precond<Scalar> *pc);

      /// Sets the verboseness.
      virtual void set_verbose_output(bool to_set);

      /// Utility.
      virtual int get_matrix_size();

      // Linear Solver creation.
      static paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* return_paralutionSolver(ParalutionSolverType type);
      
    private:
      /// Preconditioner.
      Preconditioners::ParalutionPrecond<Scalar> *preconditioner;

      /// Internal solver.
      paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* paralutionSolver;
      /// Internal solver is not reusable and will have to be re-created.
      void reset_internal_solver();

      /// Set internal solver for the current solution.
      void init_internal_solver();
      
      /// Matrix to solve.
      ParalutionMatrix<Scalar> *matrix;
      /// Right hand side vector.
      ParalutionVector<Scalar> *rhs;

      // Paralution solver type.
      ParalutionSolverType paralutionSolverType;

      // Store num_iters.
      int num_iters;

      // Store final_residual.
      double final_residual;

      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs, bool use_direct_solver);
      template<typename T> friend class AMGParalutionLinearMatrixSolver;
    };

    /// \brief Encapsulation of PARALUTION AMG linear solver.
    ///
    /// @ingroup Solvers
    template <typename Scalar>
    class HERMES_API AMGParalutionLinearMatrixSolver : public AMGSolver<Scalar>
    {
    public:
      /// Constructor of UMFPack solver.
      /// @param[in] m pointer to matrix
      /// @param[in] rhs pointer to right hand side vector
      AMGParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *m, ParalutionVector<Scalar> *rhs);
      virtual ~AMGParalutionLinearMatrixSolver();

      virtual bool solve(Scalar* initial_guess);
      virtual bool solve();

      /// Get number of iterations.
      virtual int get_num_iters();

      /// Get the residual value.
      virtual double get_residual();

      /// Set smoother (another PARALUTION linear matrix solver).
      virtual void set_smoother(typename IterativeParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType solverType, typename ParalutionPrecond<Scalar>::ParalutionPreconditionerType preconditionerType);

      /// Sets the verboseness.
      virtual void set_verbose_output(bool to_set);

      /// Utility.
      virtual int get_matrix_size();

    private:
      /// Smoother.
      typename IterativeParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType smootherSolverType;
      typename ParalutionPrecond<Scalar>::ParalutionPreconditionerType smootherPreconditionerType;

      /// Internal solver.
      paralution::AMG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* paralutionSolver;

      /// Set internal solver for the current solution.
      void init_internal_solver();
      /// Internal solver is not reusable and will have to be re-created.
      void reset_internal_solver();
      /// Matrix to solve.
      ParalutionMatrix<Scalar> *matrix;
      /// Right hand side vector.
      ParalutionVector<Scalar> *rhs;

      // Linear Solver creation.
      void init_paralutionSolver();
      
      // Store num_iters.
      int num_iters;

      // Store final_residual.
      double final_residual;

      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs, bool use_direct_solver);
    };
  }
}
#endif
#endif
