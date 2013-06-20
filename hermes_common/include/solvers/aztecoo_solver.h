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
/*! \file aztecoo_solver.h
\brief AztecOOSolver class as an interface to AztecOO.
*/
#ifndef __HERMES_COMMON_AZTECOO_SOLVER_H_
#define __HERMES_COMMON_AZTECOO_SOLVER_H_
#include "../config.h"
#ifdef HAVE_AZTECOO
#include "epetra.h"
#include "linear_matrix_solver.h"
#include "precond_ifpack.h"
#include <AztecOO.h>

namespace Hermes
{
  namespace Solvers
  {
    /// \brief Encapsulation of AztecOO linear solver.
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API AztecOOSolver : public IterSolver<Scalar>
    {
    public:

      /// Set the type of the solver
      /// @param[in] solver - name of the solver[ gmres | cg | cgs | tfqmr | bicgstab ]
      void set_solver(const char *solver);
      /// Set the convergence tolerance
      /// @param[in] tol - the tolerance to set
      void set_tolerance(double tol);
      /// Set maximum number of iterations to perform
      /// @param[in] iters - number of iterations
      void set_max_iters(int iters);

      /// Set Aztec internal preconditioner
      /// @param[in] name - name of the preconditioner[ none | jacobi | neumann | least-squares ]
      virtual void set_precond(const char *name);
      
      /// \brief Set preconditioner from IFPACK or ML.
      /// @param[in] pc - IFPACK or ML preconditioner
      virtual void set_precond(Precond<Scalar> *pc);

      AztecOOSolver(EpetraMatrix<Scalar> *m, EpetraVector<Scalar> *rhs);
      virtual ~AztecOOSolver();
      
      virtual bool solve();
      virtual bool solve(Scalar* initial_guess);
      
      virtual int get_matrix_size();
    
      virtual int get_num_iters();
      virtual double get_residual();

      /// Option setting function
      void set_option(int option, int value);

      /// Parameter setting function
      void set_param(int param, double value);

    protected:
      AztecOO aztec;    ///< Instance of the Aztec solver.
      EpetraMatrix<Scalar> *m;
      EpetraVector<Scalar> *rhs;

      EpetraPrecond<Scalar> *pc;
      
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs, bool use_direct_solver = false);
    };
  }
}
#endif
#endif