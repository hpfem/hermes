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
/*! \file matrix_solver.h
\brief General (linear/nonlinear) matrix solver functionality.
*/
#ifndef __HERMES_COMMON_MATRIX_SOLVER_H_
#define __HERMES_COMMON_MATRIX_SOLVER_H_

#include "interfaces/umfpack_solver.h"
#include "mixins.h"
#include "util/compat.h"

namespace Hermes
{
  namespace Solvers
  {
    /// \Basic interface for both linear and nonlinear algebraic solvers.
    template <typename Scalar>
    class HERMES_API MatrixSolver : public Hermes::Mixins::Loggable
    {
    public:
      MatrixSolver(bool force_use_direct_solver = false);
      virtual ~MatrixSolver();

      /// Data values (types) for UMFPACK reporting.
      enum UMFPACK_reporting_data_value
      {
        FactorizationSize = 0,
        PeakMemoryUsage = 1,
        Flops = 2
      };
      
       /// Return the solution vector.
      virtual Scalar *get_sln_vector();
     
      /// Sets the jacobian to be constant, i.e. reused whenever possible.
      void set_jacobian_constant(bool to_set = true);

      virtual Hermes::Solvers::LinearMatrixSolver<Scalar>* get_linear_matrix_solver();

      /// \TODO This is not used now.
      /// Set Reporting of UMFPACK numerical factorization data provided the used matrix solver is UMFPACK.
      virtual void set_UMFPACK_output(bool to_set = true, bool with_output = false);
      
      /// \TODO This is not used now.
      /// Get UMFPACK numerical factorization data provided the used matrix solver is UMFPACK
      virtual double get_UMFPACK_reporting_data(UMFPACK_reporting_data_value data_value);

      /// Only a shortcut for algebraic solver (->) get_matrix().
      SparseMatrix<Scalar>* get_jacobian();
      /// Only a shortcut for algebraic solver (->) get_rhs().
      Vector<Scalar>* get_residual();

      /// Verbose output.
      virtual void set_verbose_output(bool to_set);

    protected:
      /// Linear solver.
      Hermes::Solvers::LinearMatrixSolver<Scalar>* linear_matrix_solver;
      
      /// Jacobian can be reused if possible.
      bool constant_jacobian;

      /// Jacobian is ready to be reused if desirable.
      bool jacobian_reusable;
      
      /// Number of equations.
      int problem_size;

      /// \TODO This is not used now.
      /// Switch for UMFPACK reporting.
      bool do_UMFPACK_reporting;
      /// \TODO This is not used now.
      void handle_UMFPACK_reports();

      /// \TODO This is not used now.
      /// Data for UMFPACK reporting.
      double UMFPACK_reporting_data[3];
      
      /// The solution vector.
      Scalar* sln_vector;
    };
  }
}
#endif