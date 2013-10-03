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

#include "matrix_solver.h"

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    MatrixSolver<Scalar>::MatrixSolver(bool force_use_direct_solver) : Hermes::Mixins::Loggable(true)
    {
      /// Linear solver.
      SparseMatrix<Scalar>* A = create_matrix<Scalar>(force_use_direct_solver);
      Vector<Scalar>* b = create_vector<Scalar>(force_use_direct_solver);
      this->linear_matrix_solver = create_linear_solver<Scalar>(A, b, force_use_direct_solver);

      this->constant_jacobian = false;

      this->jacobian_reusable = false;

      this->problem_size = -1;

      this->do_UMFPACK_reporting = false;

      this->sln_vector = nullptr;
    }

    template<typename Scalar>
    MatrixSolver<Scalar>::~MatrixSolver()
    {
      SparseMatrix<Scalar>* temp_matrix = linear_matrix_solver->get_matrix();
      Vector<Scalar>* temp_rhs = linear_matrix_solver->get_rhs();

      delete this->linear_matrix_solver;

      if(temp_matrix)
        delete temp_matrix;
      if(temp_rhs)
        delete temp_rhs;
    }

    template<typename Scalar>
    void MatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
      this->linear_matrix_solver->set_verbose_output(to_set);
    }

    template<typename Scalar>
    Hermes::Solvers::LinearMatrixSolver<Scalar>* MatrixSolver<Scalar>::get_linear_matrix_solver()
    {
      return this->linear_matrix_solver;
    }

    template<typename Scalar>
    SparseMatrix<Scalar>* MatrixSolver<Scalar>::get_jacobian()
    {
      return this->linear_matrix_solver->get_matrix();
    }

    template<typename Scalar>
    Vector<Scalar>* MatrixSolver<Scalar>::get_residual()
    {
      return this->linear_matrix_solver->get_rhs();
    }

    template<typename Scalar>
    Scalar* MatrixSolver<Scalar>::get_sln_vector()
    {
      return sln_vector;
    }

    template<typename Scalar>
    void MatrixSolver<Scalar>::set_jacobian_constant(bool to_set)
    {
      this->constant_jacobian = to_set;
      if(!to_set)
        this->jacobian_reusable = false;
    }

    template<typename Scalar>
    double MatrixSolver<Scalar>::get_UMFPACK_reporting_data(UMFPACK_reporting_data_value data_value)
    {
      return this->UMFPACK_reporting_data[data_value];
    }

    template<typename Scalar>
    void MatrixSolver<Scalar>::set_UMFPACK_output(bool to_set, bool with_output)
    {
      if(!dynamic_cast<UMFPackLinearMatrixSolver<Scalar>*>(this->linear_matrix_solver))
      {
        this->warn("A different MatrixSolver than UMFPACK is used, ignoring the call to set_UMFPACK_reporting().");
        return;
      }

      if(with_output)
        ((UMFPackLinearMatrixSolver<Scalar>*)this->linear_matrix_solver)->set_output_level(2);
      else
        ((UMFPackLinearMatrixSolver<Scalar>*)this->linear_matrix_solver)->set_output_level(0);

      this->do_UMFPACK_reporting = to_set;
    }

    template<typename Scalar>
    void MatrixSolver<Scalar>::handle_UMFPACK_reports()
    {
      if(this->do_UMFPACK_reporting)
      {
        UMFPackLinearMatrixSolver<Scalar>* umfpack_matrix_solver = (UMFPackLinearMatrixSolver<Scalar>*)this->linear_matrix_solver;
        if(this->linear_matrix_solver->get_used_reuse_scheme() != HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY)
        {
          this->UMFPACK_reporting_data[this->FactorizationSize] = umfpack_matrix_solver->Info[UMFPACK_NUMERIC_SIZE] * umfpack_matrix_solver->Info[UMFPACK_SIZE_OF_UNIT];
          this->UMFPACK_reporting_data[this->PeakMemoryUsage] = umfpack_matrix_solver->Info[UMFPACK_PEAK_MEMORY] * umfpack_matrix_solver->Info[UMFPACK_SIZE_OF_UNIT];
          this->UMFPACK_reporting_data[this->Flops] = umfpack_matrix_solver->Info[UMFPACK_FLOPS];
        }
        else
          memset(this->UMFPACK_reporting_data, 0, 3 * sizeof(double));
      }
    }

    template class HERMES_API MatrixSolver<double>;
    template class HERMES_API MatrixSolver<std::complex<double> >;
  }
}