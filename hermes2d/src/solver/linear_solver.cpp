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
/*! \file linear_solver.cpp
\brief General linear solver functionality.
*/
#include "solver/linear_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver() : Solver()
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(DiscreteProblem<Scalar>* dp) : Solver(dp)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space) : Solver(wf, space)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::LinearSolver(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) : Solver(wf, spaces)
    {
      this->init_linear();
    }

    template<typename Scalar>
    LinearSolver<Scalar>::~LinearSolver()
    {
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::init_linear()
    {
      this->dp->nonlinear = false;
    }

    template<typename Scalar>
    bool LinearSolver<Scalar>::isOkay() const
    {
      return this->dp->isOkay();
    }

    template<typename Scalar>
    void LinearSolver<Scalar>::solve()
    {
      this->check();

      this->tick();

      this->on_initialization();

      dp->assemble(this->jacobian, this->residual);
      if(this->output_rhsOn && (this->output_rhsIterations == -1 || this->output_rhsIterations >= 1))
      {
        char* fileName = new char[this->RhsFilename.length() + 5];
        if(this->RhsFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
          sprintf(fileName, "%s%i.m", this->RhsFilename.c_str(), 1);
        else
          sprintf(fileName, "%s%i", this->RhsFilename.c_str(), 1);
        FILE* rhs_file = fopen(fileName, "w+");
        residual->dump(rhs_file, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
        fclose(rhs_file);
      }
      if(this->output_matrixOn && (this->output_matrixIterations == -1 || this->output_matrixIterations >= 1))
        {
          char* fileName = new char[this->matrixFilename.length() + 5];
          if(this->matrixFormat == Hermes::Algebra::DF_MATLAB_SPARSE)
            sprintf(fileName, "%s%i.m", this->matrixFilename.c_str(), 1);
          else
            sprintf(fileName, "%s%i", this->matrixFilename.c_str(), 1);
          FILE* matrix_file = fopen(fileName, "w+");

          jacobian->dump(matrix_file, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
          fclose(matrix_file);
        }

      this->matrix_solver->solve();

      this->sln_vector = matrix_solver->get_sln_vector();

      this->on_finish();
      
      this->tick();
      this->info("\tLinear solver solution duration: %f s.\n", this->last());
    }

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
  }
}
