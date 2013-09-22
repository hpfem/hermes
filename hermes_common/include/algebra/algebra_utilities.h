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
/*! \file algebra_utilities.h
\brief Utilities for all Algebra code.
*/
#ifndef __HERMES_COMMON_ALGEBRA_UTILITIES_H
#define __HERMES_COMMON_ALGEBRA_UTILITIES_H

namespace Hermes
{
  enum MatrixSolverType
  {
    SOLVER_UMFPACK = 0,
    SOLVER_PARALUTION_ITERATIVE = 1,
    SOLVER_PARALUTION_AMG = 2,
    SOLVER_PETSC = 3,
    SOLVER_MUMPS = 4,
    SOLVER_SUPERLU = 5,
    SOLVER_AMESOS = 6,
    SOLVER_AZTECOO = 7,
    SOLVER_EXTERNAL = 8,
    SOLVER_EMPTY = 100
  };

  enum DirectMatrixSolverType
  {
    DIRECT_SOLVER_UMFPACK = 0,
    DIRECT_SOLVER_MUMPS = 4,
    DIRECT_SOLVER_SUPERLU = 5,
    DIRECT_SOLVER_AMESOS = 6,
    // Solver external is here, because direct solvers are used in projections.
    DIRECT_SOLVER_EXTERNAL = 8
  };

  enum IterativeMatrixSolverType
  {
    ITERATIVE_SOLVER_PARALUTION = 1,
    ITERATIVE_SOLVER_PETSC = 3,
    ITERATIVE_SOLVER_AZTECOO = 7
  };

  enum AMGMatrixSolverType
  {
    AMG_SOLVER_PARALUTION = 2
  };

  namespace Algebra
  {
    /// Format of file matrix and vector output
    enum MatrixExportFormat
    {
      /// \brief Plain ascii file
      /// lines contains row column and value
      EXPORT_FORMAT_PLAIN_ASCII = 1,
      /// Binary MATio format
      EXPORT_FORMAT_MATLAB_MATIO = 4,
      /// \brief Matrix Market which can be read by pysparse library
      EXPORT_FORMAT_MATRIX_MARKET = 3
    };
  }
}
#endif
