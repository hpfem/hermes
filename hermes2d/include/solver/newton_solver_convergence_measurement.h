// This file is part of Hermes2D
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
/*! \file solver_newton.h
\brief Newton's method.
*/
#ifndef __H2D_SOLVER_NEWTON_CONVERGENCE_MEASUREMENT_H_
#define __H2D_SOLVER_NEWTON_CONVERGENCE_MEASUREMENT_H_

#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class NewtonSolver;

    /// Convergence measurement strategies.
    /// This specifies the quantity that is compared to newton_tolerance (settable by set_tolerance()).
    enum NewtonSolverConvergenceMeasurement
    {
      ResidualNormRelativeToInitial,
      ResidualNormRelativeToPrevious,
      ResidualNormRatioToInitial,
      ResidualNormRatioToPrevious,
      ResidualNormAbsolute,
      SolutionDistanceFromPreviousAbsolute,
      SolutionDistanceFromPreviousRelative
    };

    template<typename Scalar>
    /// Convergence measurement function - returns converged true/false.
    HERMES_API bool newtonConverged(NewtonSolver<Scalar>*);
  }
}
#endif