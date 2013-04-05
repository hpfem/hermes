/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_DISCRETE_PROBLEM_LINEAR_H
#define __H2D_DISCRETE_PROBLEM_LINEAR_H

#include "discrete_problem.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// SIMPLIFICATION of Discrete problem class for LINEAR PROBLEMS.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    /// @ingroup inner
    class HERMES_API DiscreteProblemLinear : public DiscreteProblem<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblemLinear(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces);

      /// Constructor for one equation.
      DiscreteProblemLinear(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space);

      /// Empty constructor for special purposes.
      DiscreteProblemLinear();
    };
  }
}
#endif
