// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_REFLINEARPROBLEM_H
#define __H2D_REFLINEARPROBLEM_H

//#include "linear_problem.h"
#include "ref_discrete_problem.h"
#include "views/order_view.h"

class Mesh;
class ExactSolution;

///
///
///
class H2D_API RefLinearProblem : public RefDiscreteProblem
{
public:
  RefLinearProblem(LinearProblem* base, int order_increase = 1, int refinement = 1);
  virtual ~RefLinearProblem();

protected:

  LinearProblem* base;
  Mesh** meshes;
  int order_increase;
  int refinement;
};

#endif
