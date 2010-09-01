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

#include "common.h"
#include "space.h"
#include "weakform.h"
#include "linear_problem.h"
#include "ref_linear_problem.h"
#include "solution.h"

RefLinearProblem::RefLinearProblem(LinearProblem* base, int order_increase, 
				   int refinement) : RefDiscreteProblem(base, order_increase, refinement) {};

RefLinearProblem::~RefLinearProblem()
{
}



