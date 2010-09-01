// This file is part of Hermes2D
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

#ifndef __H2D_ITERSOLVER_H
#define __H2D_ITERSOLVER_H

#include "common.h"

/// Abstract class for defining solver interface
///
///
/// TODO: Adjust interface to support faster update of matrix and rhs
///
/// @ingroup solvers
class IterSolver
{
public:
	IterSolver() { sln = NULL; }
	virtual ~IterSolver() { delete [] sln; }

	virtual bool solve() = 0;
	scalar *get_solution_vector() { return sln; }

	int get_error() { return error; }

protected:
	scalar *sln;
	int error;
};

#endif
