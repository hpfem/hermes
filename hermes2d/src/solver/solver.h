// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "../common.h"  // Also includes preprocessor definitions for the various 
                        // solver libraries via config.h.

/// @defgroup solvers Solvers
///
/// TODO: description

class FeProblem;

/// Abstract class for defining solver interface
///
///
/// TODO: Adjust interface to support faster update of matrix and rhs
///
/// @ingroup solvers
class Solver {
public:
	Solver() { sln = NULL; time = -1.0; }
	virtual ~Solver() { if (sln != NULL) delete [] sln; }

	virtual bool solve() = 0;
	scalar *get_solution() { return sln; }

	int get_error() { return error; }
	double get_time() { return time; }
        

protected:
	scalar *sln;
	int error;
	double time;			/// time spent on solving (in secs)
};


/// Abstract class for defining interface for LinearSolvers
///
///
/// @ingroup solvers
class LinearSolver : public Solver 
{
public:
	LinearSolver() : Solver() {}
};

/// Abstract class for defining interface for LinearSolvers
///
///
/// @ingroup solvers
class NonlinearSolver : public Solver {
public:
	NonlinearSolver() : Solver() { fp = NULL; }
	NonlinearSolver(FeProblem *fp) : Solver() { this->fp = fp; }

protected:
	FeProblem *fp;        // FE problem being solved (not NULL in case of using
	                      // NonlinearProblem(DiscreteProblem *) ctor
};

class IterSolver : public Solver
{
public:
  IterSolver() { sln = NULL; }
  int get_error() { return error; }
};

#endif
