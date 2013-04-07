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
/*! \file hermes_common.h
    \brief File containing includes of all HermesCommon functionality + solvers. Intended to be included.
*/
#include "common.h"
#include "solvers/linear_matrix_solver.h"
#include "solvers/amesos_solver.h"
#include "solvers/aztecoo_solver.h"
#include "solvers/epetra.h"
#include "solvers/mumps_solver.h"
#include "solvers/petsc_solver.h"
#include "solvers/umfpack_solver.h"
#include "solvers/superlu_solver.h"
#include "solvers/precond.h"
#include "solvers/precond_ifpack.h"
#include "solvers/precond_ml.h"
#include "solvers/eigensolver.h"
#include "hermes_function.h"
#include "compat.h"
#include "callstack.h"
#include "vector.h"
#include "tables.h"
#include "array.h"
#include "qsort.h"
#include "ord.h"
#include "mixins.h"
#include "api.h"