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
#include "hermes_common/common.h"
#include "hermes_common/solvers/linear_solver.h"
#include "hermes_common/solvers/nonlinear_solver.h"
#include "hermes_common/solvers/amesos_solver.h"
#include "hermes_common/solvers/aztecoo_solver.h"
#include "hermes_common/solvers/epetra.h"
#include "hermes_common/solvers/mumps_solver.h"
#include "hermes_common/solvers/nox_solver.h"
#include "hermes_common/solvers/petsc_solver.h"
#include "hermes_common/solvers/umfpack_solver.h"
#include "hermes_common/solvers/superlu_solver.h"
#include "hermes_common/solvers/precond.h"
#include "hermes_common/solvers/precond_ifpack.h"
#include "hermes_common/solvers/precond_ml.h"
#include "hermes_common/solvers/eigensolver.h"
