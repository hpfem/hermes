// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

/*
 * main.cc
 *
 * Test for H1 lobatto shapeset for Hex
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// forward declarations (simpler than defining a special header file for each module that exports just one function)
bool test_lin_indep(Shapeset *shapeset);
bool test_zero_values(Shapeset *shapeset);
bool test_continuity(Shapeset *shapeset);
bool test_gradients(Shapeset *shapeset);
bool test_gradients_directly(Shapeset *shapeset);

//
// main
//
int main(int argc, char *argv[]) {
	_F_
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	H1ShapesetLobattoHex shapeset;

	try {
		// I. linear independency
		if (!test_lin_indep(&shapeset)) throw ERR_FAILURE;
		// II. test zero fn. values
		if (!test_zero_values(&shapeset)) throw ERR_FAILURE;
		// III. continuity on boundaries
		if (!test_continuity(&shapeset)) throw ERR_FAILURE;
		// IV. gradients
		if (!test_gradients(&shapeset)) throw ERR_FAILURE;
		// V. computes gradients numericaly from fn values and compares
		if (!test_gradients_directly(&shapeset)) throw ERR_FAILURE;

		printf("Shapeset OK\n");
	}
	catch (int e) {
		printf("Test failed\n");
		res = e;
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}
