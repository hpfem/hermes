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
 * mesh-copy.cc
 *
 * usage: $0 <mesh file>
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

enum ECopyType {
	CT_NONE,
	CT_COPY,
	CT_COPY_BASE
};

// helpers ////////////////////////////////////////////////////////////////////////////////////////

int parse_reft(char *str) {
	if (strcasecmp(str, "x") == 0) return H3D_REFT_HEX_X;
	else if (strcasecmp(str, "y") == 0) return H3D_REFT_HEX_Y;
	else if (strcasecmp(str, "z") == 0) return H3D_REFT_HEX_Z;
	else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return H3D_H3D_REFT_HEX_XY;
	else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return H3D_H3D_REFT_HEX_XZ;
	else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return H3D_H3D_REFT_HEX_YZ;
	else if (strcasecmp(str, "xyz") == 0) return H3D_H3D_H3D_REFT_HEX_XYZ;
	else return H3D_REFT_HEX_NONE;
}

ECopyType parse_copy_type(char *str) {
	if (strcasecmp(str, "copy") == 0) return CT_COPY;
	else if (strcasecmp(str, "copy_base") == 0) return CT_COPY_BASE;
	else return CT_NONE;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 2) error("Not enough parameters");

	ECopyType copy_type = parse_copy_type(args[1]);
	if (copy_type == CT_NONE) error("Unknown copy type (%s).\n", args[1]);

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[2], &mesh)) error("Loading mesh file '%s'\n", args[2]);

	// apply refinements
	bool ok = true;
	for (int k = 3; k < argc && ok; k += 2) {
		int elem_id, reft_id;
		sscanf(args[k], "%d", &elem_id);
		reft_id = parse_reft(args[k + 1]);
		ok = mesh.refine_element(elem_id, reft_id);
	}

	if (ok) {
		Mesh dup;

		switch (copy_type) {
			case CT_COPY:      dup.copy(mesh); break;
			case CT_COPY_BASE: dup.copy_base(mesh); break;
			default: break;
		}
		dup.dump();
	}
	else {
		warning("Unable to refine a mesh.");
		res = ERR_FAILURE;
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}
