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

#include "config.h"
#include <hermes3d.h>

#define ERR_SUCCESS					0
#define ERR_FAILURE					-1


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

int main(int argc, char *argv[]) {
	set_verbose(false);

	if (argc < 1) error("Not enough parameters");

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("Loading mesh file '%s'\n", argv[1]);

	// apply refinements
	bool ok = true;
	for (int k = 2; k < argc && ok; k += 2) {
		int elem_id, reft_id;
		sscanf(argv[k], "%d", &elem_id);
		reft_id = parse_reft(argv[k + 1]);
		ok = mesh.refine_element(elem_id, reft_id);
	}

	mesh.regularize();

	mesh.dump();

#ifdef OUTPUT_DIR
	const char *of_name = OUTPUT_DIR "/mesh.gmsh";
	FILE *ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		GmshOutputEngine output(ofile);
		output.out(&mesh);

		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
#endif

	return ERR_SUCCESS;
}
