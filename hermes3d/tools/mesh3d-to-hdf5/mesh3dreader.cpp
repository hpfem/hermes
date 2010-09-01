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

//
// mesh3dreader.cc
//
// Loading mesh from mesh3d format
//

//#include "config.h"

//#include "common.h"
#include <string.h>
//#include <common/trace.h>
//#include <common/error.h>
#include "mesh3dreader.h"
#include "mesh.h"

#define ERROR(...)

// maximal row length in bytes (used for reading the mesh3d-file)
#define MAX_ROW_LEN							1024

// exception error codes
#define E_CANT_OPEN_FILE					-1
#define E_READ_ERROR						-2

// number of markers on mesd3d file
#define MARKERS								1

Mesh3DLoader::Mesh3DLoader() {
}

Mesh3DLoader::~Mesh3DLoader() {
}

static bool read_num(char *row, int &vertex_count) {
	int n;
	if (sscanf(row, "%d", &n) == 1) {
		vertex_count = n;
		return true;
	} else
		return false;
}

static bool read_n_nums(char *row, int n, double values[]) {
	int i = 0;
	char delims[] = " \t\n\r";
	char *token = strtok(row, delims);
	while (token != NULL && i < n) {
		double d;
		sscanf(token, "%lf", &d);
		values[i++] = d;

		token = strtok(NULL, delims);
	}

	return (i == n);
}

static bool read_n_nums(char *row, int n, uint values[]) {
	int i = 0;
	char delims[] = " \t\n\r";
	char *token = strtok(row, delims);
	while (token != NULL && i < n) {
		int n;
		sscanf(token, "%d", &n);
		values[i++] = n;

		token = strtok(NULL, delims);
	}

	return (i == n);
}

static bool range_check(int max_index, uint *vs, int num_vs) {
	for (int i = 0; i < num_vs; i++) {
		if (vs[i] <= 0 || vs[i] > max_index)
			return false;
	}
	return true;
}

bool Mesh3DLoader::load(const char *file_name, Mesh *mesh) {
	FILE *file = fopen(file_name, "r");
	if (file == NULL)
		return E_CANT_OPEN_FILE;

	try {
		enum EState {
			STATE_VERTICES_NUM,
			STATE_VERTICES,
			STATE_TETRAS_NUM,
			STATE_TETRAS,
			STATE_HEXES_NUM,
			STATE_HEXES,
			STATE_PRISMS_NUM,
			STATE_PRISMS,
			STATE_TRIS_NUM,
			STATE_TRIS,
			STATE_QUADS_NUM,
			STATE_QUADS,
			STATE_OK
		} state = STATE_VERTICES_NUM;

		int vertex_count = 0;
		int tetra_count = 0;
		int hex_count = 0;
		int prism_count = 0;
		int tri_count = 0;
		int quad_count = 0;
		double buffer[10];
		uint vs[10];

		int max_vertex_index = 0;

		char row[MAX_ROW_LEN] = {0};
		while (fgets(row, MAX_ROW_LEN, file) != NULL && state != STATE_OK) {
			if (row[0] == '#') continue; // comment

			// trim right
			int i = strlen(row) - 1;
			for (; i >= 0; i--) {
				if (row[i] != ' ' && row[i] != '\r' && row[i] != '\n' && row[i] != '\t')
					break;
			}
			row[i + 1] = '\0';
			if (strlen(row) <= 0) continue; // skip empty lines

			switch (state) {
				case STATE_VERTICES_NUM:
					if (read_num(row, vertex_count)) {
						state = STATE_VERTICES;
						if (vertex_count <= 0)
							throw E_READ_ERROR;
						max_vertex_index = vertex_count; //vertices are counted from 1 in mesh3d format
					}
					break;

				case STATE_VERTICES:
					if (read_n_nums(row, Vertex::NUM_COORDS, buffer)) {
						mesh->add_vertex(buffer[0], buffer[1], buffer[2]);

						vertex_count--;
						if (vertex_count == 0)
							state = STATE_TETRAS_NUM;
					}
					break;

				case STATE_TETRAS_NUM:
					if (read_num(row, tetra_count)) {
						state = STATE_TETRAS;
						if (tetra_count <= 0)
							state = STATE_HEXES_NUM;
					}
					break;

				case STATE_TETRAS:
				if (read_n_nums(row, Tetra::NUM_VERTICES, vs)) {
					if (!range_check(max_vertex_index, vs, Tetra::NUM_VERTICES)) {
						ERROR("Invalid vertex index was found in the section defining tetras.");
						throw E_READ_ERROR;
					}

					for (int i = 0; i < Tetra::NUM_VERTICES; i++) vs[i]--;
					Tetra *tetra = mesh->add_tetra(vs);

					tetra_count--;
					if (tetra_count == 0)
						state = STATE_HEXES_NUM;
				}
				break;

				case STATE_HEXES_NUM:
					if (read_num(row, hex_count)) {
						state = STATE_HEXES;
						if (hex_count <= 0)
							state = STATE_PRISMS_NUM;
					}
					break;

				case STATE_HEXES:
					if (read_n_nums(row, Hex::NUM_VERTICES, vs)) {
						if (!range_check(max_vertex_index, vs, Hex::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining hexes.");
							throw E_READ_ERROR;
						}

						for (int i = 0; i < Hex::NUM_VERTICES; i++) vs[i]--;
						Hex *hex = mesh->add_hex(vs);

						hex_count--;
						if (hex_count == 0)
							state = STATE_PRISMS_NUM;
					}
					break;

				case STATE_PRISMS_NUM:
					if (read_num(row, prism_count)) {
						state = STATE_PRISMS;
						if (prism_count <= 0)
							state = STATE_TRIS_NUM;
					}
					break;

				case STATE_PRISMS:
					if (read_n_nums(row, Prism::NUM_VERTICES, vs)) {
						if (!range_check(max_vertex_index, vs, Prism::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining prisms.");
							throw E_READ_ERROR;
						}
						// add prism
						for (int i = 0; i < Prism::NUM_VERTICES; i++) vs[i]--;
						Prism *prism = mesh->add_prism(vs);
						prism_count--;
						if (prism_count == 0)
							state = STATE_TRIS_NUM;
					}
					break;

				case STATE_TRIS_NUM:
					if (read_num(row, tri_count)) {
						state = STATE_TRIS;
						if (tri_count <= 0)
							state = STATE_QUADS_NUM;
					}
					break;

				case STATE_TRIS:
					if (read_n_nums(row, Tri::NUM_VERTICES + MARKERS, vs)) {
						if (!range_check(max_vertex_index, vs, Tri::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining tris.");
							throw E_READ_ERROR;
						}

						uint facet_idxs[Tri::NUM_VERTICES] = {vs[0] - 1, vs[1] - 1, vs[2] - 1};
						mesh->add_tri_boundary(facet_idxs, (int) vs[Tri::NUM_VERTICES]);
						tri_count--;
						if (tri_count == 0)
							state = STATE_QUADS_NUM;
					}
					else {
						ERROR("Not enough information for tris. You probably forgot to define boundary condition.");
						throw E_READ_ERROR;
					}
					break;

				case STATE_QUADS_NUM:
					if (read_num(row, quad_count)) {
						state = STATE_QUADS;
						if (quad_count <= 0)
							state = STATE_QUADS;
					}
					break;

				case STATE_QUADS:
					if (read_n_nums(row, Quad::NUM_VERTICES + MARKERS, vs)) {
						if (!range_check(max_vertex_index, vs, Quad::NUM_VERTICES)) {
							ERROR("Invalid vertex index was found in the section defining quads.");
							throw E_READ_ERROR;
						}

						uint facet_idxs[Quad::NUM_VERTICES] = { vs[0] - 1, vs[1] - 1, vs[2] - 1, vs[3] - 1 };
						mesh->add_quad_boundary(facet_idxs, (int) vs[Quad::NUM_VERTICES]);
						quad_count--;
						if (quad_count == 0)
							state = STATE_OK;
					}
					else {
						ERROR("Not enough information for quads. You probably forgot to define boundary condition.");
						throw E_READ_ERROR;
					}
					break;

				case STATE_OK:
				break;
			}
		}

		// check if all "outer" faces have defined boundary condition
		for (int i = mesh->facets.first(); i != INVALID_IDX; i = mesh->facets.next(i)) {
			Facet *facet = mesh->facets.get(i);

			if ((facet->left == INVALID_IDX) || (facet->right == INVALID_IDX)) {
				ERROR("Not all outer faces have defined boundary condition.");
				throw E_READ_ERROR;
			}
		}
	}
	catch (int e) {
		fclose(file);
		return false;
	}

	fclose(file);
	return true;
}
