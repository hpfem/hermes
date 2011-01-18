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

#include "exodusii.h"
#include <string.h>
#include "../../../hermes_common/error.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/callstack.h"
#include "../mesh.h"
#include "../refdomain.h"

#ifdef WITH_EXODUSII
#include <exodusII.h>
#endif

// TODO: error checking

ExodusIIReader::ExodusIIReader()
{
	_F_
#ifdef WITH_EXODUSII
#else
	error("hermes3d was not compiled with ExodusII support");
#endif
}

ExodusIIReader::~ExodusIIReader()
{
	_F_
}

bool ExodusIIReader::load(const char *file_name, Mesh *mesh)
{
	_F_
#ifdef WITH_EXODUSII
	int err;
	int cpu_ws = sizeof(double);		// use float or double
	int io_ws = 8;						// store variables as doubles
	float version;
	int exoid = ex_open(file_name, EX_READ, &cpu_ws, &io_ws, &version);

	// read initialization parameters
	int n_dims, n_nodes, n_elems, n_eblocks, n_nodesets, n_sidesets;
	char title[MAX_LINE_LENGTH + 1];
	err = ex_get_init(exoid, title, &n_dims, &n_nodes, &n_elems, &n_eblocks, &n_nodesets, &n_sidesets);
	if (n_dims != 3) {
		error("File '%s' does not contain 3D mesh.", file_name);
	}

	// load coordinates
	double *x = new double [n_nodes]; MEM_CHECK(x);
	double *y = new double [n_nodes]; MEM_CHECK(y);
	double *z = new double [n_nodes]; MEM_CHECK(z);
	err = ex_get_coord(exoid, x, y, z);
	for (int i = 0; i < n_nodes; i++)
		mesh->vertices.insert(std::pair<unsigned int, Vertex*> (i + 1, new Vertex(x[i], y[i], z[i])));
	delete [] x;
	delete [] y;
	delete [] z;

	// get info about element blocks
	int *eid_blocks = new int [n_eblocks]; MEM_CHECK(eid_blocks);
	err = ex_get_elem_blk_ids(exoid, eid_blocks);
	// go over all element blocks
	for (int i = 0; i < n_eblocks; i++) {
		int id = eid_blocks[i];

		// get block info
		char elem_type[MAX_STR_LENGTH + 1];
		int n_elems_in_blk, n_elem_nodes, n_attrs;
		err = ex_get_elem_block(exoid, id, elem_type, &n_elems_in_blk, &n_elem_nodes, &n_attrs);

		// read connectivity array
		int *connect = new int [n_elem_nodes * n_elems_in_blk]; MEM_CHECK(connect);
		err = ex_get_elem_conn(exoid, id, connect);

		// load into mesh
		int ic = 0;
		for (int j = 0; j < n_elems_in_blk; j++) {
			// convert connectivity ints into unsigned int
			unsigned int *vtcs = new unsigned int[n_elem_nodes];
			for (int k = 0; k < n_elem_nodes; k++, ic++) vtcs[k] = connect[ic];

			if (n_elem_nodes == Tetra::NUM_VERTICES) {
				Tetra *tet = mesh->add_tetra(vtcs);
				tet->marker = id;
			}
			else if (n_elem_nodes == Hex::NUM_VERTICES) {
				Hex *hex = mesh->add_hex(vtcs);
				hex->marker = id;
			}
			else
				error("Unknown type of element.");
      delete [] vtcs;
		}
		delete [] connect;
	}
	delete [] eid_blocks;

	// query number of side sets
	int *sid_blocks = new int [n_sidesets]; MEM_CHECK(sid_blocks);
	err = ex_get_side_set_ids(exoid, sid_blocks);
	// go over the sidesets
	for (int i = 0; i < n_sidesets; i++) {
		int n_sides_in_set, n_df_in_set;

		int sid = sid_blocks[i];

		err = ex_get_side_set_param(exoid, sid, &n_sides_in_set, &n_df_in_set);
		int num_elem_in_set = n_sides_in_set;

		int *elem_list = new int [num_elem_in_set]; MEM_CHECK(elem_list);
		int *side_list = new int [n_sides_in_set]; MEM_CHECK(side_list);

		err = ex_get_side_set(exoid, sid, elem_list, side_list);

		const int face_num[] = { 2, 1, 3, 0, 4, 5 };   // mapping from exodus to hermes face numbers
		for (int j = 0; j < num_elem_in_set; j++) {
			Element *elem = mesh->elements[elem_list[j]];
			int iface = face_num[side_list[j] - 1];
			int nv = elem->get_num_face_vertices(iface);
			unsigned int *vtcs = new unsigned int[nv];
			elem->get_face_vertices(iface, vtcs);

			switch (elem->get_face_mode(iface)) {
				case HERMES_MODE_TRIANGLE: mesh->add_tri_boundary(vtcs, sid); break;
				case HERMES_MODE_QUAD: mesh->add_quad_boundary(vtcs, sid); break;
				default: warning("Unknown type of face"); break;
			}
      delete [] vtcs;
		}

		delete [] elem_list;
		delete [] side_list;
	}
	delete [] sid_blocks;

	// we are done
	err = ex_close(exoid);

	mesh->ugh();

	return true;
#else
	return false;
#endif
}

bool ExodusIIReader::save(const char *file_name, Mesh *mesh)
{
	error(HERMES_ERR_NOT_IMPLEMENTED);
  return false;
}
