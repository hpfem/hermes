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

/*
 * main.cc
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include "mesh.h"
#include "mesh3dreader.h"

extern "C" {
#include <hdf5.h>
}

bool save_vertices(hid_t parent_group_id, Mesh *mesh) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "vertices", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = mesh->vertices.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

	hsize_t dims = 3;
	hid_t vertex_dataspace_id = H5Screate_simple(1, &dims, NULL);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// Create the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, vertex_dataspace_id, H5P_DEFAULT);
		double pt[] = { mesh->vertices[i]->x, mesh->vertices[i]->y, mesh->vertices[i]->z };
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt);
		status = H5Dclose(dataset_id);
    }

    H5Sclose(vertex_dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

// Elements ////

bool save_hex(hid_t parent_group_id, Array<Element *> &elems) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "hex", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = elems.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

    ///
	hsize_t dims = Hex::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, elems[i]->get_vertices());
		status = H5Dclose(dataset_id);
    }

    H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

bool save_tetra(hid_t parent_group_id, Array<Element *> &elems) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "tetra", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = elems.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

    ///
	hsize_t dims = Tetra::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, elems[i]->get_vertices());
		status = H5Dclose(dataset_id);
    }

    H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

bool save_prism(hid_t parent_group_id, Array<Element *> &elems) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "prism", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = elems.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

    ///
	hsize_t dims = Prism::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, elems[i]->get_vertices());
		status = H5Dclose(dataset_id);
    }

    H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

bool save_elements(hid_t parent_group_id, Mesh *mesh) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "elements", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = mesh->elements.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

    ///
	Array<Element *> tet, hex, pri;
    for (int i = 0; i < count; i++) {
    	Element *elem = mesh->elements[i];
    	switch (elem->get_mode()) {
    		case MODE_TETRAHEDRON: tet.add(elem); break;
    		case MODE_HEXAHEDRON: hex.add(elem); break;
    		case MODE_PRISM: pri.add(elem); break;
    	}
    }

    save_tetra(group_id, tet);
    save_hex(group_id, hex);
    save_prism(group_id, pri);

	status = H5Gclose(group_id);		// close the group
}

// BC ////

bool save_tri_bc(hid_t parent_group_id, Array<Boundary *> &bcs) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "tri", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = bcs.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);

    ///
	hsize_t dims = Tri::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, bcs[i]->get_vertices());

		// marker
		hid_t attr_marker = H5Acreate(dataset_id, "marker", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
		uint marker = bcs[i]->get_marker();
		status = H5Awrite(attr_marker, H5T_NATIVE_UINT32, &marker);
		H5Aclose(attr_marker);

		status = H5Dclose(dataset_id);
    }

    H5Sclose(elem_dataspace_id);
    H5Sclose(dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

bool save_quad_bc(hid_t parent_group_id, Array<Boundary *> &bcs) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "quad", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = bcs.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);

    ///
	hsize_t dims = Quad::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	hid_t merker_dataspace_id = H5Screate(H5S_SCALAR);

    // dump vertices
    for (int i = 0; i < count; i++) {
    	char name[256];
    	sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, bcs[i]->get_vertices());

		// marker
		hid_t attr_marker = H5Acreate(dataset_id, "marker", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
		uint marker = bcs[i]->get_marker();
		status = H5Awrite(attr_marker, H5T_NATIVE_UINT32, &marker);
		H5Aclose(attr_marker);

		status = H5Dclose(dataset_id);
    }

    H5Sclose(elem_dataspace_id);
    H5Sclose(dataspace_id);

	status = H5Gclose(group_id);		// close the group
}

bool save_bc(hid_t parent_group_id, Mesh *mesh) {
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "bc", 0);

	// count
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_count = H5Acreate(group_id, "count", H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	uint count = mesh->boundaries.count();
	status = H5Awrite(attr_count, H5T_NATIVE_UINT32, &count);
	H5Aclose(attr_count);
    H5Sclose(dataspace_id);

    ///
    ///
	Array<Boundary *> tri, quad;
    for (int i = 0; i < count; i++) {
    	Boundary *bnd = mesh->boundaries[i];
    	switch (bnd->get_mode()) {
    		case MODE_TRIANGLE: tri.add(bnd); break;
    		case MODE_QUAD: quad.add(bnd); break;
    	}
    }

    save_tri_bc(group_id, tri);
    save_quad_bc(group_id, quad);

	status = H5Gclose(group_id);		// close the group
}

bool dump_hdf5(const char *file_name, Mesh *mesh) {
	herr_t status;

	// init HDF5
	H5open();

	// create a file
	hid_t file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// create main group
	hid_t mesh_group_id = H5Gcreate(file_id, "/mesh3d", 0);

	// version
	hsize_t dims = 2;
	hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
	hid_t attr_ver = H5Acreate(mesh_group_id, "version", H5T_STD_I8BE, dataspace_id, H5P_DEFAULT);
	char attr_data[2] = { 1, 0 };
	status = H5Awrite(attr_ver, H5T_NATIVE_CHAR, attr_data);
	H5Aclose(attr_ver);
    H5Sclose(dataspace_id);

	// description
    hid_t type = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(type, H5T_VARIABLE);
    const char *descr = "Test mesh";
	hid_t dataspace_id2 = H5Screate(H5S_SCALAR);
	hid_t attr_descr = H5Acreate(mesh_group_id, "description", type, dataspace_id2, H5P_DEFAULT);
	status = H5Awrite(attr_descr, type, &descr);
	H5Aclose(attr_descr);
    H5Tclose(type);
    H5Sclose(dataspace_id2);


    save_vertices(mesh_group_id, mesh);
    save_elements(mesh_group_id, mesh);
    save_bc(mesh_group_id, mesh);

	status = H5Gclose(mesh_group_id);		// close the group
	status = H5Fclose(file_id);			// close the file

	// deinit HDF5
	H5close();

	return 0;
}

//

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Not enough parameters.\n");
		return -1;
	}

	Mesh3DLoader loader;
	Mesh mesh;
	if (mesh.load(argv[1], &loader)) {
		dump_hdf5(argv[2], &mesh);


		printf("done\n");
	}
	else {
		printf("Error reading mesh file '%s'.\n", argv[1]);
	}


	return 0;
}

