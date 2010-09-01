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
// hdf5reader.cc
//
// Loading mesh from HDF5 format
//

#include "../h3dconfig.h"
#ifdef WITH_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#include "hdf5.h"
#include "../mesh.h"
#include <string.h>
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

// exception error codes
#define E_CANT_OPEN_FILE					-1
#define E_READ_ERROR						-2
#define E_INVALID_VERSION					-3
#define E_NOT_HDF5_FILE						-4
#define E_ERROR								-5
#define E_WRITE_ERROR						-6

HDF5Reader::HDF5Reader() {
	_F_
#ifdef WITH_HDF5
#else
	error("hermes3d was not built with HDF5 support.");
#endif
}

HDF5Reader::~HDF5Reader() {
	_F_
#ifdef WITH_HDF5
#endif
}


// Load ///////////////////////////////////////////////////////////////////////

#ifdef WITH_HDF5

// checks the version of hdf5 file
static bool check_version(hid_t id) {
	_F_
	herr_t status;

	hid_t version_attr = H5Aopen_name(id, "version");
	if (version_attr < 0) return false;

	char attr_data[2] = { 0 };
	status = H5Aread(version_attr, H5T_NATIVE_CHAR, attr_data);

	H5Aclose(version_attr);

	// check version
	int version = attr_data[0] * 0x100 + attr_data[1];
	if (version != 0x0100) return false;

	return true;
}

/// reads the count attribute in the group 'id'
///
/// @param[in] id identifier of the group
/// @param[in] name name of the attribute
/// @param[out] count the value of the attribute
static bool read_attr(hid_t id, const char *name, unsigned int &count) {
	_F_
	herr_t status;

	hid_t attr_id = H5Aopen_name(id, name);
	if (attr_id < 0) return false;

	status = H5Aread(attr_id, H5T_NATIVE_UINT32, &count);

	H5Aclose(attr_id);

	return (status >= 0);
}

static bool read_vertices(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	// open vertices group
	hid_t group_id = H5Gopen(id, "vertices");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				double pt[3];
				if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt) >= 0) {
					mesh->add_vertex(pt[0], pt[1], pt[2]);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id);

	return ret;
}

static bool read_hexes(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	// open group with hexes
	hid_t group_id = H5Gopen(id, "hex");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				Word_t vtcs[Hex::NUM_VERTICES] = { 0 };
				if (H5Dread(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs) >= 0) {
					mesh->add_hex(vtcs);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id); // close the group

	return ret;
}

static bool read_tetras(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	hid_t group_id = H5Gopen(id, "tetra");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				Word_t vtcs[Tetra::NUM_VERTICES] = { 0 };
				if (H5Dread(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs) >= 0) {
					mesh->add_tetra(vtcs);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id); // close the group

	return ret;
}

static bool read_prisms(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	hid_t group_id = H5Gopen(id, "prism");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				Word_t vtcs[Prism::NUM_VERTICES] = { 0 };
				if (H5Dread(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs) >= 0) {
					mesh->add_prism(vtcs);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id); // close the group

	return ret;
}

static bool read_elements(hid_t id, Mesh *mesh) {
	_F_
	hid_t group_id = H5Gopen(id, "elements");
	if (group_id < 0) return false;

	bool ret =
		read_hexes(group_id, mesh) &&
		read_tetras(group_id, mesh) &&
		read_prisms(group_id, mesh);

	H5Gclose(group_id); // close the group

	return ret;
}

// BCs ////////////////////////////////////////////////////////////////////////

static bool read_tris(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	hid_t group_id = H5Gopen(id, "tri");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				Word_t vtcs[Tri::NUM_VERTICES] = { 0 };
				unsigned int marker = 0;
				if (H5Dread(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs) >= 0 && read_attr(dataset_id, "marker", marker)) {
					mesh->add_tri_boundary(vtcs, marker);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id); // close the group

	return ret;
}

static bool read_quads(hid_t id, Mesh *mesh) {
	_F_
	bool ret = true;

	hid_t group_id = H5Gopen(id, "quad");
	if (group_id < 0) return false;

	// read the number of vertices
	unsigned int count;
	if (read_attr(group_id, "count", count)) {
		for (unsigned int i = 0; i < count; i++) {
			// open data set
			char name[16] = { 0 };
			sprintf(name, "%d", i);
			hid_t dataset_id = H5Dopen(group_id, name);
			if (dataset_id >= 0) {
				Word_t vtcs[Quad::NUM_VERTICES] = { 0 };
				unsigned int marker = 0;
				if (H5Dread(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs) >= 0 && read_attr(dataset_id, "marker", marker)) {
					mesh->add_quad_boundary(vtcs, marker);
				}
				else {
					H5Dclose(dataset_id);
					ret = false;
					break;
				}
			}
			H5Dclose(dataset_id);
		}
	}
	else ret = false;

	H5Gclose(group_id); // close the group

	return ret;
}

static bool read_bcs(hid_t id, Mesh *mesh) {
	_F_
	hid_t group_id = H5Gopen(id, "bc");
	if (group_id < 0) return false;

	bool ret =
		read_tris(group_id, mesh) &&
		read_quads(group_id, mesh);

	H5Gclose(group_id); // close the group

	return ret;
}

#endif

bool HDF5Reader::load(const char *file_name, Mesh *mesh) {
	_F_
#ifdef WITH_HDF5
	bool ret = true;

	H5open();
	try {
		// check if the file is HDF5
		int err = H5Fis_hdf5(file_name);
		if (err == 0) throw E_NOT_HDF5_FILE;
		else if (err < 0) throw E_ERROR;

		hid_t file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);

		hid_t mesh_group_id = H5Gopen(file_id, "/mesh3d");
		if (mesh_group_id < 0) throw E_READ_ERROR;

		// check version
		if (!check_version(mesh_group_id)) throw E_INVALID_VERSION;
		if (!read_vertices(mesh_group_id, mesh)) throw E_READ_ERROR;
		if (!read_elements(mesh_group_id, mesh)) throw E_READ_ERROR;
		if (!read_bcs(mesh_group_id, mesh)) throw E_READ_ERROR;

		H5Gclose(mesh_group_id);

		H5Fclose(file_id);

		mesh->ugh();
	}
	catch (int e) {
		// TODO: save the error code
		// TODO: close the mesh_group_id
		ret = false;
	}

	H5close();
	return ret;
#else
	return false;
#endif
}

// Save ///////////////////////////////////////////////////////////////////////

#ifdef WITH_HDF5

// not included in the class just to hide it
// TODO: improve error handling

static bool write_attr(hid_t loc_id, const char *name, uint value) {
	_F_
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr = H5Acreate(loc_id, name, H5T_NATIVE_UINT32, dataspace_id, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_UINT32, &value);
	H5Aclose(attr);
	H5Sclose(dataspace_id);
	return true;
}

static bool save_vertices(hid_t parent_group_id, Mesh *mesh) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "vertices", 0);

	// count
	uint count = mesh->vertices.count();
	write_attr(group_id, "count", count);

	hsize_t dims = 3;
	hid_t vertex_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		// Create the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_DOUBLE, vertex_dataspace_id, H5P_DEFAULT);
		double pt[] = { mesh->vertices[i]->x, mesh->vertices[i]->y, mesh->vertices[i]->z };
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pt);
		status = H5Dclose(dataset_id);
	}

	H5Sclose(vertex_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

// Elements ////

static bool save_hex(hid_t parent_group_id, Array<Element *> &elems) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "hex", 0);

	// count
	uint count = elems.count();
	write_attr(group_id, "count", count);

	///
	hsize_t dims = Hex::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		uint vtcs[Hex::NUM_VERTICES] = {
			elems[i]->get_vertex(0), elems[i]->get_vertex(1), elems[i]->get_vertex(2), elems[i]->get_vertex(3),
			elems[i]->get_vertex(4), elems[i]->get_vertex(5), elems[i]->get_vertex(6), elems[i]->get_vertex(7)
		};
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs);
		status = H5Dclose(dataset_id);
	}

	H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

static bool save_tetra(hid_t parent_group_id, Array<Element *> &elems) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "tetra", 0);

	// count
	uint count = elems.count();
	write_attr(group_id, "count", count);

	///
	hsize_t dims = Tetra::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		uint vtcs[Tetra::NUM_VERTICES] = {
			elems[i]->get_vertex(0), elems[i]->get_vertex(1), elems[i]->get_vertex(2), elems[i]->get_vertex(3)
		};
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs);
		status = H5Dclose(dataset_id);
	}

	H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

static bool save_prism(hid_t parent_group_id, Array<Element *> &elems) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "prism", 0);

	// count
	uint count = elems.count();
	write_attr(group_id, "count", count);

	///
	hsize_t dims = Prism::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		uint vtcs[Prism::NUM_VERTICES] = {
			elems[i]->get_vertex(0), elems[i]->get_vertex(1), elems[i]->get_vertex(2),
			elems[i]->get_vertex(3), elems[i]->get_vertex(4), elems[i]->get_vertex(5)
		};
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs);
		status = H5Dclose(dataset_id);
	}

	H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

static bool save_elements(hid_t parent_group_id, Mesh *mesh) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "elements", 0);

	// count
	uint count = mesh->elements.count();
	write_attr(group_id, "count", count);

	///
	Array<Element *> tet, hex, pri;
	for (uint i = 0; i < count; i++) {
		Element *elem = mesh->elements[i];
		switch (elem->get_mode()) {
			case MODE_TETRAHEDRON: tet.add(elem); break;
			case MODE_HEXAHEDRON: hex.add(elem); break;
			case MODE_PRISM: pri.add(elem); break;
		}
	}

	bool ret =
		save_tetra(group_id, tet) &&
		save_hex(group_id, hex) &&
		save_prism(group_id, pri);

	status = H5Gclose(group_id); // close the group

	return ret;
}

// BC ////

static bool save_tri_bc(hid_t parent_group_id, Mesh *mesh, Array<Word_t> &bcs) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "tri", 0);

	// count
	uint count = bcs.count();
	write_attr(group_id, "count", count);

	///
	hsize_t dims = Tri::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		Word_t fid = bcs[i];
		Facet *facet = mesh->facets.get(fid);
		Element *elem = mesh->elements[facet->left];
		Boundary *bnd = mesh->boundaries[facet->right];

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		const int *vidx = elem->get_face_vertices(facet->left_face_num);
		uint vtcs[Tri::NUM_VERTICES] = { elem->get_vertex(vidx[0]), elem->get_vertex(vidx[1]), elem->get_vertex(vidx[2]) };
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs);

		// marker
		write_attr(dataset_id, "marker", bnd->marker);

		status = H5Dclose(dataset_id);
	}

	H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

static bool save_quad_bc(hid_t parent_group_id, Mesh *mesh, Array<Word_t> &bcs) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "quad", 0);

	// count
	uint count = bcs.count();
	write_attr(group_id, "count", count);

	///
	hsize_t dims = Quad::NUM_VERTICES;
	hid_t elem_dataspace_id = H5Screate_simple(1, &dims, NULL);

	// dump vertices
	for (uint i = 0; i < count; i++) {
		char name[256];
		sprintf(name, "%d", i);

		Word_t fid = bcs[i];
		Facet *facet = mesh->facets.get(fid);
		Element *elem = mesh->elements[facet->left];
		Boundary *bnd = mesh->boundaries[facet->right];

		// the dataset
		hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_UINT32, elem_dataspace_id, H5P_DEFAULT);
		const int *vidx = elem->get_face_vertices(facet->left_face_num);
		uint vtcs[Quad::NUM_VERTICES] = { elem->get_vertex(vidx[0]), elem->get_vertex(vidx[1]), elem->get_vertex(vidx[2]), elem->get_vertex(vidx[3]) };
		elem->get_face_vertices(facet->left_face_num);
		status = H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vtcs);

		// marker
		write_attr(dataset_id, "marker", bnd->marker);

		status = H5Dclose(dataset_id);
	}

	H5Sclose(elem_dataspace_id);

	status = H5Gclose(group_id); // close the group

	return true;
}

static bool save_bc(hid_t parent_group_id, Mesh *mesh) {
	_F_
	herr_t status;

	// create main group
	hid_t group_id = H5Gcreate(parent_group_id, "bc", 0);

	// count
	uint count = mesh->boundaries.count();
	write_attr(group_id, "count", count);

	// sort out boundaries that are triangular and quadrilateral
	Array<Word_t> tri, quad;
	FOR_ALL_FACETS(fid, mesh) {
		Facet *facet = mesh->facets.get(fid);
		if (facet->type == Facet::OUTER) {
			switch (facet->mode) {
				case MODE_TRIANGLE: tri.add(fid); break;
				case MODE_QUAD: quad.add(fid); break;
			}
		}
	}

	bool ret =
		save_tri_bc(group_id, mesh, tri) &&
		save_quad_bc(group_id, mesh, quad);

	status = H5Gclose(group_id); // close the group

	return ret;
}

#endif

bool HDF5Reader::save(const char *file_name, Mesh *mesh) {
	_F_
#ifdef WITH_HDF5
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
	if (description != NULL) {
		hid_t type = H5Tcopy(H5T_C_S1);
		status = H5Tset_size(type, H5T_VARIABLE);
		hid_t dataspace_id2 = H5Screate(H5S_SCALAR);
		hid_t attr_descr = H5Acreate(mesh_group_id, "description", type, dataspace_id2, H5P_DEFAULT);
		status = H5Awrite(attr_descr, type, &description);
		H5Aclose(attr_descr);
		H5Tclose(type);
		H5Sclose(dataspace_id2);
	}

	// mesh
	bool ret =
		save_vertices(mesh_group_id, mesh) &&
		save_elements(mesh_group_id, mesh) &&
		save_bc(mesh_group_id, mesh);

	status = H5Gclose(mesh_group_id); // close the group
	status = H5Fclose(file_id); // close the file

	// deinit HDF5
	H5close();

	return ret;
#else
	return false;
#endif
}
