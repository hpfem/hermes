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

#include "h3dconfig.h"
#include "common.h"
#include "mesh.h"
#include "refmap.h"
#include "refdomain.h"
#include <common/error.h>
#include <common/trace.h>
#include <common/callstack.h>

#include "shapeset/common.h"
#include "shapeset/refmapss.h"
#include "determinant.h"

//

#ifdef WITH_TETRA
	static RefMapShapesetTetra		ref_map_shapeset_tetra;
	static ShapeFunction 			ref_map_pss_tetra(&ref_map_shapeset_tetra);
	#define H3D_REFMAP_SHAPESET_TETRA	&ref_map_shapeset_tetra
	#define H3D_REFMAP_PSS_TETRA		&ref_map_pss_tetra
#else
	#define H3D_REFMAP_SHAPESET_TETRA	NULL
	#define H3D_REFMAP_PSS_TETRA		NULL
#endif

#ifdef WITH_HEX
	static RefMapShapesetHex 		ref_map_shapeset_hex;
	static ShapeFunction			ref_map_pss_hex(&ref_map_shapeset_hex);
	#define H3D_REFMAP_SHAPESET_HEX		&ref_map_shapeset_hex
	#define H3D_REFMAP_PSS_HEX			&ref_map_pss_hex
#else
	#define H3D_REFMAP_SHAPESET_HEX		NULL
	#define H3D_REFMAP_PSS_HEX			NULL
#endif

// TODO: prisms

static ShapeFunction *ref_map_pss[] = { H3D_REFMAP_PSS_TETRA, H3D_REFMAP_PSS_HEX, NULL };

// RefMap /////////////////////////////////////////////////////////////////////////////////////////

RefMap::RefMap() {
	_F_
	this->mesh = NULL;
	this->pss = NULL;
}

RefMap::RefMap(Mesh *mesh) {
	_F_
	this->mesh = mesh;
	this->pss = NULL;
}

RefMap::~RefMap() {
	_F_
}

void RefMap::set_active_element(Element *e) {
	_F_
	assert(e != NULL);

	EMode3D mode = e->get_mode();

	pss = ref_map_pss[mode];
	pss->set_active_element(e);

	if (e == element) return;
	element = e;

	reset_transform();

	is_const_jacobian = mode == MODE_TETRAHEDRON;

	int nvertices = element->get_num_vertices();

	// prepare the shapes and coefficients of the reference map
	Shapeset *shapeset = this->pss->get_shapeset();
	int i, k = 0;
	for (i = 0; i < nvertices; i++)
		indices[k++] = shapeset->get_vertex_index(i);

	// straight element
	for (int iv = 0; iv < nvertices; iv++)
		vertex[iv] = *mesh->vertices[e->get_vertex(iv)];
	coefs = vertex;
	n_coefs = nvertices;

	// calculate the order of the reference map
	switch (mode) {
		case MODE_TETRAHEDRON: ref_order = order3_t(0); break;
		case MODE_HEXAHEDRON:  ref_order = order3_t(1, 1, 1); break;
		case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
	}

	// calculate the order of the inverse reference map
	switch (mode) {
		case MODE_TETRAHEDRON: inv_ref_order = order3_t(0); break;
		case MODE_HEXAHEDRON:  inv_ref_order = order3_t(1, 1, 1); break;
		case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
	}

	// constant inverse reference map
	if (this->is_const_jacobian) calc_const_inv_ref_map();
	else const_jacobian = 0.0;
}


void RefMap::push_transform(int son) {
	_F_
	Transformable::push_transform(son);
	const_jacobian *= 0.125;
}

void RefMap::pop_transform() {
	_F_
	Transformable::pop_transform();
	const_jacobian *= 8;
}

void RefMap::force_transform(uint64 sub_idx, Trf *ctm) {
	_F_
	this->sub_idx = sub_idx;
	stack[top] = *ctm;
	ctm = stack + top;
	if (is_const_jacobian) calc_const_inv_ref_map();
}

double3x3 *RefMap::get_ref_map(const int np, const QuadPt3D *pt) {
	_F_

	double3x3 *m = new double3x3[np]; MEM_CHECK(m);
	memset(m, 0, np * sizeof(double3x3));

	if (is_const_jacobian) {
		for (int i = 0; i < np; i++)
			memcpy(m + i, const_ref_map, sizeof(double3x3));
	}
	else {
		pss->force_transform(sub_idx, ctm);
		for (int i = 0; i < n_coefs; i++) {
			double *dx, *dy, *dz;

			pss->set_active_shape(indices[i]);
			pss->precalculate(np, pt, FN_DEFAULT);
			pss->get_dx_dy_dz_values(dx, dy, dz);
			for (int j = 0; j < np; j++) {
				m[j][0][0] += coefs[i].x * dx[j];
				m[j][0][1] += coefs[i].x * dy[j];
				m[j][0][2] += coefs[i].x * dz[j];
				m[j][1][0] += coefs[i].y * dx[j];
				m[j][1][1] += coefs[i].y * dy[j];
				m[j][1][2] += coefs[i].y * dz[j];
				m[j][2][0] += coefs[i].z * dx[j];
				m[j][2][1] += coefs[i].z * dy[j];
				m[j][2][2] += coefs[i].z * dz[j];
			}
		}
	}

	return m;
}

double *RefMap::get_jacobian(const int np, const QuadPt3D *pt, bool trans) {
	_F_

	double *jac = new double[np]; MEM_CHECK(jac);
	if (is_const_jacobian) {
		if (trans)
			for (int i = 0; i < np; i++)
				jac[i] = const_jacobian * pt[i].w;
		else
			for (int i = 0; i < np; i++)
				jac[i] = const_jacobian;
	}
	else {
		double3x3 *m = get_ref_map(np, pt);

		double trj = get_transform_jacobian();
		if (trans)
			for (int i = 0; i < np; i++)
				jac[i] = det(m[i]) * trj * pt[i].w;
		else
			for (int i = 0; i < np; i++)
				jac[i] = det(m[i]) * trj;

		delete [] m;
	}

	return jac;
}

double3x3 *RefMap::get_inv_ref_map(const int np, const QuadPt3D *pt) {
	_F_

	double3x3 *irm = new double3x3[np]; MEM_CHECK(irm);

	if (is_const_jacobian) {
		for (int i = 0; i < np; i++)
			memcpy(irm + i, const_inv_ref_map, sizeof(double3x3));
	}
	else {
		double3x3 *m = get_ref_map(np, pt);

		double trj = get_transform_jacobian();
		double *jac = new double[np];
		MEM_CHECK(jac);
		for (int i = 0; i < np; i++) {
			jac[i] = det(m[i]);

			double ij = 1.0 / jac[i];
			irm[i][0][0] = (m[i][1][1] * m[i][2][2] - m[i][1][2] * m[i][2][1]) * ij;
			irm[i][1][0] = (m[i][0][2] * m[i][2][1] - m[i][0][1] * m[i][2][2]) * ij;
			irm[i][2][0] = (m[i][0][1] * m[i][1][2] - m[i][0][2] * m[i][1][1]) * ij;
			irm[i][0][1] = (m[i][1][2] * m[i][2][0] - m[i][1][0] * m[i][2][2]) * ij;
			irm[i][1][1] = (m[i][0][0] * m[i][2][2] - m[i][0][2] * m[i][2][0]) * ij;
			irm[i][2][1] = (m[i][0][2] * m[i][1][0] - m[i][0][0] * m[i][1][2]) * ij;
			irm[i][0][2] = (m[i][1][0] * m[i][2][1] - m[i][1][1] * m[i][2][0]) * ij;
			irm[i][1][2] = (m[i][0][1] * m[i][2][0] - m[i][0][0] * m[i][2][1]) * ij;
			irm[i][2][2] = (m[i][0][0] * m[i][1][1] - m[i][0][1] * m[i][1][0]) * ij;

			jac[i] *= trj;
		}

		delete [] m;
		delete [] jac;
	}

	return irm;
}

double *RefMap::get_phys_x(const int np, const QuadPt3D *pt) {
	_F_
	// transform all x coordinates of the integration points
	double *x = new double[np];
	MEM_CHECK(x);
	memset(x, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < n_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->precalculate(np, pt, FN_DEFAULT);
		blas_axpy(np, coefs[i].x, pss->get_fn_values(), 1, x, 1);
	}

	return x;
}

double *RefMap::get_phys_y(const int np, const QuadPt3D *pt) {
	_F_
	// transform all y coordinates of the integration points
	double *y = new double[np];
	MEM_CHECK(y);
	memset(y, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < n_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->precalculate(np, pt, FN_DEFAULT);
		blas_axpy(np, coefs[i].y, pss->get_fn_values(), 1, y, 1);
	}

	return y;
}

double *RefMap::get_phys_z(const int np, const QuadPt3D *pt) {
	_F_
	// transform all z coordinates of the integration points
	double *z = new double[np];
	MEM_CHECK(z);
	memset(z, 0, np * sizeof(double));
	pss->force_transform(sub_idx, ctm);
	for (int i = 0; i < n_coefs; i++) {
		pss->set_active_shape(indices[i]);
		pss->precalculate(np, pt, FN_DEFAULT);
		blas_axpy(np, coefs[i].z, pss->get_fn_values(), 1, z, 1);
	}

	return z;
}

// TODO: refactor me, this approach is insane
double *RefMap::get_face_jacobian(int iface, const int np, const QuadPt3D *pt, bool trans) {
	_F_
	assert(mesh != NULL);

	double *jac = new double[np]; MEM_CHECK(jac);
	if (is_const_jacobian) {
		double face_jac = calc_face_const_jacobian(iface);
		if (trans)
			for (int i = 0; i < np; i++)
				jac[i] = face_jac * pt[i].w;
		else
			for (int i = 0; i < np; i++)
				jac[i] = face_jac;
	}
	else {
		double3x3 *m = new double3x3[np]; MEM_CHECK(m);
		memset(m, 0, np * sizeof(double3x3));
		const int *face_vertices = RefHex::get_face_vertices(iface);

		for (int i = 0; i < RefHex::get_num_face_vertices(iface); i++) {
			double *dx, *dy, *dz;
			pss->set_active_shape(indices[face_vertices[i]]);
			pss->precalculate(np, pt, FN_DEFAULT);
			pss->get_dx_dy_dz_values(dx, dy, dz);
			for (int j = 0; j < np; j++) {
				m[j][0][0] += vertex[face_vertices[i]].x * dx[j];
				m[j][0][1] += vertex[face_vertices[i]].x * dy[j];
				m[j][0][2] += vertex[face_vertices[i]].x * dz[j];
				m[j][1][0] += vertex[face_vertices[i]].y * dx[j];
				m[j][1][1] += vertex[face_vertices[i]].y * dy[j];
				m[j][1][2] += vertex[face_vertices[i]].y * dz[j];
				m[j][2][0] += vertex[face_vertices[i]].z * dx[j];
				m[j][2][1] += vertex[face_vertices[i]].z * dy[j];
				m[j][2][2] += vertex[face_vertices[i]].z * dz[j];
			}
		}

		// what vectors to take to make vector product
		int dir1, dir2;
		switch (iface) {
			case 0:
			case 1:
				dir1 = 1; dir2 = 2;
				break;

			case 2:
			case 3:
				dir1 = 0; dir2 = 2;
				break;

			case 4:
			case 5:
				dir1 = 0; dir2 = 1;
				break;
		}

		// in fact, it is not jacobian, but vector product
		for (int i = 0; i < np; i++) {
			Point3D vec1 = { m[i][0][dir1], m[i][1][dir1], m[i][2][dir1] };
			Point3D vec2 = { m[i][0][dir2], m[i][1][dir2], m[i][2][dir2] };
			Point3D vec_product = cross_product(vec1, vec2);
			jac[i] = norm(vec_product);
			if (trans) jac[i] *= pt[i].w;
		}

		delete [] m;
	}

	return jac;
}



void RefMap::calc_face_normal(int iface, const int np, const QuadPt3D *pt, double *&nx, double *&ny, double *&nz) {
	_F_
	assert(mesh != NULL);

	double3x3 *m = get_ref_map(np, pt);

	nx = new double[np]; MEM_CHECK(nx);
	ny = new double[np]; MEM_CHECK(ny);
	nz = new double[np]; MEM_CHECK(nz);

	// FIXME: refactor this!
	// - it would be nice if calculation of normals would be based on the same algorithm for all elements
	//   e.g. to take a normal from ref. domain and transform it to the physical one
	int t_dir_1, t_dir_2; //directions of tangents ot the reference face such that t_dir_1 x t_dir_2 = outer normal
	switch (element->get_mode()) {
		case MODE_TETRAHEDRON: {
			const int *face_vtx = element->get_face_vertices(iface);
			Vertex vtx[Tri::NUM_VERTICES];
			for (int i = 0; i < element->get_num_vertices(); i++)
				vtx[i] = vertex[face_vtx[i]];

			Point3D v1 = { vtx[1].x - vtx[0].x, vtx[1].y - vtx[0].y, vtx[1].z - vtx[0].z };
			Point3D v2 = { vtx[2].x - vtx[0].x, vtx[2].y - vtx[2].y, vtx[2].z - vtx[0].z };
			Point3D n = normalize(cross_product(v1, v2));

			for (int i = 0; i < np; i++) {
				nx[i] = n.x;
				ny[i] = n.y;
				nz[i] = n.z;
			}
			} break;

		case MODE_HEXAHEDRON:
			switch (iface) {
				case 0: t_dir_1 = 2; t_dir_2 = 1; break;
				case 1: t_dir_1 = 1; t_dir_2 = 2; break;
				case 2: t_dir_1 = 0; t_dir_2 = 2; break;
				case 3: t_dir_1 = 2; t_dir_2 = 0; break;
				case 4: t_dir_1 = 1; t_dir_2 = 0; break;
				case 5: t_dir_1 = 0; t_dir_2 = 1; break;
			}
			for (int i = 0; i < np; i++) {
				Point3D tangent1 = { m[i][0][t_dir_1], m[i][1][t_dir_1], m[i][2][t_dir_1] };
				Point3D tangent2 = { m[i][0][t_dir_2], m[i][1][t_dir_2], m[i][2][t_dir_2] };
				Point3D normal = normalize(cross_product(tangent1, tangent2));

				nx[i] = normal.x;
				ny[i] = normal.y;
				nz[i] = normal.z;
			}
			break;

		case MODE_PRISM:
			EXIT(H3D_ERR_NOT_IMPLEMENTED);
	}

	delete [] m;
}


void RefMap::calc_const_inv_ref_map() {
	_F_
	// for linear tetrahedra only
	// TODO: does not take in account the transformation (we do not have it for tetras)

	double3x3 m = {
		{ (vertex[1].x - vertex[0].x) / 2, (vertex[2].x - vertex[0].x) / 2, (vertex[3].x - vertex[0].x) / 2 },
		{ (vertex[1].y - vertex[0].y) / 2, (vertex[2].y - vertex[0].y) / 2, (vertex[3].y - vertex[0].y) / 2 },
		{ (vertex[1].z - vertex[0].z) / 2, (vertex[2].z - vertex[0].z) / 2, (vertex[3].z - vertex[0].z) / 2 }
	};
	memcpy(&const_ref_map, &m, sizeof(double3x3));

	const_jacobian =
		m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1] -
		m[2][0] * m[1][1] * m[0][2] - m[2][1] * m[1][2] * m[0][0] - m[2][2] * m[1][0] * m[0][1];

	double ij = 1.0 / const_jacobian;

	const_inv_ref_map[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * ij;
	const_inv_ref_map[1][0] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * ij;
	const_inv_ref_map[2][0] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * ij;
	const_inv_ref_map[0][1] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * ij;
	const_inv_ref_map[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * ij;
	const_inv_ref_map[2][1] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * ij;
	const_inv_ref_map[0][2] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * ij;
	const_inv_ref_map[1][2] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * ij;
	const_inv_ref_map[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * ij;
}


// TODO: rewrite in similar way for non-constant functions
double RefMap::calc_face_const_jacobian(int face) {
	_F_
	// this calculates jacobian on faces of tetrahedron
	// physical triangle
	const int *face_vtx = element->get_face_vertices(face);
	Vertex vtx[Tri::NUM_VERTICES];
	for (int i = 0; i < element->get_num_vertices(); i++)
		vtx[i] = vertex[face_vtx[i]];

	double3x3 m = {
		{ (vtx[1].x - vtx[0].x), (vtx[2].x - vtx[0].x), 1.0 },
		{ (vtx[1].y - vtx[0].y), (vtx[2].y - vtx[0].y), 1.0 },
		{ (vtx[1].z - vtx[0].z), (vtx[2].z - vtx[0].z), 1.0 }
	};
	double phys_triangle_surface = 0.5 *
		sqrt(
			sqr(m[1][0] * m[2][1] - m[1][1] * m[2][0]) +
			sqr(m[0][0] * m[2][1] - m[0][1] * m[2][0]) +
			sqr(m[0][0] * m[1][1] - m[0][1] * m[1][0]));

	// reference triangle from a tetrahedron
	const int *ref_idx = RefTetra::get_face_vertices(face);
	const Point3D *ref_vtx = RefTetra::get_vertices();

	double3x3 n = {
		{ ref_vtx[ref_idx[1]].x - ref_vtx[ref_idx[0]].x, ref_vtx[ref_idx[2]].x - ref_vtx[ref_idx[0]].x, 1.0 },
		{ ref_vtx[ref_idx[1]].y - ref_vtx[ref_idx[0]].y, ref_vtx[ref_idx[2]].y - ref_vtx[ref_idx[0]].y, 1.0 },
		{ ref_vtx[ref_idx[1]].z - ref_vtx[ref_idx[0]].z, ref_vtx[ref_idx[2]].z - ref_vtx[ref_idx[0]].z, 1.0 }
	};

	double ref_triangle_surface = 0.5 *
		sqrt(
			sqr(n[1][0] * n[2][1] - n[1][1] * n[2][0]) +
			sqr(n[0][0] * n[2][1] - n[0][1] * n[2][0]) +
			sqr(n[0][0] * n[1][1] - n[0][1] * n[1][0]));


	return phys_triangle_surface / ref_triangle_surface;
}
