// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
// - Pavel Kus
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
// testing correctness of a gradients
//

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

void hcurl_int_vol(RealFunction *fu, double3 result) {
	_F_
/*	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);

	order3_t o = quad->get_max_order();

	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);

	fu->precalculate(np, pt, FN_DX | FN_DY | FN_DZ);

	double *dx = fu->get_dx_values();
	double *dy = fu->get_dy_values();
	double *dz = fu->get_dz_values();

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	result[0] = result[1] = result[2] = 0.0;
	for (int i = 0; i < np; i++) {
		result[0] += pt[i].w * dx[i];
		result[1] += pt[i].w * dy[i];
		result[2] += pt[i].w * dz[i];
	}
*/
}

void hcurl_int_surf(RealFunction *fu, double3 result) {
	_F_
/*	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);

	Point3D norm[] = {
		{ -1.0,  0.0,  0.0 },
		{  1.0,  0.0,  0.0 },
		{  0.0, -1.0,  0.0 },
		{  0.0,  1.0,  0.0 },
		{  0.0,  0.0, -1.0 },
		{  0.0,  0.0,  1.0 },
	};

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	result[0] = result[1] = result[2] = 0.0;
	for (int face = 0; face < Hex::NUM_FACES; face++) {
		order2_t face_order = quad->get_face_max_order(face);

		QuadPt3D *pt = quad->get_face_points(face, face_order);
		int np = quad->get_face_num_points(face, face_order);

		fu->precalculate(np, pt, FN_VAL);
		double *val = fu->get_fn_values();

		for (int i = 0; i < np; i++) {
			result[0] += pt[i].w * norm[face].x * val[i];
			result[1] += pt[i].w * norm[face].y * val[i];
			result[2] += pt[i].w * norm[face].z * val[i];
		}
	}
*/
}

bool test_grad(int fn_idx, Shapeset *shapeset) {
	_F_
	ShapeFunction pss_u(shapeset);

	pss_u.set_active_shape(fn_idx);

	printf(".");
	fflush(stdout);			// prevent caching of output (to see that it did not freeze)

	double3 vol_val = { 0.0, 0.0, 0.0 };
	double3 surf_val = { 0.0, 0.0, 0.0 };

	hcurl_int_vol(&pss_u, vol_val);
	hcurl_int_surf(&pss_u, surf_val);

	if (fabs(vol_val[0] - surf_val[0]) > EPS || fabs(vol_val[1] - surf_val[1]) > EPS || fabs(vol_val[2] - surf_val[2]) > EPS) {
		printf("\n");
		warning("Gradient values for fn #%d do not match", fn_idx);
		return false;
	}

	return true;
}

bool test_gradients(Shapeset *shapeset) {
	_F_
	printf("IV. gradients\n");

	// edge fns
	printf("\n* Edge functions\n");
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		order1_t order = H3D_MAX_ELEMENT_ORDER;
		int *edge_idx = shapeset->get_edge_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_edge_fns(order); j++) {
			if (!test_grad(edge_idx[j], shapeset))
				return false;
		}
	}

	// face fns
	printf("\n* Face functions\n");
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		order2_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
		int *face_idx = shapeset->get_face_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_face_fns(order); j++) {
			if (!test_grad(face_idx[j], shapeset))
				return false;
		}
	}

	// bubble
	printf("\n* Bubble functions\n");
	order3_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
	int *bubble_idx = shapeset->get_bubble_indices(order);
	for (int j = 0; j < shapeset->get_num_bubble_fns(order); j++) {
		if (!test_grad(bubble_idx[j], shapeset))
			return false;
	}

	printf("\n");

	return true;
}

