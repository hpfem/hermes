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
// testing linear indepenedence of a shape functions
//

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// TODO: Mostly the same as for H1 => review the code

// l2 product
double l2_product(RealFunction *fu, RealFunction *fv) {
	_F_
	Quad3D *quad = get_quadrature(MODE);

	// integrate with maximum order
	order3_t o = fu->get_fn_order() + fv->get_fn_order() + order3_t(2, 2, 2);
	o.limit();

	int np = quad->get_num_points(o);
	QuadPt3D *pt = quad->get_points(o);
	fu->precalculate(np, pt, FN_DEFAULT);
	fv->precalculate(np, pt, FN_DEFAULT);

	scalar *u0, *u1, *u2;
	u0 = fu->get_fn_values(0);
	u1 = fu->get_fn_values(1);
	u2 = fu->get_fn_values(2);

	scalar *v0, *v1, *v2;
	v0 = fv->get_fn_values(0);
	v1 = fv->get_fn_values(1);
	v2 = fv->get_fn_values(2);

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	double result = 0.0;
	for (int i = 0; i < np; i++)
		result += pt[i].w * (REAL(sqr(u0[i] - v0[i]) + sqr(u1[i] - v1[i]) + sqr(u2[i] - v2[i])));

	return result;
}

bool test_lin_indep(Shapeset *shapeset) {
	_F_
	printf("I. linear independency\n");

	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(&mat, &rhs);

	ShapeFunction pss_u(shapeset), pss_v(shapeset);
	int n = Hex::NUM_EDGES * shapeset->get_num_edge_fns(H3D_MAX_ELEMENT_ORDER)
		+ Hex::NUM_FACES * shapeset->get_num_face_fns(order2_t(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER))
	    + shapeset->get_num_bubble_fns(order3_t(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER));

	printf("number of functions = %d\n", n);

	int *fn_idx = new int[n];
	int m = 0;
	// edge fns
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		int order = H3D_MAX_ELEMENT_ORDER;
		int *edge_idx = shapeset->get_edge_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_edge_fns(order); j++, m++)
			fn_idx[m] = edge_idx[j];
	}
	// face fns
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		order2_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
		int *face_idx = shapeset->get_face_indices(i, 0, order);
		for (int j = 0; j < shapeset->get_num_face_fns(order); j++, m++)
			fn_idx[m] = face_idx[j];
	}
	// bubble
	order3_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
	int *bubble_idx = shapeset->get_bubble_indices(order);
	for (int j = 0; j < shapeset->get_num_bubble_fns(order); j++, m++)
		fn_idx[m] = bubble_idx[j];


	// precalc structure
	mat.prealloc(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			mat.pre_add_ij(i, j);
	mat.alloc();
	rhs.alloc(n);

	printf("assembling matrix ");

	for (int i = 0; i < n; i++) {
		pss_u.set_active_shape(fn_idx[i]);

		printf(".");
		fflush(stdout); // prevent caching of output (to see that it did not freeze)

		for (int j = 0; j < n; j++) {
			pss_v.set_active_shape(fn_idx[j]);

			double value = l2_product(&pss_u, &pss_v);

			mat.add(i, j, value);
		}
	}
	printf("\n");

	for (int i = 0; i < n; i++)
		rhs.add(i, 0.0);

	printf("solving matrix\n");


	// solve the system
	if (solver.solve()) {
		double *sln = solver.get_solution();
		bool indep = true;
		for (int i = 1; i < n + 1; i++) {
			if (sln[i] >= EPS) {
				indep = false;
				break;
			}
		}

		if (indep)
			printf("ok\n");
		else
			printf("Shape functions are not linearly independent\n");
	}
	else {
		printf("Shape functions are not linearly independent\n");
	}

	delete[] fn_idx;

	return true;
}
