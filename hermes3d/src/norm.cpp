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
#include "norm.h"
#include "quad.h"
#include "refmap.h"
#include "traverse.h"
#include <common/callstack.h>

#define H3D_TINY			10e-15

#define INTEGRATE_EXPRESSION(exp) \
	double *jac = ru->get_jacobian(np, pt); \
	for (int i = 0; i < np; i++) \
		result += jac[i] * (exp); \
	delete [] jac;

/// Calculates the absolute error between sln1 and sln2 using function fn
double calc_error(double (*fn)(MeshFunction*, MeshFunction*, int, QuadPt3D*), MeshFunction *sln1, MeshFunction *sln2) {
	_F_
	Mesh *meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
	Transformable *tr[2] = { sln1, sln2 };
	Traverse trav;
	trav.begin(2, meshes, tr);

	double error = 0.0;
	Element **ee;
	while ((ee = trav.get_next_state(NULL, NULL)) != NULL) {
		EMode3D mode = ee[0]->get_mode();

		RefMap *ru = sln1->get_refmap();
		order3_t order = max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
		order.limit();

		Quad3D *quad = get_quadrature(mode);
		int np = quad->get_num_points(order);
		QuadPt3D *pt = quad->get_points(order);

		error += fn(sln1, sln2, np, pt);
	}
	trav.finish();

	return error > H3D_TINY ? sqrt(error) : error;		// do not ruin the precision by taking the sqrt
}

/// Calculates the norm of sln using function fn
double calc_norm(double (*fn)(MeshFunction*, int, QuadPt3D*), MeshFunction *sln) {
	_F_
	double norm = 0.0;
	Mesh *mesh = sln->get_mesh();

	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];
		sln->set_active_element(e);

		RefMap *ru = sln->get_refmap();
		order3_t o = sln->get_fn_order() + ru->get_inv_ref_order();
		o.limit();

		Quad3D *quad = get_quadrature(e->get_mode());
		int np = quad->get_num_points(o);
		QuadPt3D *pt = quad->get_points(o);

		norm += fn(sln, np, pt);
	}

	return norm > H3D_TINY ? sqrt(norm) : norm;			// do not ruin the precision by taking the sqrt
}

// H1 space /////////////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in H1 norm
double error_fn_h1(MeshFunction *sln1, MeshFunction *sln2, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln1->get_refmap();
	sln1->precalculate(np, pt, FN_DEFAULT);
	sln2->precalculate(np, pt, FN_DEFAULT);

	scalar *uval, *vval, *dudx, *dudy, *dudz, *dvdx, *dvdy, *dvdz;
	uval = sln1->get_fn_values();
	vval = sln2->get_fn_values();
	sln1->get_dx_dy_dz_values(dudx, dudy, dudz);
	sln2->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	double result = 0.0;
	INTEGRATE_EXPRESSION(sqr(uval[i] - vval[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]) + sqr(dudz[i] - dvdz[i]));
	return result;
}

// function used to calculate H1 norm of the solution
double norm_fn_h1(MeshFunction *sln, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln->get_refmap();
	sln->precalculate(np, pt, FN_DEFAULT);

	scalar *uval = sln->get_fn_values();
	scalar *dudx = sln->get_dx_values();
	scalar *dudy = sln->get_dy_values();
	scalar *dudz = sln->get_dz_values();

	double result = 0.0;
	INTEGRATE_EXPRESSION(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]) + sqr(dudz[i]));
	return result;
}


double h1_error(MeshFunction *sln1, MeshFunction *sln2) {
	_F_
	double error = calc_error(error_fn_h1, sln1, sln2);
	double norm = calc_norm(norm_fn_h1, sln2);
	if (norm > H3D_TINY) return error / norm;
	else return error;
}

double h1_norm(MeshFunction *sln) {
	_F_
	return calc_norm(norm_fn_h1, sln);
}

// function used to calculate error in L2 norm
double error_fn_l2(MeshFunction *sln1, MeshFunction *sln2, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln1->get_refmap();
	sln1->precalculate(np, pt, FN_DEFAULT);
	sln2->precalculate(np, pt, FN_DEFAULT);

	scalar *uval, *vval;
	uval = sln1->get_fn_values();
	vval = sln2->get_fn_values();

	double result = 0.0;
	INTEGRATE_EXPRESSION(sqr(uval[i] - vval[i]));
	return result;
}


// function used to calculate L2 norm of the solution
double norm_fn_l2(MeshFunction *sln, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln->get_refmap();
	sln->precalculate(np, pt, FN_DEFAULT);

	scalar *uval = sln->get_fn_values();

	double result = 0.0;
	INTEGRATE_EXPRESSION(sqr(uval[i]));
	return result;
}


double l2_error(MeshFunction *sln1, MeshFunction *sln2) {
	_F_
	double error = calc_error(error_fn_l2, sln1, sln2);
	double norm = calc_norm(norm_fn_l2, sln2);
	if (norm > H3D_TINY) return error / norm;
	else return error;
}

double l2_norm(MeshFunction *sln) {
	_F_
	return calc_norm(norm_fn_l2, sln);
}


// Hcurl space /////////////////////////////////////////////////////////////////////////////////////////

#define U_CURL_0 (du2dy[i] - du1dz[i])
#define U_CURL_1 (du0dz[i] - du2dx[i])
#define U_CURL_2 (du1dx[i] - du0dy[i])

#define V_CURL_0 (dv2dy[i] - dv1dz[i])
#define V_CURL_1 (dv0dz[i] - dv2dx[i])
#define V_CURL_2 (dv1dx[i] - dv0dy[i])


// function used to calculate error in HCurl norm
double error_fn_hcurl(MeshFunction *sln1, MeshFunction *sln2, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln1->get_refmap();

	sln1->precalculate(np, pt, FN_DEFAULT);
	sln2->precalculate(np, pt, FN_DEFAULT);

	scalar *u0, *u1, *u2, *du0dx, *du0dy, *du0dz, *du1dx, *du1dy, *du1dz, *du2dx, *du2dy, *du2dz;
	u0 = sln1->get_fn_values(0);
	u1 = sln1->get_fn_values(1);
	u2 = sln1->get_fn_values(2);
	sln1->get_dx_dy_dz_values(du0dx, du0dy, du0dz, 0);
	sln1->get_dx_dy_dz_values(du1dx, du1dy, du1dz, 1);
	sln1->get_dx_dy_dz_values(du2dx, du2dy, du2dz, 2);

	scalar *v0, *v1, *v2, *dv0dx, *dv0dy, *dv0dz, *dv1dx, *dv1dy, *dv1dz, *dv2dx, *dv2dy, *dv2dz;
	v0 = sln2->get_fn_values(0);
	v1 = sln2->get_fn_values(1);
	v2 = sln2->get_fn_values(2);
	sln2->get_dx_dy_dz_values(dv0dx, dv0dy, dv0dz, 0);
	sln2->get_dx_dy_dz_values(dv1dx, dv1dy, dv1dz, 1);
	sln2->get_dx_dy_dz_values(dv2dx, dv2dy, dv2dz, 2);

	scalar result = 0.0;
	INTEGRATE_EXPRESSION(sqr(u0[i] - v0[i]) + sqr(u1[i] - v1[i]) + sqr(u2[i] - v2[i]) +
	                     sqr(U_CURL_0 - V_CURL_0) + sqr(U_CURL_1 - V_CURL_1) + sqr(U_CURL_2 - V_CURL_2));
	return REAL(result);
}

// function used to calculate HCurl norm of the solution
double norm_fn_hcurl(MeshFunction *sln, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln->get_refmap();

	sln->precalculate(np, pt, FN_DEFAULT);

	scalar *u0, *u1, *u2, *du0dx, *du0dy, *du0dz, *du1dx, *du1dy, *du1dz, *du2dx, *du2dy, *du2dz;
	u0 = sln->get_fn_values(0);
	u1 = sln->get_fn_values(1);
	u2 = sln->get_fn_values(2);
	sln->get_dx_dy_dz_values(du0dx, du0dy, du0dz, 0);
	sln->get_dx_dy_dz_values(du1dx, du1dy, du1dz, 1);
	sln->get_dx_dy_dz_values(du2dx, du2dy, du2dz, 2);

	scalar result = 0.0;
	INTEGRATE_EXPRESSION(sqr(u0[i]) + sqr(u1[i]) + sqr(u2[i]) + sqr(U_CURL_0) + sqr(U_CURL_1) + sqr(U_CURL_2));
	return REAL(result);
}


double hcurl_error(MeshFunction *sln1, MeshFunction *sln2) {
	_F_
	double error = calc_error(error_fn_hcurl, sln1, sln2);
	double norm = calc_norm(norm_fn_hcurl, sln2);
	if (norm > H3D_TINY) return error / norm;
	else return error;
}

double hcurl_norm(MeshFunction *sln) {
	_F_
	return calc_norm(norm_fn_hcurl, sln);
}

// function used to calculate error in L2 norm
double error_fn_l2_hcurl(MeshFunction *sln1, MeshFunction *sln2, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln1->get_refmap();

	sln1->precalculate(np, pt, FN_DEFAULT);
	sln2->precalculate(np, pt, FN_DEFAULT);

	scalar *u0 = sln1->get_fn_values(0);
	scalar *u1 = sln1->get_fn_values(1);
	scalar *u2 = sln1->get_fn_values(2);

	scalar *v0 = sln2->get_fn_values(0);
	scalar *v1 = sln2->get_fn_values(1);
	scalar *v2 = sln2->get_fn_values(2);

	scalar result = 0.0;
	INTEGRATE_EXPRESSION(sqr(u0[i] - v0[i]) + sqr(u1[i] - v1[i]) + sqr(u2[i] - v2[i]));
	return REAL(result);
}


// function used to calculate L2 norm of the solution
double norm_fn_l2_hcurl(MeshFunction *sln, int np, QuadPt3D *pt) {
	_F_

	RefMap *ru = sln->get_refmap();

	sln->precalculate(np, pt, FN_DEFAULT);

	scalar *u0 = sln->get_fn_values(0);
	scalar *u1 = sln->get_fn_values(1);
	scalar *u2 = sln->get_fn_values(2);

	scalar result = 0.0;
	INTEGRATE_EXPRESSION(sqr(u0[i]) + sqr(u1[i]) + sqr(u2[i]));
	return REAL(result);
}


double l2_error_hcurl(MeshFunction *sln1, MeshFunction *sln2) {
	_F_
	double error = calc_error(error_fn_l2_hcurl, sln1, sln2);
	double norm = calc_norm(norm_fn_l2_hcurl, sln2);
	if (norm > H3D_TINY) return error / norm;
	else return error;
}

double l2_norm_hcurl(MeshFunction *sln) {
	_F_
	return calc_norm(norm_fn_l2_hcurl, sln);
}
