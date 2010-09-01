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

#include "../h3dconfig.h"
#include "../function.h"
#include "../solution.h"
#include "h1proj.h"
#include "../matrix.h"
#include "../quad.h"
#include "../refdomain.h"
#include "../transform.h"
#include "../shapeset/common.h"
#include "../shapeset/lobatto.h"
#include <common/callstack.h>

#ifdef DEBUG
	#define PRINTF			printf
#else
	#define PRINTF(...)
#endif

bool H1Projection::has_prods = false;
double H1Projection::prod_fn[N_FNS][N_FNS];
double H1Projection::prod_dx[N_FNS][N_FNS];

void H1Projection::precalc_fn_prods(double fn[N_FNS][N_FNS])
{
	Quad1D *quad_1d = get_quadrature_1d();

	int order = H3D_MAX_QUAD_ORDER;
	QuadPt1D *pt = quad_1d->get_points(order);
	int np = quad_1d->get_num_points(order);

	for (int i = 0; i < N_FNS; i++) {
		shape_fn_1d_t fni = lobatto_fn_tab_1d[i];
		for (int j = 0; j < N_FNS; j++) {
			shape_fn_1d_t fnj = lobatto_fn_tab_1d[j];
			double val = 0.0;
			for (int k = 0; k < np; k++)
				val += pt[k].w * fni(pt[k].x) * fnj(pt[k].x);
			fn[i][j] = val;
		}
	}
}

void H1Projection::precalc_dx_prods(double dx[N_FNS][N_FNS])
{
	Quad1D *quad_1d = get_quadrature_1d();

	int order = H3D_MAX_QUAD_ORDER;
	QuadPt1D *pt = quad_1d->get_points(order);
	int np = quad_1d->get_num_points(order);

	for (int i = 0; i < N_FNS; i++) {
		shape_fn_1d_t fni = lobatto_der_tab_1d[i];
		for (int j = 0; j < N_FNS; j++) {
			shape_fn_1d_t fnj = lobatto_der_tab_1d[j];
			double val = 0.0;
			for (int k = 0; k < np; k++)
				val += pt[k].w * fni(pt[k].x) * fnj(pt[k].x);
			dx[i][j] = val;
		}
	}
}

H1Projection::H1Projection(Solution *afn, Element *e, Shapeset *ss) : Projection(afn, e, ss)
{
	if (!has_prods) {
		precalc_fn_prods(prod_fn);
		precalc_dx_prods(prod_dx);
		has_prods = true;
	}
}

double H1Projection::get_error(int split, int son, const order3_t &order)
{
	_F_
	sln->enable_transform(false);

	calc_projection(split, son + 1, order);

	order3_t order_rhs = order;
	QuadPt3D *pt = quad->get_points(order_rhs);
	int np = quad->get_num_points(order_rhs);

	double error = 0.0;
	for (int i = 0; i < int_ns[split]; i++) {
		Trf *tr = get_trf(int_trf[split][i]);

		Word_t son_idx = base_elem->get_son(int_son[son + 1][i]);
		sln->set_active_element(mesh->elements[son_idx]);
		sln->precalculate(np, pt, FN_DEFAULT);
		scalar *rval = sln->get_fn_values();
		scalar *rdx, *rdy, *rdz;
		sln->get_dx_dy_dz_values(rdx, rdy, rdz);

		QuadPt3D tpt[np];
		transform_points(np, pt, tr, tpt);
		scalar prfn[np], prdx[np], prdy[np], prdz[np];
		memset(prfn, 0, np * sizeof(double));
		memset(prdx, 0, np * sizeof(double));
		memset(prdy, 0, np * sizeof(double));
		memset(prdz, 0, np * sizeof(double));

		for (int i = 0; i < n_fns; i++) {
#ifndef H3D_COMPLEX
			double tmp[np];
			ss->get_fn_values(fn_idx[i], np, tpt, 0, tmp);
			blas_axpy(np, proj_coef[i], tmp, 1, prfn, 1);
			ss->get_dx_values(fn_idx[i], np, tpt, 0, tmp);
			blas_axpy(np, proj_coef[i], tmp, 1, prdx, 1);
			ss->get_dy_values(fn_idx[i], np, tpt, 0, tmp);
			blas_axpy(np, proj_coef[i], tmp, 1, prdy, 1);
			ss->get_dz_values(fn_idx[i], np, tpt, 0, tmp);
			blas_axpy(np, proj_coef[i], tmp, 1, prdz, 1);
#else
			double tmp[np];
			scalar sctmp[np];
			ss->get_fn_values(fn_idx[i], np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj_coef[i], sctmp, 1, prfn, 1);
			ss->get_dx_values(fn_idx[i], np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj_coef[i], sctmp, 1, prdx, 1);
			ss->get_dy_values(fn_idx[i], np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj_coef[i], sctmp, 1, prdy, 1);
			ss->get_dz_values(fn_idx[i], np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj_coef[i], sctmp, 1, prdz, 1);
#endif
		}

		for (int k = 0; k < np; k++)
			error += pt[k].w *
				(sqr(rval[k] - prfn[k]) +
				 sqr(rdx[k] * mdx[split] - prdx[k]) +
				 sqr(rdy[k] * mdy[split] - prdy[k]) +
				 sqr(rdz[k] * mdz[split] - prdz[k]));
	}

	sln->enable_transform(true);

	return error;
}

void H1Projection::calc_projection(int split, int son, const order3_t &order)
{
	_F_

	n_fns = (order.x + 1) * (order.y + 1) * (order.z + 1);

	delete [] fn_idx;
	fn_idx = new int [n_fns];
	int mm = 0;
	// vertex functions
	for (int vtx = 0; vtx < Hex::NUM_VERTICES; vtx++, mm++)
		fn_idx[mm] = ss->get_vertex_index(vtx);
	// edge functions
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		order1_t edge_order = order.get_edge_order(iedge);
		int n_edge_fns = ss->get_num_edge_fns(edge_order);
		if (n_edge_fns > 0) {
			const int *edge_fn_idx = ss->get_edge_indices(iedge, 0, edge_order);
			for (int i = 0; i < n_edge_fns; i++, mm++)
				fn_idx[mm] = edge_fn_idx[i];
		}
	}
	// face functions
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		order2_t face_order = order.get_face_order(iface);
		int n_face_fns = ss->get_num_face_fns(face_order);
		if (n_face_fns > 0) {
			const int *face_fn_idx = ss->get_face_indices(iface, 0, face_order);
			for (int i = 0; i < n_face_fns; i++, mm++)
				fn_idx[mm] = face_fn_idx[i];
		}
	}
	{
		// bubble functions
		int n_bubble_fns = ss->get_num_bubble_fns(order);
		if (n_bubble_fns > 0) {
			const int *bubble_fn_idx = ss->get_bubble_indices(order);
			for (int i = 0; i < n_bubble_fns; i++, mm++)
				fn_idx[mm] = bubble_fn_idx[i];
		}
	}

	double **proj_mat = new_matrix<double>(n_fns, n_fns);
	scalar *proj_rhs = new scalar[n_fns];
	memset(proj_rhs, 0, sizeof(scalar) * n_fns);

	// proj matrix
	for (int i = 0; i < n_fns; i++) {
		int iidx = fn_idx[i];
		order3_t oi = ss->get_dcmp(iidx);

		for (int j = 0; j < n_fns; j++) {
			int jidx = fn_idx[j];

			order3_t oj = ss->get_dcmp(jidx);
			double val =
				prod_fn[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_dx[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_fn[oi.x][oj.x] * prod_dx[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_fn[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_dx[oi.z][oj.z];
			proj_mat[i][j] += val;
		}
	}

	// rhs
	for (int e = 0; e < int_ns[split]; e++) {
		Word_t son_idx = base_elem->get_son(int_son[son][e]);
		sln->set_active_element(mesh->elements[son_idx]);

		Trf *tr = get_trf(int_trf[split][e]);
		for (int i = 0; i < n_fns; i++) {
			int iidx = fn_idx[i];
			fu->set_active_shape(iidx);

			order3_t order_rhs = ss->get_order(iidx) + order;
			QuadPt3D *pt = quad->get_points(order_rhs);
			int np = quad->get_num_points(order_rhs);

			if (int_trf[split][e] != -1) fu->push_transform(int_trf[split][e]);
			fu->precalculate(np, pt, FN_DEFAULT);
			sln->precalculate(np, pt, FN_DEFAULT);

			double *uval = fu->get_fn_values();
			scalar *rval = sln->get_fn_values();

			double *dudx, *dudy, *dudz;
			scalar *drdx, *drdy, *drdz;

			fu->get_dx_dy_dz_values(dudx, dudy, dudz);
			sln->get_dx_dy_dz_values(drdx, drdy, drdz);

			QuadPt3D tpt[np];
			transform_points(np, pt, tr, tpt);

			scalar value = 0.0;
			for (int k = 0; k < np; k++) {
				value += pt[k].w * (uval[k] * rval[k] +
						dudx[k] * drdx[k] * mdx[split] +
						dudy[k] * drdy[k] * mdy[split] +
						dudz[k] * drdz[k] * mdz[split]);
			}
			proj_rhs[i] += value * (1 / (double) int_ns[split]);

			if (int_trf[split][e] != -1) fu->pop_transform();
		}
	}

	double d;
	int iperm[n_fns];
	ludcmp(proj_mat, n_fns, iperm, &d);
	lubksb(proj_mat, n_fns, iperm, proj_rhs);

	proj_coef = new double [n_fns];
	memcpy(proj_coef, proj_rhs, n_fns * sizeof(double));

	delete [] proj_mat;
	delete [] proj_rhs;
}
