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
#include "h1projipol.h"
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

bool H1ProjectionIpol::has_prods = false;
double H1ProjectionIpol::prod_fn[N_FNS][N_FNS];
double H1ProjectionIpol::prod_dx[N_FNS][N_FNS];

H1ProjectionIpol::H1ProjectionIpol(Solution *afn, Element *e, Shapeset *ss) : ProjectionIpol(afn, e, ss)
{
	if (!has_prods) {
		H1Projection::precalc_fn_prods(prod_fn);
		H1Projection::precalc_dx_prods(prod_dx);
		has_prods = true;
	}
}

double H1ProjectionIpol::get_error(int split, int son, const order3_t &order)
{
	_F_
	sln->enable_transform(false);

	calc_projection(split, son, order);

	// error
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

		for (int i = 0; i < proj_fns; i++) {
#ifndef H3D_COMPLEX
			double tmp[np];
			ss->get_fn_values(proj[i]->idx, np, tpt, 0, tmp);
			blas_axpy(np, proj[i]->coef, tmp, 1, prfn, 1);
			ss->get_dx_values(proj[i]->idx, np, tpt, 0, tmp);
			blas_axpy(np, proj[i]->coef, tmp, 1, prdx, 1);
			ss->get_dy_values(proj[i]->idx, np, tpt, 0, tmp);
			blas_axpy(np, proj[i]->coef, tmp, 1, prdy, 1);
			ss->get_dz_values(proj[i]->idx, np, tpt, 0, tmp);
			blas_axpy(np, proj[i]->coef, tmp, 1, prdz, 1);
#else
			double tmp[np];
			scalar sctmp[np];
			ss->get_fn_values(proj[i]->idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj[i]->coef, sctmp, 1, prfn, 1);
			ss->get_dx_values(proj[i]->idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj[i]->coef, sctmp, 1, prdx, 1);
			ss->get_dy_values(proj[i]->idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj[i]->coef, sctmp, 1, prdy, 1);
			ss->get_dz_values(proj[i]->idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, proj[i]->coef, sctmp, 1, prdz, 1);
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

void H1ProjectionIpol::calc_vertex_proj(int split, int son)
{
	_F_
	int nv = base_elem->get_num_vertices();

	vertex_proj = new ProjItem[nv];
	MEM_CHECK(vertex_proj);

	// FIXME: hex specific
	const Point3D *vtx_pt = RefHex::get_vertices();
	for (int ivtx = 0; ivtx < nv; ivtx++) {
		Word_t son_idx = base_elem->get_son(vtx_son[son][ivtx]);
		sln->set_active_element(mesh->elements[son_idx]);
		QuadPt3D pt(vtx_pt[ivtx].x, vtx_pt[ivtx].y, vtx_pt[ivtx].z, 1.0);
		sln->precalculate(1, &pt, FN_VAL);

		scalar *proj = sln->get_fn_values();
		vertex_proj[ivtx].coef = proj[0];
		vertex_proj[ivtx].idx = ss->get_vertex_index(ivtx);
	}
}

void H1ProjectionIpol::calc_edge_proj(int iedge, int split, int son, const order3_t &order)
{
	_F_
	order1_t edge_order = order.get_edge_order(iedge);
	int edge_fns = edge_order - 1;
	if (edge_fns <= 0) return;

	scalar *proj_rhs = new scalar[edge_fns];
	MEM_CHECK(proj_rhs);
	memset(proj_rhs, 0, sizeof(scalar) * edge_fns);
	double **proj_mat = new_matrix<double>(edge_fns, edge_fns);
	MEM_CHECK(proj_rhs);

	// local edge vertex numbers
	const int *edge_vtx = RefHex::get_edge_vertices(iedge);
	ProjItem vtxp[] = { vertex_proj[edge_vtx[0]], vertex_proj[edge_vtx[1]] };

	int *edge_fn_idx = ss->get_edge_indices(iedge, 0, edge_order);	// indices of edge functions
	for (int i = 0; i < edge_fns; i++) {
		int iidx = edge_fn_idx[i];
		order3_t oi = ss->get_dcmp(iidx);
		for (int j = 0; j < edge_fns; j++) {
			int jidx = edge_fn_idx[j];
			order3_t oj = ss->get_dcmp(jidx);
			double val = 0.0;
			if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) {
				val = prod_fn[oi.x][oj.x] + prod_dx[oi.x][oj.x];
			}
			else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) {
				val = prod_fn[oi.y][oj.y] + prod_dx[oi.y][oj.y];
			}
			else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) {
				val = prod_fn[oi.z][oj.z] + prod_dx[oi.z][oj.z];
			}
			else
				EXIT("Local edge number out of range.");
			proj_mat[i][j] += val;
		}
	}

	for (int e = 0; e < edge_ns[split][iedge]; e++) {
		edge_fn_idx = ss->get_edge_indices(iedge, 0, edge_order);	// indices of edge functions

		Word_t son_idx = base_elem->get_son(edge_son[son][iedge][e]);
		sln->set_active_element(mesh->elements[son_idx]);

		Trf *tr = get_trf(edge_trf[split][iedge][e]);
		for (int i = 0; i < edge_fns; i++) {
			int iidx = edge_fn_idx[i];
			fu->set_active_shape(iidx);

			order1_t ord = (ss->get_order(iidx) + order).get_edge_order(iedge);
			QuadPt3D *pt = quad->get_edge_points(iedge, ord);
			int np = quad->get_edge_num_points(iedge, ord);

			if (edge_trf[split][iedge][e] != -1) fu->push_transform(edge_trf[split][iedge][e]);
			fu->precalculate(np, pt, FN_DEFAULT);
			sln->precalculate(np, pt, FN_DEFAULT);

			double *uval = fu->get_fn_values();
			scalar *rval = sln->get_fn_values();

			double *du, md;
			scalar *dr;
			if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) {
				du = fu->get_dx_values();
				dr = sln->get_dx_values();
				md = mdx[split];
			}
			else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) {
				du = fu->get_dy_values();
				dr = sln->get_dy_values();
				md = mdy[split];
			}
			else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) {
				du = fu->get_dz_values();
				dr = sln->get_dz_values();
				md = mdz[split];
			}
			else
				EXIT("Local edge number out of range.");

			QuadPt3D tpt[np];
			transform_points(np, pt, tr, tpt);

			double tmp[np];
			scalar sctmp[np];
			scalar g[np];						// interpolant
			memset(g, 0, np * sizeof(scalar));
#ifndef H3D_COMPLEX
			ss->get_fn_values(vtxp[0].idx, np, tpt, 0, tmp);
			blas_axpy(np, vtxp[0].coef, tmp, 1, g, 1);
			ss->get_fn_values(vtxp[1].idx, np, tpt, 0, tmp);
			blas_axpy(np, vtxp[1].coef, tmp, 1, g, 1);
#else
			ss->get_fn_values(vtxp[0].idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, vtxp[0].coef, sctmp, 1, g, 1);
			ss->get_fn_values(vtxp[1].idx, np, tpt, 0, tmp);
			for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
			blas_axpy(np, vtxp[1].coef, sctmp, 1, g, 1);
#endif

			scalar dg[np];
			memset(dg, 0, np * sizeof(scalar));
			if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) {
				ss->get_dx_values(vtxp[0].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[0].coef, sctmp, 1, dg, 1);
				ss->get_dx_values(vtxp[1].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[1].coef, sctmp, 1, dg, 1);
			}
			else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) {
				ss->get_dy_values(vtxp[0].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[0].coef, sctmp, 1, dg, 1);
				ss->get_dy_values(vtxp[1].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[1].coef, sctmp, 1, dg, 1);
			}
			else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) {
				ss->get_dz_values(vtxp[0].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[0].coef, sctmp, 1, dg, 1);
				ss->get_dz_values(vtxp[1].idx, np, tpt, 0, tmp);
				for (int ii = 0; ii < np; ii++) sctmp[ii] = tmp[ii];
				blas_axpy(np, vtxp[1].coef, sctmp, 1, dg, 1);
			}
			else
				EXIT("Local edge number out of range.");


			scalar value = 0.0;
			for (int k = 0; k < np; k++)
				value += pt[k].w * (uval[k] * (rval[k] - g[k]) + du[k] * ((dr[k] * md) - dg[k]));
			proj_rhs[i] += value * (1 / (double) edge_ns[split][iedge]);

			if (edge_trf[split][iedge][e] != -1) fu->pop_transform();
		}
	}

	double d;
	int iperm[edge_fns];
	ludcmp(proj_mat, edge_fns, iperm, &d);
	lubksb(proj_mat, edge_fns, iperm, proj_rhs);

	// copy functions and coefficients to the basis
	edge_proj[iedge] = new ProjItem[edge_fns];
	for (int i = 0; i < edge_fns; i++) {
		edge_proj[iedge][i].coef = proj_rhs[i];
		edge_proj[iedge][i].idx = edge_fn_idx[i];
	}

	delete [] proj_mat;
	delete [] proj_rhs;
}

void H1ProjectionIpol::calc_face_proj(int iface, int split, int son, const order3_t &order)
{
	_F_
	order2_t face_order = order.get_face_order(iface);
	int face_fns = (face_order.x - 1) * (face_order.y - 1);
	if (face_fns <= 0) return;

	scalar *proj_rhs = new scalar[face_fns];
	MEM_CHECK(proj_rhs);
	memset(proj_rhs, 0, sizeof(scalar) * face_fns);
	double **proj_mat = new_matrix<double>(face_fns, face_fns);
	MEM_CHECK(proj_mat);

	const int *face_vertex = RefHex::get_face_vertices(iface);
	const int *face_edge = RefHex::get_face_edges(iface);

	// get total number of functions for interpolant (vertex + edge functions)
	int ipol_fns = RefHex::get_num_face_vertices(iface);
	for (int iedge = 0; iedge < RefHex::get_num_face_edges(iface); iedge++)
		ipol_fns += order.get_edge_order(face_edge[iedge]) - 1;

	// interpolant
	ProjItem ipol[ipol_fns];
	int mm = 0;
	for (int vtx = 0; vtx < RefHex::get_num_face_vertices(iface); vtx++, mm++)
		ipol[mm] = vertex_proj[face_vertex[vtx]];
	for (int iedge = 0; iedge < RefHex::get_num_face_edges(iface); iedge++) {
		order1_t edge_order = order.get_edge_order(face_edge[iedge]);
		int edge_fns = edge_order - 1;
		for (int i = 0; i < edge_fns; i++, mm++)
			ipol[mm] = edge_proj[face_edge[iedge]][i];
	}

	int face_ori = 0;
	int *face_fn_idx = ss->get_face_indices(iface, face_ori, face_order);
	for (int i = 0; i < face_fns; i++) {
		int iidx = face_fn_idx[i];
		order3_t oi = ss->get_dcmp(iidx);
		for (int j = 0; j < face_fns; j++) {
			int jidx = face_fn_idx[j];
			order3_t oj = ss->get_dcmp(jidx);
			double val = 0.0;
			if (iface == 0 || iface == 1) {
				val =
					prod_fn[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
					prod_dx[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
					prod_fn[oi.y][oj.y] * prod_dx[oi.z][oj.z];
			}
			else if (iface == 2 || iface == 3) {
				val =
					prod_fn[oi.x][oj.x] * prod_fn[oi.z][oj.z] +
					prod_dx[oi.x][oj.x] * prod_fn[oi.z][oj.z] +
					prod_fn[oi.x][oj.x] * prod_dx[oi.z][oj.z];
			}
			else if (iface == 4 || iface == 5) {
				val =
					prod_fn[oi.x][oj.x] * prod_fn[oi.y][oj.y] +
					prod_dx[oi.x][oj.x] * prod_fn[oi.y][oj.y] +
					prod_fn[oi.x][oj.x] * prod_dx[oi.y][oj.y];
			}
			else
				EXIT("Local face number out of range.");
			proj_mat[i][j] += val;
		}
	}

	for (int e = 0; e < face_ns[split][iface]; e++) {
		Word_t son_idx = base_elem->get_son(face_son[son][iface][e]);
		sln->set_active_element(mesh->elements[son_idx]);

		Trf *tr = get_trf(face_trf[split][iface][e]);
		for (int i = 0; i < face_fns; i++) {
			int iidx = face_fn_idx[i];
			fu->set_active_shape(iidx);

			order2_t ord = (ss->get_order(iidx) + order).get_face_order(iface);
			QuadPt3D *pt = quad->get_face_points(iface, ord);
			int np = quad->get_face_num_points(iface, ord);

			if (face_trf[split][iface][e] != -1) fu->push_transform(face_trf[split][iface][e]);
			fu->precalculate(np, pt, FN_DEFAULT);
			sln->precalculate(np, pt, FN_DEFAULT);

			double *uval = fu->get_fn_values();
			scalar *rval = sln->get_fn_values();

			double *dudx, *dudy;
			scalar *drdx, *drdy;
			double md, me;

			if (iface == 0 || iface == 1) {
				dudx = fu->get_dy_values();
				drdx = sln->get_dy_values();
				dudy = fu->get_dz_values();
				drdy = sln->get_dz_values();
				md = mdy[split];
				me = mdz[split];
			}
			else if (iface == 2 || iface == 3) {
				dudx = fu->get_dx_values();
				drdx = sln->get_dx_values();
				dudy = fu->get_dz_values();
				drdy = sln->get_dz_values();
				md = mdx[split];
				me = mdz[split];
			}
			else if (iface == 4 || iface == 5) {
				dudx = fu->get_dx_values();
				drdx = sln->get_dx_values();
				dudy = fu->get_dy_values();
				drdy = sln->get_dy_values();
				md = mdx[split];
				me = mdy[split];
			}
			else
				EXIT("Local face number out of range.");

			QuadPt3D tpt[np];
			transform_points(np, pt, tr, tpt);

			scalar g[np], dgdx[np], dgdy[np];
			memset(g, 0, np * sizeof(scalar));
			memset(dgdx, 0, np * sizeof(scalar));
			memset(dgdy, 0, np * sizeof(scalar));

			for (int l = 0; l < ipol_fns; l++) {
				double h[np];
				scalar sch[np];
				ss->get_fn_values(ipol[l].idx, np, tpt, 0, h);
				for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
				blas_axpy(np, ipol[l].coef, sch, 1, g, 1);

				if (iface == 0 || iface == 1) {
					ss->get_dy_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdx, 1);
					ss->get_dz_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdy, 1);
				}
				else if (iface == 2 || iface == 3) {
					ss->get_dx_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdx, 1);
					ss->get_dz_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdy, 1);
				}
				else if (iface == 4 || iface == 5) {
					ss->get_dx_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdx, 1);
					ss->get_dy_values(ipol[l].idx, np, tpt, 0, h);
					for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
					blas_axpy(np, ipol[l].coef, sch, 1, dgdy, 1);
				}
				else
					EXIT("Local face number out of range.");
			}

			scalar value = 0.0;
			for (int k = 0; k < np; k++)
				value += pt[k].w * (uval[k] * (rval[k] - g[k]) + dudx[k] * ((drdx[k] * md) - dgdx[k]) + dudy[k] * ((drdy[k] * me) - dgdy[k]));
			proj_rhs[i] += value * (1 / (double) face_ns[split][iface]);

			if (face_trf[split][iface][e] != -1) fu->pop_transform();
		}
	}

	double d;
	int iperm[face_fns];
	ludcmp(proj_mat, face_fns, iperm, &d);
	lubksb(proj_mat, face_fns, iperm, proj_rhs);

	face_proj[iface] = new ProjItem [face_fns];
	for (int i = 0; i < face_fns; i++) {
		face_proj[iface][i].coef = proj_rhs[i];
		face_proj[iface][i].idx = face_fn_idx[i];
	}

	delete [] proj_mat;
	delete [] proj_rhs;
}

void H1ProjectionIpol::calc_bubble_proj(int split, int son, const order3_t &order) {
	_F_
	int bubble_fns = (order.x - 1) * (order.y - 1) * (order.z - 1);
	if (bubble_fns <= 0) return;

	scalar *proj_rhs = new scalar[bubble_fns];
	MEM_CHECK(proj_rhs);
	memset(proj_rhs, 0, sizeof(scalar) * bubble_fns);
	double **proj_mat = new_matrix<double>(bubble_fns, bubble_fns);
	MEM_CHECK(proj_mat);

	// get total number of functions (vertex + edge + face)
	int ipol_fns = Hex::NUM_VERTICES;
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		ipol_fns += order.get_edge_order(iedge) - 1;
	}
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		order2_t face_order = order.get_face_order(iface);
		ipol_fns += (face_order.x - 1) * (face_order.y - 1);
	}

	ProjItem ipol[ipol_fns];
	int mm = 0;
	// vertex projection coefficients
	for (int vtx = 0; vtx < Hex::NUM_VERTICES; vtx++, mm++)
		ipol[mm] = vertex_proj[vtx];
	// edge projection coefficients
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		order1_t edge_order = order.get_edge_order(iedge);
		int edge_fns = edge_order - 1;
		for (int i = 0; i < edge_fns; i++, mm++)
			ipol[mm] = edge_proj[iedge][i];
	}
	// face projection coefficients
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		order2_t face_order = order.get_face_order(iface);
		int face_fns = (face_order.x - 1) * (face_order.y - 1);
		for (int i = 0; i < face_fns; i++, mm++)
			ipol[mm] = face_proj[iface][i];
	}

	// do it //
	int *bubble_fn_idx = ss->get_bubble_indices(order);
	for (int i = 0; i < bubble_fns; i++) {
		int iidx = bubble_fn_idx[i];
		order3_t oi = ss->get_dcmp(iidx);
		for (int j = 0; j < bubble_fns; j++) {
			int jidx = bubble_fn_idx[j];
			order3_t oj = ss->get_dcmp(jidx);
			double val =
				prod_fn[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_dx[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_fn[oi.x][oj.x] * prod_dx[oi.y][oj.y] * prod_fn[oi.z][oj.z] +
				prod_fn[oi.x][oj.x] * prod_fn[oi.y][oj.y] * prod_dx[oi.z][oj.z];
			proj_mat[i][j] += val;
		}
	}

	for (int e = 0; e < int_ns[split]; e++) {
		Word_t son_idx = base_elem->get_son(int_son[son][e]);
		sln->set_active_element(mesh->elements[son_idx]);

		Trf *tr = get_trf(int_trf[split][e]);
		for (int i = 0; i < bubble_fns; i++) {
			int iidx = bubble_fn_idx[i];
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

			scalar g[np], dgdx[np], dgdy[np], dgdz[np];
			memset(g, 0, np * sizeof(scalar));
			memset(dgdx, 0, np * sizeof(scalar));
			memset(dgdy, 0, np * sizeof(scalar));
			memset(dgdz, 0, np * sizeof(scalar));

			for (int l = 0; l < ipol_fns; l++) {
				double h[np];
				scalar sch[np];
				ss->get_fn_values(ipol[l].idx, np, tpt, 0, h);
				for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
				blas_axpy(np, ipol[l].coef, sch, 1, g, 1);
				ss->get_dx_values(ipol[l].idx, np, tpt, 0, h);
				for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
				blas_axpy(np, ipol[l].coef, sch, 1, dgdx, 1);
				ss->get_dy_values(ipol[l].idx, np, tpt, 0, h);
				for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
				blas_axpy(np, ipol[l].coef, sch, 1, dgdy, 1);
				ss->get_dz_values(ipol[l].idx, np, tpt, 0, h);
				for (int ii = 0; ii < np; ii++) sch[ii] = h[ii];
				blas_axpy(np, ipol[l].coef, sch, 1, dgdz, 1);
			}

			scalar value = 0.0;
			for (int k = 0; k < quad->get_num_points(order_rhs); k++) {
				value += pt[k].w * (uval[k] * (rval[k] - g[k]) +
					dudx[k] * ((drdx[k] * mdx[split]) - dgdx[k]) +
					dudy[k] * ((drdy[k] * mdy[split]) - dgdy[k]) +
					dudz[k] * ((drdz[k] * mdz[split]) - dgdz[k]));
			}
			proj_rhs[i] += value * (1 / (double) int_ns[split]);

			if (int_trf[split][e] != -1) fu->pop_transform();
		}
	}

	double d;
	int iperm[bubble_fns];
	ludcmp(proj_mat, bubble_fns, iperm, &d);
	lubksb(proj_mat, bubble_fns, iperm, proj_rhs);

	bubble_proj = new ProjItem [bubble_fns];
	for (int i = 0; i < bubble_fns; i++) {
		bubble_proj[i].coef = proj_rhs[i];
		bubble_proj[i].idx = bubble_fn_idx[i];
	}

	delete [] proj_mat;
	delete [] proj_rhs;
}
