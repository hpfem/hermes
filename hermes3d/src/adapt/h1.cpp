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
#include "../common.h"
#include "../solution.h"
#include "../refmap.h"
#include "../quad.h"
#include "../matrix.h"
#include "../traverse.h"
#include "../norm.h"
#include "../forms.h"
#include "h1.h"
#include "h1projipol.h"
#include <common/timer.h>
#include <common/callstack.h>

//#define DEBUG_PRINT

template<typename f_t, typename res_t>
res_t h1_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<res_t> *u, fn_t<res_t> *v, geom_t<f_t> *e,
              user_data_t<res_t> *ext)
{
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (
				u->fn[i] * CONJ(v->fn[i]) +
				u->dx[i] * CONJ(v->dx[i]) +
				u->dy[i] * CONJ(v->dy[i]) +
				u->dz[i] * CONJ(v->dz[i]));
	return result;
}

// H1 adapt ///////////////////////////////////////////////////////////////////////////////////////

void H1Adapt::init(Tuple<Space *> sp)
{
	_F_
	this->num = sp.size();

	this->spaces = new Space *[this->num];
	for (int i = 0; i < this->num; i++) spaces[i] = sp[i];

	this->sln = new Solution *[this->num];
	this->rsln = new Solution *[this->num];
	this->errors = new double *[this->num];
	this->norms = new double [this->num];
	for (int i = 0; i < this->num; i++) {
		this->sln[i] = NULL;
		this->rsln[i] = NULL;
		this->errors[i] = NULL;
		this->norms[i] = 0.0;
	}

	form = new_matrix<biform_val_t>(num, num);
	ord = new_matrix<biform_ord_t>(num, num);
	for (int i = 0; i < num; i++)
		for (int j = 0; j < num; j++) {
			if (i == j) {
				form[i][j] = h1_form<double, scalar>;
				ord[i][j]  = h1_form<ord_t, ord_t>;
			}
			else {
				form[i][j] = NULL;
				ord[i][j]  = NULL;
			}
		}

	esort = NULL;
	have_errors = false;

	// default parameters
	h_only = false;
	strategy = 0;
	max_order = H3D_MAX_ELEMENT_ORDER;
	aniso = false;			// FIXME: when implementing aniso, change this to true
	exponent = 1.0 / 3.0;

	log_file = NULL;
}

H1Adapt::~H1Adapt()
{
	_F_
	for (int i = 0; i < num; i++)
		delete [] errors[i];
	delete [] sln;
	delete [] rsln;
	delete [] errors;
	delete [] norms;
	delete [] spaces;
	delete [] form;
	delete [] ord;

	delete [] esort;
}

void H1Adapt::set_error_form(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord)
{
	if (i < 0 || i >= num || j < 0 || j >= num)
		error("Invalid equation number.");

	form[i][j] = bi_form;
	ord[i][j] = bi_ord;
}

double H1Adapt::get_projection_error(Element *e, int split, int son, const order3_t &order, Solution *rsln,
                                     Shapeset *ss)
{
	_F_
	ProjKey key(split, son, order);
	double err;
	if (proj_err.lookup(key, err))
		return err;
	else {
		H1ProjectionIpol proj(rsln, e, ss);
		err = proj.get_error(split, son, order);
		proj_err.set(key, err);
		return err;
	}


}

//// optimal refinement search /////////////////////////////////////////////////////////////////////

static inline int ndofs_elem(const order3_t &order)
{
	assert(order.type == MODE_HEXAHEDRON);
	return (order.x + 1) * (order.y + 1) * (order.z + 1);
}

static inline int ndofs_bubble(const order3_t &order)
{
	return (order.x - 1) * (order.y - 1) * (order.z - 1);
}

static inline int ndofs_face(int face, const order3_t &order1, const order3_t &order2)
{
	order2_t forder[] = { order1.get_face_order(face), order2.get_face_order(face) };
	return (
		(std::min(forder[0].x, forder[1].x) - 1) *
		(std::min(forder[0].y, forder[1].y) - 1));
}

static inline int ndofs_edge(int edge, const order3_t &o)
{
	return o.get_edge_order(edge) - 1;
}

static inline int ndofs_edge(int edge, const order3_t &o1, const order3_t &o2)
{
	return std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)) - 1;
}

static inline int ndofs_edge(int edge, const order3_t &o1, const order3_t &o2, const order3_t &o3,
                             const order3_t &o4)
{
	return
		std::min(
			std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)),
			std::min(o3.get_edge_order(edge), o4.get_edge_order(edge))
		) - 1;
}

int H1Adapt::get_dof_count(int split, order3_t order[])
{
	_F_
	int dofs = 0;
	switch (split) {
		case H3D_REFT_HEX_NONE:
			dofs = ndofs_elem(order[0]);
			break;

		case H3D_REFT_HEX_X:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(0, order[0], order[1]);			// face
			dofs -= 2 * ndofs_edge(3, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(4, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case H3D_REFT_HEX_Y:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(2, order[0], order[1]);
			dofs -= 2 * ndofs_edge(0, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(5, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case H3D_REFT_HEX_Z:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(4, order[0], order[1]);
			dofs -= 2 * ndofs_edge(0, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(1, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case H3D_H3D_REFT_HEX_XY:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2])
				+ ndofs_elem(order[3]);
			dofs -= ndofs_face(1, order[0], order[1]);			// faces
			dofs -= ndofs_face(3, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]);
			dofs -= ndofs_face(2, order[3], order[0]);
			dofs -= 5 * ndofs_edge(4, order[0]);				// edge
			dofs -= 2 * ndofs_edge(0, order[0], order[3]) + 2 * ndofs_edge(0, order[1], order[2]);
			dofs -= 2 * ndofs_edge(1, order[0], order[1]) + 2 * ndofs_edge(1, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case H3D_H3D_REFT_HEX_XZ:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2])
				+ ndofs_elem(order[3]);
			dofs -= ndofs_face(1, order[0], order[1]);			// faces
			dofs -= ndofs_face(5, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]);
			dofs -= ndofs_face(4, order[3], order[0]);
			dofs -= 5 * ndofs_edge(1, order[0]);				// edge
			dofs -= 2 * ndofs_edge(0, order[0], order[3]) + 2 * ndofs_edge(0, order[1], order[2]);
			dofs -= 2 * ndofs_edge(4, order[0], order[1]) + 2 * ndofs_edge(4, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case H3D_H3D_REFT_HEX_YZ:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2])
				+ ndofs_elem(order[3]);
			dofs -= ndofs_face(3, order[0], order[1]);			// faces
			dofs -= ndofs_face(5, order[1], order[2]);
			dofs -= ndofs_face(2, order[2], order[3]);
			dofs -= ndofs_face(4, order[3], order[0]);
			dofs -= 5 * ndofs_edge(0, order[0]);				// edge
			dofs -= 2 * ndofs_edge(1, order[0], order[3]) + 2 * ndofs_edge(1, order[1], order[2]);
			dofs -= 2 * ndofs_edge(4, order[0], order[1]) + 2 * ndofs_edge(4, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case H3D_H3D_H3D_REFT_HEX_XYZ:
			for (int i = 0; i < 8; i++) dofs += ndofs_elem(order[i]);
			dofs -= 15;			// vertex fns
			dofs -= ndofs_edge(0, order[0], order[4]) + ndofs_edge(0, order[3], order[7])
				+ ndofs_edge(0, order[0], order[3], order[4], order[7]);
			dofs -= ndofs_edge(0, order[1], order[5]) + ndofs_edge(0, order[2], order[6])
				+ ndofs_edge(0, order[1], order[2], order[5], order[6]);

			dofs -= ndofs_edge(1, order[0], order[4]) + ndofs_edge(1, order[1], order[5])
				+ ndofs_edge(1, order[0], order[1], order[4], order[5]);
			dofs -= ndofs_edge(1, order[3], order[7]) + ndofs_edge(1, order[2], order[6])
				+ ndofs_edge(1, order[2], order[3], order[6], order[7]);

			dofs -= ndofs_edge(4, order[0], order[1]) + ndofs_edge(4, order[2], order[3])
				+ ndofs_edge(4, order[0], order[1], order[2], order[3]);
			dofs -= ndofs_edge(4, order[4], order[5]) + ndofs_edge(4, order[6], order[7])
				+ ndofs_edge(4, order[4], order[5], order[6], order[7]);

			dofs -= ndofs_edge(4, order[0], order[3]) + ndofs_edge(4, order[4], order[7]);
			dofs -= ndofs_edge(4, order[1], order[2]) + ndofs_edge(4, order[5], order[6]);

			dofs -= ndofs_edge(1, order[0], order[1]) + ndofs_edge(1, order[2], order[3]);
			dofs -= ndofs_edge(1, order[4], order[5]) + ndofs_edge(1, order[6], order[7]);

			dofs -= ndofs_edge(0, order[1], order[2]) + ndofs_edge(0, order[0], order[3]);
			dofs -= ndofs_edge(0, order[5], order[6]) + ndofs_edge(0, order[4], order[7]);


			dofs -= ndofs_face(1, order[0], order[1]) + ndofs_face(3, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]) + ndofs_face(2, order[0], order[3]);

			dofs -= ndofs_face(5, order[0], order[4]) + ndofs_face(5, order[1], order[5]);
			dofs -= ndofs_face(5, order[2], order[6]) + ndofs_face(5, order[3], order[7]);

			dofs -= ndofs_face(1, order[4], order[5]) + ndofs_face(3, order[5], order[6]);
			dofs -= ndofs_face(0, order[6], order[7]) + ndofs_face(2, order[7], order[4]);
			break;

		default: assert(false);
	}

	return dofs;
}

void H1Adapt::get_optimal_refinement(Mesh *mesh, Element *e, const order3_t &order, Solution *rsln,
                                     Shapeset *ss, int &split, order3_t p[8])
{
	_F_
	int i, k, n = 0;
	const int MAX_CAND = 5000;

	struct Cand {
		double error;
		int dofs, split;
		order3_t p[Hex::NUM_SONS];				// polynomial degree

		Cand() {
			error = 0.0;
			dofs = 0;
			split = -1;
			memset(p, 0, sizeof(p));
		}
	};
	Cand cand[MAX_CAND];

#define MAKE_P_CAND(q) { \
    assert(n < MAX_CAND);   \
    cand[n].split = H3D_REFT_HEX_NONE; \
    cand[n].p[1] = cand[n].p[2] = cand[n].p[3] = cand[n].p[4] =\
    cand[n].p[5] = cand[n].p[6] = cand[n].p[7] = 0; \
    cand[n].p[0] = (q); \
    n++; }

#define MAKE_HP_CAND(q0, q1, q2, q3, q4, q5, q6, q7) { \
    assert(n < MAX_CAND);  \
    cand[n].split = H3D_H3D_H3D_REFT_HEX_XYZ; \
    cand[n].p[0] = (q0); \
    cand[n].p[1] = (q1); \
    cand[n].p[2] = (q2); \
    cand[n].p[3] = (q3); \
    cand[n].p[4] = (q4); \
    cand[n].p[5] = (q5); \
    cand[n].p[6] = (q6); \
    cand[n].p[7] = (q7); \
    n++; }

#define MAKE_ANI2_CAND(s, q0, q1) { \
	if (mesh->can_refine_element(e->id, s)) {\
		assert(n < MAX_CAND);  \
		cand[n].split = s; \
		cand[n].p[2] = cand[n].p[3] = cand[n].p[4] = cand[n].p[5] = cand[n].p[6] =\
		cand[n].p[7] = 0; \
		cand[n].p[0] = (q0); \
		cand[n].p[1] = (q1); \
		n++; }}

#define MAKE_ANI4_CAND(s, q0, q1, q2, q3) { \
	if (mesh->can_refine_element(e->id, s)) {\
		assert(n < MAX_CAND);  \
		cand[n].split = s; \
		cand[n].p[4] = cand[n].p[5] = cand[n].p[6] = cand[n].p[7] = 0; \
		cand[n].p[0] = (q0); \
		cand[n].p[1] = (q1); \
		cand[n].p[2] = (q2); \
		cand[n].p[3] = (q3); \
		n++; }}

	order3_t pp[] = {
		order3_t(order.x, order.y, order.z),
		order3_t(std::min(max_order, order.x + 1),
				 std::min(max_order, order.y + 1),
				 std::min(max_order, order.z + 1)),
		order3_t(std::min(max_order, order.x + 2),
				 std::min(max_order, order.y + 2),
				 std::min(max_order, order.z + 2))
	};

	if (h_only) {
		MAKE_P_CAND(pp[0]);
		MAKE_HP_CAND(pp[0], pp[0], pp[0], pp[0], pp[0], pp[0], pp[0], pp[0]);

		if (aniso) {
			MAKE_ANI2_CAND(H3D_REFT_HEX_X, pp[0], pp[0]);
			MAKE_ANI2_CAND(H3D_REFT_HEX_Y, pp[0], pp[0]);
			MAKE_ANI2_CAND(H3D_REFT_HEX_Z, pp[0], pp[0]);

			MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_XY, pp[0], pp[0], pp[0], pp[0]);
			MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_YZ, pp[0], pp[0], pp[0], pp[0]);
			MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_XZ, pp[0], pp[0], pp[0], pp[0]);
		}
	}
	else {
		// prepare p-candidates
		if (aniso) {
			for (unsigned int q0 = 0; q0 < countof(pp); q0++)
				for (unsigned int q1 = 0; q1 < countof(pp); q1++)
					for (unsigned int q2 = 0; q2 < countof(pp); q2++)
						MAKE_P_CAND(order3_t(pp[q0].x, pp[q1].y, pp[q2].z));
		}
		else {
			MAKE_P_CAND(pp[0]);
			MAKE_P_CAND(pp[1]);
			MAKE_P_CAND(pp[2]);
		}

		// prepare hp-candidates
		{
			order3_t hpp[] = {
				order3_t((order.x + 1) / 2, (order.y + 1) / 2, (order.z + 1) / 2),
				order3_t(std::min(((order.x + 1) / 2) + 1, max_order),
				         std::min(((order.y + 1) / 2) + 1, max_order),
				         std::min(((order.z + 1) / 2) + 1, max_order))
			};

			for (unsigned int q0 = 0; q0 < countof(hpp); q0++)
				for (unsigned int q1 = 0; q1 < countof(hpp); q1++)
					for (unsigned int q2 = 0; q2 < countof(hpp); q2++)
						for (unsigned int q3 = 0; q3 < countof(hpp); q3++)
							for (unsigned int q4 = 0; q4 < countof(hpp); q4++)
								for (unsigned int q5 = 0; q5 < countof(hpp); q5++)
									for (unsigned int q6 = 0; q6 < countof(hpp); q6++)
										for (unsigned int q7 = 0; q7 < countof(hpp); q7++)
											MAKE_HP_CAND(hpp[q0], hpp[q1], hpp[q2], hpp[q3],
											             hpp[q4], hpp[q5], hpp[q6], hpp[q7]);

		}

		if (aniso) {
			// X
			order3_t ppx[] = {
				order3_t((order.x + 1) / 2, order.y, order.z),
				order3_t(std::min(((order.x + 1) / 2) + 1, max_order), order.y, order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppx); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppx); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_X, ppx[q0], ppx[q1]);
			// Y
			order3_t ppy[] = {
				order3_t(order.x, (order.y + 1) / 2, order.z),
				order3_t(order.x, std::min(((order.y + 1) / 2) + 1, max_order), order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppy); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppy); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_Y, ppy[q0], ppy[q1]);
			// Z
			order3_t ppz[] = {
				order3_t(order.x, order.y, (order.z + 1) / 2),
				order3_t(order.x, order.y, std::min(((order.z + 1) / 2) + 1, max_order)),
			};
			for (unsigned int q0 = 0; q0 < countof(ppz); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppz); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_Z, ppz[q0], ppz[q1]);

			// XY
			order3_t ppxy[] = {
				order3_t((order.x + 1) / 2, (order.y + 1) / 2, order.z),
				order3_t(std::min(((order.x + 1) / 2) + 1, max_order),
				         std::min(((order.y + 1) / 2) + 1, max_order), order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppxy); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppxy); q1++)
					for (unsigned int q2 = 0; q2 < countof(ppxy); q2++)
						for (unsigned int q3 = 0; q3 < countof(ppxy); q3++)
							MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_XY, ppxy[q0], ppxy[q1], ppxy[q2], ppxy[q3]);
			// YZ
			order3_t ppyz[] = {
					order3_t(order.x, (order.y + 1) / 2, (order.z + 1) / 2),
					order3_t(order.x, std::min(((order.y + 1) / 2) + 1, max_order),
					         std::min(((order.z + 1) / 2) + 1, max_order)),
			};
			for (unsigned int q0 = 0; q0 < countof(ppyz); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppyz); q1++)
					for (unsigned int q2 = 0; q2 < countof(ppyz); q2++)
						for (unsigned int q3 = 0; q3 < countof(ppyz); q3++)
							MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_YZ, ppyz[q0], ppyz[q1], ppyz[q2], ppyz[q3]);
			// XZ
			order3_t ppxz[] = {
				order3_t((order.x + 1) / 2, order.y, (order.z + 1) / 2),
				order3_t(std::min(((order.x + 1) / 2) + 1, max_order), order.y,
				         std::min(((order.z + 1) / 2) + 1, max_order)),
			};
			for (unsigned int q0 = 0; q0 < countof(ppxz); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppxz); q1++)
					for (unsigned int q2 = 0; q2 < countof(ppxz); q2++)
						for (unsigned int q3 = 0; q3 < countof(ppxz); q3++)
							MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_XZ, ppxz[q0], ppxz[q1], ppxz[q2], ppxz[q3]);
		}
	}

#ifdef DEBUG_PRINT
	const char *split_str[] = {
		"NONE", "X   ", "Y   ", "Z   ", "XY  ", "XZ  ", "YZ  ", "XYZ "
	};
#endif

	// calculate their errors
	for (i = k = 0; i < n; i++) {
		Cand *c = cand + i;

		c->error = 0.0;
		switch (c->split) {
			case H3D_REFT_HEX_NONE:
				c->error += get_projection_error(e, c->split, -1, c->p[0], rsln, ss);
				break;

			case H3D_H3D_H3D_REFT_HEX_XYZ:
				for (int j = 0; j < 8; j++)
					c->error += get_projection_error(e, c->split, j, c->p[j], rsln, ss);
				break;

			case H3D_REFT_HEX_X:
				c->error += get_projection_error(e, c->split, 20, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 21, c->p[1], rsln, ss);
				break;

			case H3D_REFT_HEX_Y:
				c->error += get_projection_error(e, c->split, 22, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 23, c->p[1], rsln, ss);
				break;

			case H3D_REFT_HEX_Z:
				c->error += get_projection_error(e, c->split, 24, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 25, c->p[1], rsln, ss);
				break;

			case H3D_H3D_REFT_HEX_XY:
				c->error += get_projection_error(e, c->split,  8, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split,  9, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 10, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 11, c->p[3], rsln, ss);
				break;

			case H3D_H3D_REFT_HEX_XZ:
				c->error += get_projection_error(e, c->split, 12, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 13, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 14, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 15, c->p[3], rsln, ss);
				break;

			case H3D_H3D_REFT_HEX_YZ:
				c->error += get_projection_error(e, c->split, 16, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 17, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 18, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 19, c->p[3], rsln, ss);
				break;

			default:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;
		}
		c->error = sqrt(c->error);
		c->dofs = get_dof_count(c->split, c->p);
	}

#ifdef DEBUG_PRINT
	printf("\n");
	printf("- cand: #0: dofs = %d", cand[0].dofs);
	printf(" | order = (%d, %d, %d)", order.x, order.y, order.z);
	printf(" | err = % .15e", log10(cand[0].error));
	printf("\n");
#endif

	// select an above-average candidate with the steepest error decrease
	int imax = 0;
	double score, maxscore = 0.0;
	for (i = 1; i < n; i++) {
		if (cand[i].dofs - cand[0].dofs > 0) {
			score = (log10(cand[0].error) - log10(cand[i].error)) / pow(cand[i].dofs - cand[0].dofs, exponent);

#ifdef DEBUG_PRINT
			printf("- cand: #%d, split = %s", i, split_str[cand[i].split]);
			for (int ii = 0; ii < 8; ii++)
				printf(", (%d, %d, %d)", cand[i].p[ii].x, cand[i].p[ii].y, cand[i].p[ii].z);
			printf(" | dofs = %d", cand[i].dofs);
			printf(" | err = % .15e", log10(cand[i].error));
			printf(" | derr = % e", log10(cand[0].error) - log10(cand[i].error));
			printf(" | score = % e ", score);
			printf("\n");
#endif

			if (score > maxscore) {
				maxscore = score;
				imax = i;
			}
		}
	}

	// return result
	split = cand[imax].split;
	memcpy(p, cand[imax].p, Hex::NUM_SONS * sizeof(order3_t));

#ifdef DEBUG_PRINT
	printf(": best cand: #%d, split = %s", imax, split_str[cand[imax].split]);
	for (int i = 0; i < 8; i++)
		printf(", (%d, %d, %d)", p[i].x, p[i].y, p[i].z);
	printf("\n");
#endif
}

//// adapt /////////////////////////////////////////////////////////////////////////////////////////

void H1Adapt::adapt(double thr)
{
	_F_
	if (!have_errors)
		EXIT("Element errors have to be calculated first, see calc_error().");

	Timer tmr;
	tmr.start();

	Mesh *mesh[num];
	for (int j = 0; j < num; j++) {
		mesh[j] = spaces[j]->get_mesh();
		rsln[j]->enable_transform(false);
	}

	if (log_file != NULL) fprintf(log_file, "--\n");

	double err0 = 1000.0;
	double processed_error = 0.0;
	int i = 0;
	for (i = 0; i < nact; i++) {
		int comp = esort[i][1];
		int id = esort[i][0];
		double err = errors[comp][id - 1];

		// first refinement strategy:
		// refine elements until prescribed amount of error is processed
		// if more elements have similar error refine all to keep the mesh symmetric
		if ((strategy == 0) && (processed_error > sqrt(thr) * total_err) &&
				fabs((err - err0) / err0) > 1e-3)
			break;

		// second refinement strategy:
		// refine all elements whose error is bigger than some portion of maximal error
		if ((strategy == 1) && (err < thr * errors[esort[0][1]][esort[0][0] - 1]))
			break;

		assert(mesh[comp]->elements.exists(id));
		Element *e = mesh[comp]->elements[id];
#ifdef DEBUG_PRINT
		printf("  - element #%d", id);
#endif

		int split = 0;
		order3_t p[8];													// polynomial order of sons
		for (int k = 0; k < 8; k++) p[k] = order3_t(0, 0, 0);
		order3_t cur_order = spaces[comp]->get_element_order(id);

		if (h_only && !aniso) {
			p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = p[6] = p[7] = cur_order;
			split = H3D_H3D_H3D_REFT_HEX_XYZ;
#ifdef DEBUG_PRINT
			printf("\n");			// new-line
#endif
		}
		else
			get_optimal_refinement(mesh[comp], e, cur_order, rsln[comp],
			                       spaces[comp]->get_shapeset(), split, p);

		if (log_file != NULL)
			fprintf(log_file, "%ld %d %d %d %d %d %d %d %d %d\n", e->id, split,
				p[0].get_idx(), p[1].get_idx(), p[2].get_idx(), p[3].get_idx(),
				p[4].get_idx(), p[5].get_idx(), p[6].get_idx(), p[7].get_idx());

		switch (split) {
			case H3D_REFT_HEX_NONE:
				spaces[comp]->set_element_order(id, p[0]);
				break;

			case H3D_H3D_H3D_REFT_HEX_XYZ:
				mesh[comp]->refine_element(id, H3D_H3D_H3D_REFT_HEX_XYZ);
				for (int j = 0; j < Hex::NUM_SONS; j++)
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			case H3D_REFT_HEX_X:
			case H3D_REFT_HEX_Y:
			case H3D_REFT_HEX_Z:
				mesh[comp]->refine_element(id, split);
				for (int j = 0; j < 2; j++)
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			case H3D_H3D_REFT_HEX_XY:
			case H3D_H3D_REFT_HEX_XZ:
			case H3D_H3D_REFT_HEX_YZ:
				mesh[comp]->refine_element(id, split);
				for (int j = 0; j < 4; j++)
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			default: assert(false);
		}

		err0 = err;
		processed_error += err;

		proj_err.remove_all();
	}

	for (int j = 0; j < num; j++)
		rsln[j]->enable_transform(true);

	have_errors = false;

	reft_elems = i;

	tmr.stop();
	adapt_time = tmr.get_seconds();
}

//// Unrefinements /////////////////////////////////////////////////////////////////////////////////

// TODO: unrefts

//// error calculation /////////////////////////////////////////////////////////////////////////////

static double **cmp_err;
static int compare(const void* p1, const void* p2)
{
	const int2 (*e1) = ((const int2 *) p1);
	const int2 (*e2) = ((const int2 *) p2);
	return cmp_err[(*e1)[1]][(*e1)[0] - 1] < cmp_err[(*e2)[1]][(*e2)[0] - 1] ? 1 : -1;
}

order3_t H1Adapt::get_form_order(int marker, const order3_t &ordu, const order3_t &ordv, RefMap *ru,
                                 matrix_form_ord_t mf_ord)
{
	_F_
	// determine the integration order
	fn_t<ord_t> ou = init_fn(ordu);
	fn_t<ord_t> ov = init_fn(ordv);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(marker);
	ord_t o = mf_ord(1, &fake_wt, NULL, &ou, &ov, &fake_e, NULL);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();

	free_fn(&ou);
	free_fn(&ov);

	return order;
}

scalar H1Adapt::eval_error(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *sln1,
                           MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2)
{
	_F_
	RefMap *rv1 = sln1->get_refmap();
	RefMap *rv2 = sln1->get_refmap();
	RefMap *rrv1 = rsln1->get_refmap();
	RefMap *rrv2 = rsln1->get_refmap();

	order3_t order = get_form_order(marker, rsln1->get_fn_order(), rsln2->get_fn_order(), rrv1,
	                                bi_ord);

	// eval the form
	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);		// FIXME: hex-specific
	QuadPt3D *pt = quad->get_points(order);
	int np = quad->get_num_points(order);

	double *jwt = rrv1->get_jacobian(np, pt);
	geom_t<double> e = init_geom(marker, rrv1, np, pt);

	fn_t<scalar> *err1 = init_fn(sln1, rv1, np, pt);
	fn_t<scalar> *err2 = init_fn(sln2, rv2, np, pt);
	fn_t<scalar> *v1 = init_fn(rsln1, rrv1, np, pt);
	fn_t<scalar> *v2 = init_fn(rsln2, rrv2, np, pt);

	for (int i = 0; i < np; i++) {
		err1->fn[i] = err1->fn[i] - v1->fn[i];
		err1->dx[i] = err1->dx[i] - v1->dx[i];
		err1->dy[i] = err1->dy[i] - v1->dy[i];
		err1->dz[i] = err1->dz[i] - v1->dz[i];
		err2->fn[i] = err2->fn[i] - v2->fn[i];
		err2->dx[i] = err2->dx[i] - v2->dx[i];
		err2->dy[i] = err2->dy[i] - v2->dy[i];
		err2->dz[i] = err2->dz[i] - v2->dz[i];
	}

	scalar res = bi_fn(np, jwt, NULL, err1, err2, &e, NULL);

	delete [] jwt;
	free_geom(&e);
	free_fn(err1);
	free_fn(err2);
	free_fn(v1);
	free_fn(v2);

	return res;
}


scalar H1Adapt::eval_norm(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1,
                          MeshFunction *rsln2)
{
	_F_
	RefMap *rv1 = rsln1->get_refmap();
	RefMap *rv2 = rsln1->get_refmap();

	order3_t order = get_form_order(marker, rsln1->get_fn_order(), rsln2->get_fn_order(), rv1,
	                                bi_ord);

	// eval the form
	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON); 	// FIXME: hex_specific
	QuadPt3D *pt = quad->get_points(order);
	int np = quad->get_num_points(order);

	double *jwt = rv1->get_jacobian(np, pt);
	geom_t<double> e = init_geom(marker, rv1, np, pt);

	fn_t<scalar> *v1 = init_fn(rsln1, rv1, np, pt);
	fn_t<scalar> *v2 = init_fn(rsln2, rv2, np, pt);

	scalar res = bi_fn(np, jwt, NULL, v1, v2, &e, NULL);

	delete [] jwt;
	free_geom(&e);
	free_fn(v1);
	free_fn(v2);

	return res;
}

double H1Adapt::calc_error_n(Tuple<Solution *> slns, Tuple<Solution *> rslns)
{
	_F_
	int i, j, k;

	int n = slns.size();
	if (n != this->num) EXIT("Wrong number of solutions.");

	Timer tmr;
	tmr.start();

	for (i = 0; i < n; i++) {
	  this->sln[i] = slns[i];
	  this->sln[i]->enable_transform(true);
	}
	for (i = 0; i < n; i++) {
	  this->rsln[i] = rslns[i];
	  this->rsln[i]->enable_transform(true);
	}

	// prepare multi-mesh traversal and error arrays
	Mesh *meshes[2 * num];
	Transformable *tr[2 * num];
	Traverse trav;
	nact = 0;
	for (i = 0; i < num; i++) {
		meshes[i] = sln[i]->get_mesh();
		meshes[i + num] = rsln[i]->get_mesh();
		tr[i] = sln[i];
		tr[i + num] = rsln[i];

		nact += sln[i]->get_mesh()->get_num_active_elements();

		int max = meshes[i]->get_max_element_id();
		if (errors[i] != NULL) delete [] errors[i];
		errors[i] = new double[max];
		memset(errors[i], 0, sizeof(double) * max);
	}

	double total_norm = 0.0;
	double norms[num];
	memset(norms, 0, num * sizeof(double));
	double total_error = 0.0;
	if (esort != NULL) delete [] esort;
	esort = new int2[nact];

	Element **ee;
	trav.begin(2 * num, meshes, tr);
	while ((ee = trav.get_next_state(NULL, NULL)) != NULL) {
		Element *e0 = NULL;
		for (i = 0; i < 2 * num; i++)
			if ((e0 = ee[i]) != NULL) break;
		assert(e0 != NULL);

		for (i = 0; i < num; i++) {
			for (j = 0; j < num; j++) {
				if (form[i][j] != NULL) {
					double err, nrm;
#ifndef H3D_COMPLEX
					err = fabs(eval_error(e0->marker, form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j]));
					nrm = fabs(eval_norm(e0->marker, form[i][j], ord[i][j], rsln[i], rsln[j]));
#else
					err = std::abs(eval_error(e0->marker, form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j]));
					nrm = std::abs(eval_norm(e0->marker, form[i][j], ord[i][j], rsln[i], rsln[j]));
#endif
					norms[i] += nrm;
					total_norm  += nrm;
					total_error += err;
					errors[i][ee[i]->id - 1] += err;
				}

			}
		}
	}
	trav.finish();

	k = 0;
	for (i = 0; i < num; i++)
		FOR_ALL_ACTIVE_ELEMENTS(eid, meshes[i]) {
			Element *e = meshes[i]->elements[eid];
			esort[k][0] = e->id;
			esort[k++][1] = i;
			// FIXME: when norms of 2 components are very different it can help (microwave heating)
			// navier-stokes on different meshes work only without it
			errors[i][e->id - 1] /= norms[i];
		}

	assert(k == nact);
	cmp_err = errors;
	qsort(esort, nact, sizeof(int2), compare);

	have_errors = true;

	total_err = total_error / total_norm;		// FIXME: comment out the denominator when
												// commenting out the above fixme
	tmr.stop();
	error_time = tmr.get_seconds();

#ifdef DEBUG_PRINT
	printf("\n");
	for (int i = 0; i < std::min(10, nact); i++) {
		printf("  - element error #%d = % e\n", esort[i][0], errors[0][esort[i][0] - 1]);
	}
#endif

	return sqrt(total_error / total_norm);
}
