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

#include "../h3d_common.h"
#include "../solution.h"
#include "../refmap.h"
#include "../quad.h"
#include "../traverse.h"
#include "../norm.h"
#include "../weakform/forms.h"
#include "adapt.h"
#include "h1projipol.h"
#include "h1proj.h"
#include "hcurlproj.h"
#include "../../../hermes_common/matrix.h"

//#define DEBUG_PRINT

template<typename f_t, typename res_t>
res_t h1_form(int n, double *wt, Func<res_t> *u_ext[], Func<res_t> *u, Func<res_t> *v, Geom<f_t> *e,
              ExtData<res_t> *ext)
{
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (
				u->val[i] * conj(v->val[i]) +
				u->dx[i] * conj(v->dx[i]) +
				u->dy[i] * conj(v->dy[i]) +
				u->dz[i] * conj(v->dz[i]));
	return result;
}

template<typename f_t, typename res_t>
res_t hcurl_form(int n, double *wt, Func<res_t> *u_ext[], Func<res_t> *u, Func<res_t> *v, Geom<f_t> *e,
              ExtData<res_t> *ext)
{
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->curl0[i] * conj(v->curl0[i]) +
                       u->curl1[i] * conj(v->curl1[i]) +
                       u->curl2[i] * conj(v->curl2[i]) +
                       u->val0[i] * conj(v->val0[i]) +
                       u->val1[i] * conj(v->val1[i]) +
                       u->val2[i] * conj(v->val2[i]));
	return result;
}

// H1 adapt ///////////////////////////////////////////////////////////////////////////////////////

void Adapt::init(Hermes::vector<Space *> sp, Hermes::vector<ProjNormType> proj_norms)
{
	_F_
	this->num = sp.size();
  this->proj_norms = proj_norms;

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
			if (i == j && proj_norms.size() > 0) {
				switch (proj_norms[i]) {
        case HERMES_H1_NORM: form[i][i] = h1_form<double, scalar>; ord[i][i]  = h1_form<Ord, Ord>; break;
        case HERMES_HCURL_NORM: form[i][i] = hcurl_form<double, scalar>; ord[i][i]  = hcurl_form<Ord, Ord>; break;
        default: error("Unknown projection type in Adapt::Adapt().");
        }  
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
	aniso = true;
	exponent = 1.0 / 3.0;

	log_file = NULL;
}

Adapt::~Adapt()
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

void Adapt::set_error_form(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord)
{
	if (i < 0 || i >= num || j < 0 || j >= num)
		error("Invalid equation number.");

	form[i][j] = bi_form;
	ord[i][j] = bi_ord;
}

double Adapt::get_projection_error(Element *e, int split, int son, const Ord3 &order, Solution *rsln,
                                     Shapeset *ss)
{
	_F_
	ProjKey key(split, son, order);
	double err;
  Projection *proj;
  if (proj_err.find(key) != proj_err.end())
    return proj_err.find(key)->second;
	else {
    switch (ss->get_type()) {
      case HERMES_H1_SPACE:
        proj = new H1ProjectionIpol(rsln, e, ss);
		    err = proj->get_error(split, son, order);
		    proj_err[key] = err;
        delete proj;
		    return err;
        break;
      case HERMES_HCURL_SPACE:
        proj = new HCurlProjection(rsln, e, ss);
		    err = proj->get_error(split, son, order);
		    proj_err[key] = err;
        delete proj;
		    return err;
        break;
      default:
        error("Adaptivity only implemented for H1 and HCurl spaces.");
        return 0.0;
    }
  }
  return 0.0;
}

//// optimal refinement search /////////////////////////////////////////////////////////////////////

static inline int ndofs_elem(const Ord3 &order)
{
	assert(order.type == HERMES_MODE_HEX);
	return (order.x + 1) * (order.y + 1) * (order.z + 1);
}

static inline int ndofs_bubble(const Ord3 &order)
{
	return (order.x - 1) * (order.y - 1) * (order.z - 1);
}

static inline int ndofs_face(int face, const Ord3 &order1, const Ord3 &order2)
{
	Ord2 forder[] = { order1.get_face_order(face), order2.get_face_order(face) };
	return (
		(std::min(forder[0].x, forder[1].x) - 1) *
		(std::min(forder[0].y, forder[1].y) - 1));
}

static inline int ndofs_edge(int edge, const Ord3 &o)
{
	return o.get_edge_order(edge) - 1;
}

static inline int ndofs_edge(int edge, const Ord3 &o1, const Ord3 &o2)
{
	return std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)) - 1;
}

static inline int ndofs_edge(int edge, const Ord3 &o1, const Ord3 &o2, const Ord3 &o3,
                             const Ord3 &o4)
{
	return
		std::min(
			std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)),
			std::min(o3.get_edge_order(edge), o4.get_edge_order(edge))
		) - 1;
}

int Adapt::get_dof_count(int split, Ord3 order[])
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

void Adapt::get_optimal_refinement(Mesh *mesh, Element *e, const Ord3 &order, Solution *rsln,
                                     Shapeset *ss, int &split, Ord3 p[8])
{
	_F_
	int i, k, n = 0;
	const int MAX_CAND = 5000;

	struct Cand {
		double error;
		int dofs, split;
		Ord3 p[Hex::NUM_SONS];				// polynomial degree

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

	Ord3 pp[] = {
    Ord3(order.x, order.y, order.z),
    Ord3(std::min((int)max_order, (int)order.x + 1),
				 std::min((int)max_order, (int)order.y + 1),
				 std::min((int)max_order, (int)order.z + 1)),
    Ord3(std::min((int)max_order, (int)order.x + 2),
				 std::min((int)max_order, (int)order.y + 2),
				 std::min((int)max_order, (int)order.z + 2))
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
						MAKE_P_CAND(Ord3(pp[q0].x, pp[q1].y, pp[q2].z));
		}
		else {
			MAKE_P_CAND(pp[0]);
			MAKE_P_CAND(pp[1]);
			MAKE_P_CAND(pp[2]);
		}

		// prepare hp-candidates
		{
			Ord3 hpp[] = {
				Ord3((order.x + 1) / 2, (order.y + 1) / 2, (order.z + 1) / 2),
				Ord3(std::min((int)((order.x + 1) / 2) + 1, (int)max_order),
				         std::min((int)((order.y + 1) / 2) + 1, (int)max_order),
				         std::min((int)((order.z + 1) / 2) + 1, (int)max_order))
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
			Ord3 ppx[] = {
				Ord3((order.x + 1) / 2, order.y, order.z),
				Ord3(std::min((int)((order.x + 1) / 2) + 1, (int)max_order), order.y, order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppx); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppx); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_X, ppx[q0], ppx[q1]);
			// Y
			Ord3 ppy[] = {
				Ord3(order.x, (order.y + 1) / 2, order.z),
				Ord3(order.x, std::min((int)((order.y + 1) / 2) + 1, (int)max_order), order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppy); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppy); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_Y, ppy[q0], ppy[q1]);
			// Z
			Ord3 ppz[] = {
				Ord3(order.x, order.y, (order.z + 1) / 2),
				Ord3(order.x, order.y, std::min((int)((order.z + 1) / 2) + 1, (int)max_order)),
			};
			for (unsigned int q0 = 0; q0 < countof(ppz); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppz); q1++)
					MAKE_ANI2_CAND(H3D_REFT_HEX_Z, ppz[q0], ppz[q1]);

			// XY
			Ord3 ppxy[] = {
				Ord3((order.x + 1) / 2, (order.y + 1) / 2, order.z),
				Ord3(std::min((int)((order.x + 1) / 2) + 1, (int)max_order),
				         std::min((int)((order.y + 1) / 2) + 1, (int)max_order), order.z),
			};
			for (unsigned int q0 = 0; q0 < countof(ppxy); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppxy); q1++)
					for (unsigned int q2 = 0; q2 < countof(ppxy); q2++)
						for (unsigned int q3 = 0; q3 < countof(ppxy); q3++)
							MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_XY, ppxy[q0], ppxy[q1], ppxy[q2], ppxy[q3]);
			// YZ
			Ord3 ppyz[] = {
					Ord3(order.x, (order.y + 1) / 2, (order.z + 1) / 2),
					Ord3(order.x, std::min((int)((order.y + 1) / 2) + 1, (int)max_order),
					         std::min((int)((order.z + 1) / 2) + 1, (int)max_order)),
			};
			for (unsigned int q0 = 0; q0 < countof(ppyz); q0++)
				for (unsigned int q1 = 0; q1 < countof(ppyz); q1++)
					for (unsigned int q2 = 0; q2 < countof(ppyz); q2++)
						for (unsigned int q3 = 0; q3 < countof(ppyz); q3++)
							MAKE_ANI4_CAND(H3D_H3D_REFT_HEX_YZ, ppyz[q0], ppyz[q1], ppyz[q2], ppyz[q3]);
			// XZ
			Ord3 ppxz[] = {
				Ord3((order.x + 1) / 2, order.y, (order.z + 1) / 2),
				Ord3(std::min((int)((order.x + 1) / 2) + 1, (int)max_order), order.y,
				         std::min((int)((order.z + 1) / 2) + 1, (int)max_order)),
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
				EXIT(HERMES_ERR_NOT_IMPLEMENTED);
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
	memcpy(p, cand[imax].p, Hex::NUM_SONS * sizeof(Ord3));

#ifdef DEBUG_PRINT
	printf(": best cand: #%d, split = %s", imax, split_str[cand[imax].split]);
	for (int i = 0; i < 8; i++)
		printf(", (%d, %d, %d)", p[i].x, p[i].y, p[i].z);
	printf("\n");
#endif
}

//// adapt /////////////////////////////////////////////////////////////////////////////////////////

void Adapt::adapt(double thr)
{
	_F_
	if (!have_errors)
		EXIT("Element errors have to be calculated first, see calc_err_est().");

  TimePeriod tmr;

	Mesh ** mesh = new Mesh*[num];
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

		assert(mesh[comp]->elements[id] != NULL);
		Element *e = mesh[comp]->elements[id];
#ifdef DEBUG_PRINT
		printf("  - element #%d", id);
#endif

		int split = 0;
		Ord3 p[8];													// polynomial order of sons
		for (int k = 0; k < 8; k++) p[k] = Ord3(0, 0, 0);
		Ord3 cur_order = spaces[comp]->get_element_order(id);

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
			fprintf(log_file, "%u %d %d %d %d %d %d %d %d %d\n", e->id, split,
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

		proj_err.clear();
	}

	for (int j = 0; j < num; j++)
		rsln[j]->enable_transform(true);

	have_errors = false;

	reft_elems = i;

  tmr.tick();
  adapt_time = tmr.accumulated();

  for (int j = 0; j < num; j++)
    spaces[j]->assign_dofs();

  delete [] mesh;
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

Ord3 Adapt::get_form_order(int marker, const Ord3 &ordu, const Ord3 &ordv, RefMap *ru,
                                 matrix_form_ord_t mf_ord)
{
	_F_
	// determine the integration order
	Func<Ord> *ou = init_fn_ord(ordu);
	Func<Ord> *ov = init_fn_ord(ordv);

	double fake_wt = 1.0;
	Geom<Ord> fake_e = init_geom(marker);
	Ord o = mf_ord(1, &fake_wt, NULL, ou, ov, &fake_e, NULL);
	Ord3 order = ru->get_inv_ref_order();
	switch (order.type) {
		case HERMES_MODE_TET: order += Ord3(o.get_order()); break;
		case HERMES_MODE_HEX: order += Ord3(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();

	free_fn(ou);
	free_fn(ov);
  delete ou;
  delete ov;

	return order;
}

scalar Adapt::eval_error(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *sln1,
                           MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2)
{
	_F_
	RefMap *rv1 = sln1->get_refmap();
	RefMap *rv2 = sln1->get_refmap();
	RefMap *rrv1 = rsln1->get_refmap();
	RefMap *rrv2 = rsln1->get_refmap();

	Ord3 order = get_form_order(marker, rsln1->get_fn_order(), rsln2->get_fn_order(), rrv1,
	                                bi_ord);

	// eval the form
	Quad3D *quad = get_quadrature(sln1->get_active_element()->get_mode());
	QuadPt3D *pt = quad->get_points(order);
	int np = quad->get_num_points(order);

	double *jwt = rrv1->get_jacobian(np, pt);
	Geom<double> e = init_geom(marker, rrv1, np, pt);

	Func<scalar> *err1 = init_fn(sln1, rv1, np, pt);
	Func<scalar> *err2 = init_fn(sln2, rv2, np, pt);
	Func<scalar> *v1 = init_fn(rsln1, rrv1, np, pt);
	Func<scalar> *v2 = init_fn(rsln2, rrv2, np, pt);

	err1->subtract(*v1);
  err2->subtract(*v2);

	scalar res = bi_fn(np, jwt, NULL, err1, err2, &e, NULL);

	delete [] jwt;
	free_geom(&e);
	free_fn(err1);
	free_fn(err2);
	free_fn(v1);
	free_fn(v2);

	return res;
}


scalar Adapt::eval_norm(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1,
                          MeshFunction *rsln2)
{
	_F_
	RefMap *rv1 = rsln1->get_refmap();
	RefMap *rv2 = rsln1->get_refmap();

	Ord3 order = get_form_order(marker, rsln1->get_fn_order(), rsln2->get_fn_order(), rv1,
	                                bi_ord);

	// eval the form
	Quad3D *quad = get_quadrature(rsln1->get_active_element()->get_mode());
	QuadPt3D *pt = quad->get_points(order);
	int np = quad->get_num_points(order);

	double *jwt = rv1->get_jacobian(np, pt);
	Geom<double> e = init_geom(marker, rv1, np, pt);

	Func<scalar> *v1 = init_fn(rsln1, rv1, np, pt);
	Func<scalar> *v2 = init_fn(rsln2, rv2, np, pt);

	scalar res = bi_fn(np, jwt, NULL, v1, v2, &e, NULL);

	delete [] jwt;
	free_geom(&e);
	free_fn(v1);
	free_fn(v2);

	return res;
}

double Adapt::calc_err_internal(Hermes::vector<Solution *> slns, Hermes::vector<Solution *> rslns, unsigned int error_flags, Hermes::vector<double>* component_errors, bool solutions_for_adapt)
{
	_F_
	int i, j, k;

	int n = slns.size();
	if (n != this->num) EXIT("Wrong number of solutions.");

	TimePeriod tmr;

  Solution* rslns_original[10];
  Solution* slns_original[10];

	for (i = 0; i < n; i++) {
    slns_original[i] = this->sln[i];
	  this->sln[i] = slns[i];
	  this->sln[i]->enable_transform(true);
	}
	for (i = 0; i < n; i++) {
    rslns_original[i] = this->rsln[i];
	  this->rsln[i] = rslns[i];
	  this->rsln[i]->enable_transform(true);
	}

	// prepare multi-mesh traversal and error arrays
	Mesh **meshes = new Mesh *[2 * num];
	Transformable **tr = new Transformable *[2 * num];
	Traverse trav;
	nact = 0;
	for (i = 0; i < num; i++) {
		meshes[i] = sln[i]->get_mesh();
		meshes[i + num] = rsln[i]->get_mesh();
		tr[i] = sln[i];
		tr[i + num] = rsln[i];

		nact += sln[i]->get_mesh()->get_num_active_elements();

		int max = meshes[i]->get_max_element_id();
    if(solutions_for_adapt) {
      if (errors[i] != NULL) delete [] errors[i];
		  errors[i] = new double[max];
		  memset(errors[i], 0, sizeof(double) * max);
    }
	}

	double total_norm = 0.0;
	double *norms = new double[num];
	memset(norms, 0, num * sizeof(double));
	double *errors_components = new double[num];
	memset(errors_components, 0, num * sizeof(double));
  if(solutions_for_adapt) this->total_err = 0.0;
  double total_error = 0.0;
  if(solutions_for_adapt)
  {
    if (esort != NULL) delete [] esort;
  	esort = new int2[nact];
  }

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
          errors_components[i] += err;
					total_error += err;
          if(solutions_for_adapt)
          {
            this->errors[i][ee[i]->id - 1] += err;
	          this->total_err += err;
          }
				}

			}
		}
	}
	trav.finish();

  // Store the calculation for each solution component separately.
  if(component_errors != NULL) {
    component_errors->clear();
    for (int i = 0; i < num; i++) {
      if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
        component_errors->push_back(sqrt(errors_components[i]));
      else if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
        component_errors->push_back(sqrt(errors_components[i]/norms[i]));
      else {
        error("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
        return -1.0;
      }
    }
  }

  tmr.tick();
  error_time = tmr.accumulated();

  if(solutions_for_adapt)
  {
    k = 0;
    for (i = 0; i < num; i++)
      for(std::map<unsigned int, Element*>::iterator it = meshes[i]->elements.begin(); it != meshes[i]->elements.end(); it++)
		    if (it->second->used)
			    if (it->second->active) {
            Element *e = it->second;
		        esort[k][0] = e->id;
		        esort[k++][1] = i;
		        if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
              errors[i][e->id - 1] /= norms[i];
	        }
    assert(k == nact);
    cmp_err = errors;
    qsort(esort, nact, sizeof(int2), compare);
    // Element error mask is used here, because this variable is used in the adapt()
    // function, where the processed error (sum of errors of processed element errors)
    // is matched to this variable.
    if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
      total_err = total_err / total_norm;

    have_errors = true;
  }
  else {
    for (i = 0; i < n; i++) {
      this->sln[i] = slns_original[i];
      this->rsln[i] = rslns_original[i];
    }
  }  
  delete [] meshes;
  delete [] tr;
  delete [] norms;
  delete [] errors_components;

  // Return error value.
  if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
    return sqrt(total_error);
  else if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
    return sqrt(total_error / total_norm);
  else {
    error("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
    return -1.0;
  }
}
