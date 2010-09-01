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
#include "space.h"
#include "matrix.h"
#include <common/error.h>
#include <common/timer.h>
#include <common/callstack.h>

#define PRINTF(...)
//#define PRINTF printf

#define H3D_INVALID_EDGE_ORDER						-1

#define CHECK_ELEMENT_ID(id) \
	if ((id) < 1 || (id) > mesh->elements.count())\
		EXIT("Invalid element id (eid = %d).", id);\
	assert(mesh->elements.exists(id));


void Space::VertexData::dump(int id) {
	printf("vtx #%d: ced = %d, ", id, ced);
	if (ced) {
		printf("ncomp = %d ", ncomponents);
		for (int i = 0; i < ncomponents; i++) {
			if (i > 0) printf(", ");
			printf("(dof = %d, coef = " SCALAR_FMT ")", baselist[i].dof, SCALAR(baselist[i].coef));
		}
		printf(" ");
	}
	else {
		printf("dof = %d, n = %d", dof, n);
		if (dof == H3D_DIRICHLET_DOF) printf(", bc_proj = " SCALAR_FMT, SCALAR(bc_proj));
	}
	printf("\n");
}

void Space::EdgeData::dump(int id) {
	printf("edge #%d: ced = %d, ", id, ced);
	if (ced) {
		printf("edge_comp = %d", edge_ncomponents);
		for (int i = 0; i < edge_ncomponents; i++) {
			if (i > 0) printf(",");
			printf(" (edge = %ld, ori = %d, part = %d, coef = " SCALAR_FMT ")", edge_baselist[i].edge_id, edge_baselist[i].ori,
				edge_baselist[i].part.part, SCALAR(edge_baselist[i].coef));
		}
		printf(", ");

		printf("face_comp = %d", face_ncomponents);
		for (int i = 0; i < face_ncomponents; i++) {
			if (i > 0) printf(",");
			printf(" (face = %ld, ori = %d, iface = %d, part = (horz = %d, vert = %d), dir = %d, coef = " SCALAR_FMT ")",
				face_baselist[i].face_id, face_baselist[i].ori, face_baselist[i].iface,
				face_baselist[i].part.horz, face_baselist[i].part.vert, face_baselist[i].dir,
				SCALAR(face_baselist[i].coef));
		}
	}
	else {
		printf("order = %d, dof = %d, n = %d", order, dof, n);
		if (bc_proj != NULL) {
			printf(", bc_proj = (");
			for (int i = 0; i < n; i++) {
				if (i > 0) printf(", ");
				printf(SCALAR_FMT, SCALAR(bc_proj[i]));
			}
			printf(")");
		}
	}
	printf("\n");
}

void Space::FaceData::dump(int id) {
	printf("face #%d: ced = %d, ", id, ced);
	if (ced) {
		printf("part = (%d, %d), ori = %d, facet_id = %ld", part.horz, part.vert, ori, facet_id);
	}
	else {
		printf("order = %s, dof = %d, n = %d", order.str(), dof, n);
		if (bc_proj != NULL) {
			printf(", bc_proj = (");
			for (int i = 0; i < n; i++) {
				if (i > 0) printf(", ");
				printf(SCALAR_FMT, SCALAR(bc_proj[i]));
			}
			printf(")");
		}
	}
	printf("\n");
}

void Space::ElementData::dump(int id) {
	printf("elem #%d: ", id);
	printf("order = %s, dof = %d, n = %d", order.str(), dof, n);
	printf("\n");
}

Space::Space(Mesh *mesh, Shapeset *shapeset) :
	mesh(mesh), shapeset(shapeset)
{
	_F_
	set_bc_types(NULL);
	set_essential_bc_values((scalar(*)(int, double, double, double)) NULL);
	set_essential_bc_values((scalar3 &(*)(int, double, double, double)) NULL);
	mesh_seq = -1;
	seq = 0;
	was_assigned = false;

	init_data_tables();
}


Space::~Space() {
	_F_
	free_data_tables();

	for (Word_t i = fi_data.first(); i != INVALID_IDX; i = fi_data.next(i))
		delete fi_data[i];
	fi_data.remove_all();
}

void Space::init_data_tables() {
	_F_
	assert(mesh != NULL);

	FOR_ALL_ELEMENTS(idx, mesh) {
		if (mesh->elements[idx]->active) {
			elm_data[idx] = new ElementData;
			MEM_CHECK(elm_data[idx]);
		}
	}
}

void Space::free_data_tables() {
	_F_
	FOR_ALL_VERTEX_NODES(i)
		delete vn_data[i];
	vn_data.remove_all();

	FOR_ALL_EDGE_NODES(i)
		delete en_data[i];
	en_data.remove_all();

	FOR_ALL_FACE_NODES(i)
		delete fn_data[i];
	fn_data.remove_all();

	FOR_ALL_ELEMENT_NODES(i)
		delete elm_data[i];
	elm_data.remove_all();
}

// element orders ///////////////////////////////////////////////////////////////////////////////


void Space::set_element_order(Word_t eid, order3_t order) {
	_F_
	CHECK_ELEMENT_ID(eid);

	// TODO: check for validity of order
	if (!elm_data.exists(eid)) {
		elm_data[eid] = new ElementData;
		MEM_CHECK(elm_data[eid]);
	}

	assert(mesh->elements[eid]->get_mode() == order.type);
	elm_data[eid]->order = order;
	seq++;
}

order3_t Space::get_element_order(Word_t eid) const {
	_F_
	CHECK_ELEMENT_ID(eid);
	assert(elm_data.exists(eid));
	assert(elm_data[eid] != NULL);
	return elm_data[eid]->order;
}

void Space::set_uniform_order(order3_t order) {
	_F_
	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		assert(elm_data.exists(eid));
		assert(elm_data[eid] != NULL);
		assert(mesh->elements[eid]->get_mode() == order.type);
		elm_data[eid]->order = order;
	}
	seq++;
}

void Space::set_order_recurrent(Word_t eid, order3_t order) {
	_F_
	Element *e = mesh->elements[eid];
	if (e->active) {
		assert(elm_data.exists(e->id));
		assert(elm_data[e->id] != NULL);
		assert(mesh->elements[eid]->get_mode() == order.type);
		elm_data[e->id]->order = order;
	}
	else {
		for (int i = 0; i < e->get_num_sons(); i++) {
			Word_t son = e->get_son(i);
			if (son != INVALID_IDX)
				set_order_recurrent(son, order);
		}
	}
}

inline int LIMIT_ELEMENT_ORDER(int a) {
	_F_
	if (a > H3D_MAX_ELEMENT_ORDER) return H3D_MAX_ELEMENT_ORDER;
	else return a;
}

void Space::copy_orders(const Space &space, int inc) {
	_F_
	Mesh *cmesh = space.get_mesh();
	FOR_ALL_ACTIVE_ELEMENTS(eid, cmesh) {
		order3_t oo = space.get_element_order(eid);
		assert(cmesh->elements[eid]->get_mode() == mesh->elements[eid]->get_mode());

		order3_t order;
		switch (cmesh->elements[eid]->get_mode()) {
			case MODE_TETRAHEDRON: order = oo + order3_t(inc); break;
			case MODE_HEXAHEDRON: order = oo + order3_t(inc, inc, inc); break;
			default: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
		}
		order.limit();

		set_order_recurrent(eid, order);
	}
	seq++;
}

void Space::enforce_minimum_rule() {
	_F_
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *elem = mesh->elements[idx];
		ElementData *elem_node = elm_data[idx];
		order3_t elm_order = elem_node->order;

		switch (elem->get_mode()) {
			case MODE_TETRAHEDRON:
				// on faces
				for (int iface = 0; iface < elem->get_num_faces(); iface++) {
					Word_t fidx = mesh->get_facet_id(elem, iface);
					assert(fn_data.exists(fidx));
					FaceData *fnode = fn_data[fidx];
					order2_t forder = elem_node->order.get_face_order(iface);
					if (fnode->order.invalid() || forder.order < fnode->order.order)
						fnode->order = forder;
				}

				// on edges
				for (int iedge = 0; iedge < elem->get_num_edges(); iedge++) {
					Word_t eidx = mesh->get_edge_id(elem, iedge);
					assert(en_data.exists(eidx));

					EdgeData *enode = en_data[eidx];
					order1_t eorder = elem_node->order.get_edge_order(iedge);
					if (enode->order == H3D_INVALID_EDGE_ORDER || eorder < enode->order)
						enode->order = eorder;
				}
				break;

			case MODE_HEXAHEDRON:
				// on faces
				for (int iface = 0; iface < elem->get_num_faces(); iface++) {
					Word_t fidx = mesh->get_facet_id(elem, iface);
					FaceData *fnode = fn_data[fidx];

					if (!fnode->ced) {
						order2_t forder = elem_node->order.get_face_order(iface);
						if (elem->get_face_orientation(iface) >= 4) forder = order2_t(forder.y, forder.x);		// switch h- and v- order

						if (fnode->order.invalid())
							fnode->order = forder;
						else
							fnode->order = order2_t(std::min(fnode->order.x, forder.x), std::min(fnode->order.y, forder.y));
					}
				}

				// on edges
				for (int iedge = 0; iedge < elem->get_num_edges(); iedge++) {
					Word_t eidx = mesh->get_edge_id(elem, iedge);
					assert(eidx != INVALID_IDX);
					if (mesh->edges[eidx].is_active()) {
						EdgeData *enode = en_data[eidx];
						if (!enode->ced) {
							order1_t eorder = elem_node->order.get_edge_order(iedge);
							if (enode->order == H3D_INVALID_EDGE_ORDER || eorder < enode->order)
								enode->order = eorder;
						}
					}
				}
				break;

			default:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
		}
	}
}

//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void Space::assign_vertex_dofs(Word_t vid) {
	_F_
	VertexData *node = vn_data[vid];
	int ndofs = get_vertex_ndofs();
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = H3D_DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_edge_dofs(Word_t idx) {
	_F_
	EdgeData *node = en_data[idx];
	int ndofs = get_edge_ndofs(node->order);
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = H3D_DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_face_dofs(Word_t idx) {
	_F_
	FaceData *node = fn_data[idx];
	int ndofs = get_face_ndofs(node->order);
	if (node->bc_type == BC_ESSENTIAL) {
		node->dof = H3D_DIRICHLET_DOF;
	}
	else {
		node->dof = next_dof;
		next_dof += ndofs * stride;
	}
	node->n = ndofs;
}

void Space::assign_bubble_dofs(Word_t idx) {
	_F_
	ElementData *enode = elm_data[idx];
	int ndofs = get_element_ndofs(enode->order);
	enode->n = ndofs;
	enode->dof = next_dof;
	next_dof += ndofs * stride;
}

// assembly lists ////

void Space::get_vertex_assembly_list(Element *e, int ivertex, AsmList *al) {
	_F_
	Word_t vtx = e->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	int index = shapeset->get_vertex_index(ivertex);

	if (vnode->ced) {
		for (int i = 0; i < vnode->ncomponents; i++) {
			int dof = vnode->baselist[i].dof;
			assert(dof == H3D_DIRICHLET_DOF || (dof >= first_dof && dof < next_dof));
			al->add(index, dof, vnode->baselist[i].coef);
		}
	}
	else {
		scalar coef = vnode->dof >= 0 ? 1.0 : vnode->bc_proj;
		assert(vnode->dof == H3D_DIRICHLET_DOF || (vnode->dof >= first_dof && vnode->dof < next_dof));
		al->add(index, vnode->dof, coef);
	}
}

void Space::get_edge_assembly_list(Element *elem, int iedge, AsmList *al) {
	_F_
	Word_t edge_id = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge_id];
	int ori = elem->get_edge_orientation(iedge);

	if (enode->ced) {
		// an edge constrained by another edge
		for (int i = 0; i < enode->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = enode->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
			assert(cng_enode->ced == false);

			if (cng_enode->n > 0) {
				int *indices = shapeset->get_edge_indices(iedge, 0, cng_enode->order);		// iedge bude 0 (?)
				if (cng_enode->dof >= 0) {
					for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, dof += stride) {
						order1_t order = shapeset->get_order(indices[j]).get_edge_order(iedge);
						int idx = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
						assert(dof >= first_dof && dof < next_dof);
						al->add(idx, dof, ecomp->coef);
					}
				}
				else {
					for (int j = 0; j < cng_enode->n; j++) {
						order1_t order = shapeset->get_order(indices[j]).get_edge_order(iedge);
						int idx = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
						al->add(idx, H3D_DIRICHLET_DOF, ecomp->coef * cng_enode->bc_proj[j]);
					}
				}
			}
		}
		// an edge constrained by a face
		for (int i = 0; i < enode->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = enode->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining edge node
			assert(cng_fnode->ced == false);

			if (cng_fnode->n > 0) {
				int *indices = shapeset->get_face_indices(fcomp->iface, 0, cng_fnode->order);
				if (cng_fnode->dof >= 0) {
					for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, dof += stride) {
						order2_t order = shapeset->get_order(indices[j]).get_face_order(fcomp->iface);
						int idx = shapeset->get_constrained_edge_face_index(iedge, fcomp->ori, order, fcomp->part, fcomp->dir, shapeset->get_face_fn_variant(indices[j]));
						assert(dof >= first_dof && dof < next_dof);
						al->add(idx, dof, fcomp->coef);
					}
				}
				else {
					for (int j = 0; j < cng_fnode->n; j++) {
						order2_t order = shapeset->get_order(indices[j]).get_face_order(fcomp->iface);
						int idx = shapeset->get_constrained_edge_face_index(iedge, fcomp->ori, order, fcomp->part, fcomp->dir, shapeset->get_face_fn_variant(indices[j]));
						al->add(idx, H3D_DIRICHLET_DOF, fcomp->coef * cng_fnode->bc_proj[j]);
					}
				}
			}
		}
	}
	else {
		if (enode->n > 0) {
			int *indices = shapeset->get_edge_indices(iedge, ori, enode->order);
			if (enode->dof >= 0) {
				for (int j = 0, dof = enode->dof; j < enode->n; j++, dof += stride) {
					assert(dof >= first_dof && dof < next_dof);
					al->add(indices[j], dof, 1.0);
				}
			}
			else if (enode->bc_proj != NULL) {
				for (int j = 0; j < enode->n; j++) {
					scalar coef = enode->bc_proj[j];
					al->add(indices[j], H3D_DIRICHLET_DOF, coef);
				}
			}
		}
	}
}

void Space::get_face_assembly_list(Element *elem, int iface, AsmList *al) {
	_F_
	Word_t face_id = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[face_id];

	if (fnode->ced) {
		if (fnode->facet_id != INVALID_IDX) {
			FaceData *cng_fnode = fn_data[fnode->facet_id];

			if (cng_fnode->n > 0) {
				int *indices = shapeset->get_face_indices(iface, 0, cng_fnode->order);
				if (cng_fnode->dof >= 0) {
					for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, dof += stride) {
						order2_t order = shapeset->get_order(indices[j]).get_face_order(iface);
						int idx = shapeset->get_constrained_face_index(iface, fnode->ori, order, fnode->part, shapeset->get_face_fn_variant(indices[j]));
						assert(dof == H3D_DIRICHLET_DOF || (dof >= first_dof && dof < next_dof));
						al->add(idx, dof, 1.0);
					}
				}
				else {
					assert(false);
				}
			}
		}
	}
	else {
		int ori = elem->get_face_orientation(iface);

		if (fnode->n > 0) {
			int *indices = shapeset->get_face_indices(iface, ori, fnode->order);
			if (fnode->dof >= 0) {
				for (int j = 0, dof = fnode->dof; j < fnode->n; j++, dof += stride) {
					assert(dof >= first_dof && dof < next_dof);
					al->add(indices[j], dof, 1.0);
				}
			}
			else if (fnode->bc_proj != NULL) {
				for (int j = 0; j < fnode->n; j++) {
					scalar coef = fnode->bc_proj[j];
					al->add(indices[j], H3D_DIRICHLET_DOF, coef);
				}
			}
		}
	}
}

void Space::get_bubble_assembly_list(Element *e, AsmList *al) {
	_F_
	ElementData *enode = elm_data[e->id];

	if (enode->n > 0) {
		int *indices = shapeset->get_bubble_indices(enode->order);
		for (int j = 0, dof = enode->dof; j < enode->n; j++, dof += stride) {
			assert(dof >= first_dof && dof < next_dof);
			al->add(indices[j], dof, 1.0);
		}
	}
}

// BC ////

void Space::set_bc_info(NodeData *node, BCType bc, int marker) {
	_F_
	if (bc == BC_ESSENTIAL || (bc == BC_NATURAL && node->bc_type == BC_NONE)) {
		node->bc_type = bc;
		node->marker = marker;
	}
}

void Space::set_bc_information() {
	_F_
	FOR_ALL_FACETS(idx, mesh) {
		Facet *facet = mesh->facets.get(idx);
		assert(facet != NULL);

		if (facet->type == Facet::OUTER) {
			Boundary *bdr = mesh->boundaries[facet->right];
			BCType bc_type = bc_type_callback(bdr->marker);

			int marker = bdr->marker;
			// set boundary condition for face
			assert(fn_data.exists(idx));
			fn_data[idx]->bc_type = bc_type;
			if (fn_data[idx]->marker == H3D_MARKER_UNDEFINED)
				fn_data[idx]->marker = marker;

			Element *elem = mesh->elements[facet->left];
			int iface = facet->left_face_num;

			// set boundary condition for vertices on the face
			int vtx_num = elem->get_num_face_vertices(iface);
			Word_t vtcs[vtx_num];
			elem->get_face_vertices(iface, vtcs);
			for (int i = 0; i < vtx_num; i++) {
				assert(vn_data.exists(vtcs[i]));
				set_bc_info(vn_data[vtcs[i]], bc_type, marker);
			}

			// set boundary condition for edges on the face
			const int *face_edges = elem->get_face_edges(iface);
			for (int i = 0; i < elem->get_num_face_edges(iface); i++) {
				int edge_id = mesh->get_edge_id(elem, face_edges[i]);
				if (mesh->edges[edge_id].bnd) {
					assert(en_data.exists(edge_id));
					set_bc_info(en_data[edge_id], bc_type, marker);
				}
			}
		}
	}
}

// find constraints ///////////////////////////////////////////////////////////

Space::VertexData *Space::create_vertex_node_data(Word_t vid, bool ced) {
	_F_
	VertexData *vd = vn_data[vid];
	if (vd == NULL) {
		vd = vn_data[vid] = new VertexData;
		MEM_CHECK(vd);
		vd->ced = ced;
		if (ced) {
			vd->baselist = NULL;
			vd->ncomponents = 0;
		}
		else {
			vd->dof = H3D_DOF_UNASSIGNED;
			vd->n = -1;
		}
	}
	else {
		// not CED and should be => change it (but ont the other way)
		if (!vd->ced && ced) {
			vd->ced = ced;
			vd->baselist = NULL;
			vd->ncomponents = 0;
		}
	}

	return vd;
}

Space::EdgeData *Space::create_edge_node_data(Word_t eid, bool ced) {
	_F_
	EdgeData *ed = en_data[eid];
	if (ed == NULL) {
		ed = en_data[eid] = new EdgeData;
		MEM_CHECK(ed);
		ed->ced = ced;
		if (ced) {
			ed->edge_baselist = NULL;
			ed->edge_ncomponents = 0;
			ed->face_baselist = NULL;
			ed->face_ncomponents = 0;
		}
		else {
			ed->order = H3D_INVALID_EDGE_ORDER;
			ed->dof = H3D_DOF_UNASSIGNED;
			ed->n = -1;
		}
	}
	else {
		// not CED and should be => change it (but ont the other way)
		if (!ed->ced && ced) {
			ed->ced = ced;
			ed->edge_baselist = NULL;
			ed->edge_ncomponents = 0;
			ed->face_baselist = NULL;
			ed->face_ncomponents = 0;
		}
	}

	return ed;
}

Space::FaceData *Space::create_face_node_data(Word_t fid, bool ced) {
	_F_
	FaceData *fd = fn_data[fid];
	if (fd == NULL) {
		fd = fn_data[fid] = new FaceData;
		MEM_CHECK(fd);
		fd->ced = ced;
		if (ced) {
			fd->facet_id = INVALID_IDX;
			fd->ori = 0;
			fd->part.horz = 0;
			fd->part.vert = 0;
		}
		else {
			fd->dof = H3D_DOF_UNASSIGNED;
			fd->n = -1;
		}
	}
	else {
		if (!fd->ced && ced) {
			fd->ced = ced;
			fd->facet_id = INVALID_IDX;
			fd->ori = 0;
			fd->part.horz = 0;
			fd->part.vert = 0;
		}
	}

	return fd;
}

void Space::fc_face(Word_t eid, int iface, bool ced) {
	_F_

	if (eid == INVALID_IDX) return;

	Element *elem = mesh->elements[eid];
	// vertices
	int nv = elem->get_num_face_vertices(iface);
	Word_t vtcs[nv];
	elem->get_face_vertices(iface, vtcs);

	Word_t fid = mesh->get_facet_id(elem, iface);
	Facet *facet = mesh->facets[fid];

	if (ced) face_ced.set(fid);

	// set CEDs
	Word_t emp[4], fmp;
	switch (facet->ref_mask) {
		case H3D_REFT_QUAD_HORZ:
			emp[0] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[1] = mesh->peek_midpoint(vtcs[3], vtcs[0]);

			// vertices
			create_vertex_node_data(emp[0], ced);
			create_vertex_node_data(emp[1], ced);
			// edges
			create_edge_node_data(mesh->get_edge_id(vtcs[1], emp[0]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[0], vtcs[2]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[3], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[1], vtcs[0]), ced);

			create_edge_node_data(mesh->get_edge_id(vtcs[0], vtcs[1]), false);
			create_edge_node_data(mesh->get_edge_id(emp[0], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[2], vtcs[3]), false);
			break;

		case H3D_REFT_QUAD_VERT:
			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[2], vtcs[3]);

			// vertices
			create_vertex_node_data(emp[0], ced);
			create_vertex_node_data(emp[1], ced);
			// edges by edges
			create_edge_node_data(mesh->get_edge_id(vtcs[0], emp[0]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[0], vtcs[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[2], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(emp[1], vtcs[3]), ced);

			create_edge_node_data(mesh->get_edge_id(vtcs[0], vtcs[3]), false);
			create_edge_node_data(mesh->get_edge_id(emp[0], emp[1]), ced);
			create_edge_node_data(mesh->get_edge_id(vtcs[1], vtcs[2]), false);
			break;

		case H3D_REFT_QUAD_BOTH:
			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[2] = mesh->peek_midpoint(vtcs[2], vtcs[3]);
			emp[3] = mesh->peek_midpoint(vtcs[3], vtcs[0]);
			fmp = mesh->peek_midpoint(emp[0], emp[2]);

			// vertices
			for (int iv = 0; iv < Quad::NUM_VERTICES; iv++)
				create_vertex_node_data(emp[iv], ced);
			create_vertex_node_data(fmp, ced);
			// edges
			for (int i = 0; i < Quad::NUM_VERTICES; i++) {
				int i1 = (i + 1) % Quad::NUM_VERTICES;
				create_edge_node_data(mesh->get_edge_id(vtcs[i], emp[i]), ced);
				create_edge_node_data(mesh->get_edge_id(emp[i], vtcs[i1]), ced);
				create_edge_node_data(mesh->get_edge_id(emp[i], fmp), ced);
			}
			break;
	}

	// faces (common for all types of refinements)
	for (int i = 0; i < Facet::MAX_SONS; i++) {
		Word_t sid = facet->sons[i];
		if (sid != INVALID_IDX) create_face_node_data(sid, ced);
	}
}

void Space::fc_face_left(Word_t fid) {
	_F_
	if (fid == INVALID_IDX) return;

	Facet *facet = mesh->facets[fid];
	fc_face(facet->left, facet->left_face_num, true);
	// recur to sons
	for (int i = 0; i < Facet::MAX_SONS; i++)
		fc_face_left(facet->sons[i]);
}

void Space::fc_face_right(Word_t fid) {
	_F_
	if (fid == INVALID_IDX) return;

	Facet *facet = mesh->facets[fid];
	fc_face(facet->right, facet->right_face_num, true);
	// recur to sons
	for (int i = 0; i < Facet::MAX_SONS; i++)
		fc_face_right(facet->sons[i]);
}

void Space::fc_element(Word_t idx) {
	_F_
	if (idx == INVALID_IDX) return;

	Element *elem = mesh->elements[idx];
	for (int iface = 0; iface < elem->get_num_faces(); iface++) {
		Word_t fid = mesh->get_facet_id(elem, iface);
		Facet *facet = mesh->facets[fid];
		assert(facet != NULL);

		// vertices
		int nv = elem->get_num_face_vertices(iface);
		Word_t vtcs[nv];
		elem->get_face_vertices(iface, vtcs);
		for (int iv = 0; iv < nv; iv++)
			create_vertex_node_data(vtcs[iv], false);
		// edges
		int ne = elem->get_num_face_edges(iface);
		const int *edge_idx = elem->get_face_edges(iface);
		for (int ie = 0; ie < ne; ie++)
			create_edge_node_data(mesh->get_edge_id(elem, edge_idx[ie]), false);

		// face
		create_face_node_data(fid, false);

		// handle possible CEDs
		if (facet->type == Facet::INNER) {
			if (facet->lactive && !facet->ractive) {
				fc_face(facet->left, facet->left_face_num, true);
				fc_face_right(fid);
			}
			else if (!facet->lactive && facet->ractive) {
				fc_face(facet->right, facet->right_face_num, true);
				fc_face_left(fid);
			}
			else if (!facet->lactive && !facet->ractive) {
				// facet sons
				for (int i = 0; i < Facet::MAX_SONS; i++) {
					Word_t son = facet->sons[i];
					if (son != INVALID_IDX) {
						Facet *facet = mesh->facets[son];
						if (son != INVALID_IDX) {
							Facet *son_facet = mesh->facets[son];

							if (son_facet->lactive && !son_facet->ractive) {
								fc_face(facet->left, facet->left_face_num, true);
							}
							else if (!son_facet->lactive && son_facet->ractive) {
								fc_face(facet->right, facet->right_face_num, true);
							}
						}
					}
				}
			}
		}
	}
}

void Space::fc_base(Word_t eid, int iface)
{
	if (eid == INVALID_IDX) return;

	Element *elem = mesh->elements[eid];

	// vertices
	int nv = elem->get_num_face_vertices(iface);
	Word_t vtcs[nv];
	elem->get_face_vertices(iface, vtcs);
	for (int iv = 0; iv < nv; iv++)
		create_vertex_node_data(vtcs[iv], false);
	// edges
	int ne = elem->get_num_face_edges(iface);
	const int *edge_idx = elem->get_face_edges(iface);
	for (int ie = 0; ie < ne; ie++)
		create_edge_node_data(mesh->get_edge_id(elem, edge_idx[ie]), false);

	//
	Word_t fid = mesh->get_facet_id(elem, iface);
	create_face_node_data(fid, false);

}

void Space::find_constraints()
{
	_F_
	face_ced.free();

	// modified breadth-first search
	Array<Word_t> open;
	BitArray elms;

	// first include all base elements
	FOR_ALL_BASE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			if (!elms.is_set(fid)) {
				open.add(fid);
				elms.set(fid);
			}
		}
	}

	for (Word_t idx = open.first(); idx != INVALID_IDX; idx = open.next(idx)) {
		Word_t fid = open[idx];
		Facet *facet = mesh->facets[fid];

		if (facet->left != INVALID_IDX) {
			Element *e = mesh->elements[facet->left];
			for (int iface = 0; iface < e->get_num_faces(); iface++) {
				Word_t fid = mesh->get_facet_id(e, iface);
				if (!elms.is_set(fid)) {
					open.add(fid);
					elms.set(fid);
				}
			}
		}

		if (facet->type == Facet::INNER && facet->right != INVALID_IDX) {
			Element *e = mesh->elements[facet->right];
			for (int iface = 0; iface < e->get_num_faces(); iface++) {
				Word_t fid = mesh->get_facet_id(e, iface);
				if (!elms.is_set(fid)) {
					open.add(fid);
					elms.set(fid);
				}
			}
		}

		for (int i = 0; i < Facet::MAX_SONS; i++) {
			Word_t son = facet->sons[i];
			if (son != INVALID_IDX) {
				if (!elms.is_set(son)) {
					open.add(son);
					elms.set(son);
				}
			}
		}
	}

	for (Word_t idx = open.first(); idx != INVALID_IDX; idx = open.next(idx)) {
		Word_t fid = open[idx];
		Facet *facet = mesh->facets[fid];
		assert(facet != NULL);

		fc_base(facet->left, facet->left_face_num);
		if (facet->type == Facet::INNER)
			fc_base(facet->right, facet->right_face_num);

		// handle possible CEDs
		if (facet->type == Facet::INNER) {
			if (facet->lactive && !facet->ractive) {
				fc_face(facet->left, facet->left_face_num, true);
				fc_face_right(fid);
			}
			else if (!facet->lactive && facet->ractive) {
				fc_face(facet->right, facet->right_face_num, true);
				fc_face_left(fid);
			}
		}
	}
}

inline void Space::output_component(BaseVertexComponent *&current, BaseVertexComponent *&last, BaseVertexComponent *min, bool add) {
	_F_
	// if the edge is already in the list, just add half of the other coef
	if (last != NULL && last->dof == min->dof) {
		PRINTF(" * dof already in list (%d, last->coef = % lf, min->coef = % lf.\n", min->dof, last->coef, min->coef);
		if (add) last->coef += min->coef;
		else last->coef += min->coef * 0.5;
		return;
	}

	// output new vertex component
	current->dof = min->dof;
	if (add) current->coef = min->coef;
	else current->coef = min->coef * 0.5;
	last = current++;
}

/// Duplicate base list
///
/// @param[in] l - list
/// @param[in] n - number of elements in the list
template<typename T>
static T *duplicate_baselist(T *l, int n) {
	_F_
	T *dup = (T *) malloc(n * sizeof(T));
	memcpy(dup, l, n * sizeof(T));
	return dup;
}

Space::BaseVertexComponent *Space::merge_baselist(BaseVertexComponent *l1, int n1, BaseVertexComponent *l2, int n2, int &ncomponents, bool add) {
	_F_
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }
	if (l1 == NULL) { ncomponents = n2; return duplicate_baselist<BaseVertexComponent>(l2, n2); }
	if (l2 == NULL) { ncomponents = n1; return duplicate_baselist<BaseVertexComponent>(l1, n1); }

	// estimate the upper bound of the result size
	int max_result = n1 + n2;

	BaseVertexComponent *result = (BaseVertexComponent *) malloc(max_result * sizeof(BaseVertexComponent));
	BaseVertexComponent *current = result;
	BaseVertexComponent *last = NULL;

	// main loop - always output the component with smaller dof so that we get a sorted array
	int i1 = 0, i2 = 0;
	while (i1 < n1 && i2 < n2) {
		if (l1[i1].dof < l2[i2].dof) output_component(current, last, l1 + i1++, add);
		else output_component(current, last, l2 + i2++, add);
	}

	// finish the longer baselist
	while (i1 < n1) output_component(current, last, l1 + i1++, add);
	while (i2 < n2) output_component(current, last, l2 + i2++, add);

	// if we produced less components than we expected, reallocate the resulting array
	// ...this should be OK as we are always shrinking the array so no copying should occur
	ncomponents = current - result;
	if (ncomponents < max_result) return (BaseVertexComponent *) realloc(result, ncomponents * sizeof(BaseVertexComponent));
	else return result;
}

inline void Space::output_component(BaseEdgeComponent *&current, BaseEdgeComponent *&last, BaseEdgeComponent *min, bool add) {
	_F_
	// if the edge is already in the list, just add half of the other coef
	if (last != NULL && last->edge_id == min->edge_id) {
		PRINTF(" * edge already in list (%ld, last->coef = % lf, min->coef = % lf.\n", min->edge_id, last->coef, min->coef);
		if (add) last->coef += min->coef;
		else last->coef += min->coef * 0.5;
		return;
	}

	// output new edge component
	current->edge_id = min->edge_id;
	current->ori = min->ori;
	current->part = min->part;
	if (add) current->coef = min->coef;
	else current->coef = min->coef * 0.5;
	last = current++;
}

Space::BaseEdgeComponent *Space::merge_baselist(BaseEdgeComponent *l1, int n1, BaseEdgeComponent *l2, int n2, int &ncomponents, bool add) {
	_F_
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }
	if (l1 == NULL) { ncomponents = n2; return duplicate_baselist<BaseEdgeComponent>(l2, n2); }
	if (l2 == NULL) { ncomponents = n1; return duplicate_baselist<BaseEdgeComponent>(l1, n1); }

	int max_result = n1 + n2;
	BaseEdgeComponent *result = (BaseEdgeComponent *) malloc(max_result * sizeof(BaseEdgeComponent));
	BaseEdgeComponent *current = result;
	BaseEdgeComponent *last = NULL;

	// main loop - always output the component with smaller edge_id so that we get a sorted array
	int i1 = 0, i2 = 0;
	while (i1 < n1 && i2 < n2) {
		if (l1[i1].edge_id < l2[i2].edge_id) output_component(current, last, l1 + i1++, add);
		else output_component(current, last, l2 + i2++, add);
	}

	// finish the longer baselist
	while (i1 < n1) output_component(current, last, l1 + i1++, add);
	while (i2 < n2) output_component(current, last, l2 + i2++, add);

	ncomponents = current - result;
	if (ncomponents < max_result) return (BaseEdgeComponent *) realloc(result, ncomponents * sizeof(BaseEdgeComponent));
	else return result;
}

inline void Space::output_component_over(BaseFaceComponent *&current, BaseFaceComponent *min, BaseFaceComponent *m) {
	_F_
	if (min != NULL) {
		current->face_id = min->face_id;
		current->ori = min->ori;
		current->iface = min->iface;
		current->part = min->part;
		current->coef = min->coef;
		current->dir = min->dir;
	}
	else if (m != NULL) {
		current->face_id = m->face_id;
		current->ori = m->ori;
		current->iface = m->iface;
		current->part = m->part;
		current->coef = m->coef;
		current->dir = m->dir;
	}
	current++;
}

inline void Space::output_component(BaseFaceComponent *&current, BaseFaceComponent *&last, BaseFaceComponent *min, bool add) {
	_F_
	if (last != NULL && last->face_id == min->face_id &&
		last->part.vert == min->part.vert && last->part.horz == min->part.horz && last->dir == min->dir) {
		PRINTF(" * face already in list (%ld, last->coef = % lf, min->coef = % lf, last->part = %d, min->part = %d)\n",
			min->face_id, last->coef, min->coef, last->part.horz, min->part.horz);
		last->coef += min->coef * 0.5;
		return;
	}

	// output new face component
	current->face_id = min->face_id;
	current->ori = min->ori;
	current->iface = min->iface;
	current->part = min->part;
	current->dir = min->dir;
	if (add) current->coef = min->coef;
	else current->coef = min->coef * 0.5;
	last = current++;
}

Space::BaseFaceComponent *Space::merge_baselist(BaseFaceComponent *l1, int n1, BaseFaceComponent *l2, int n2, int &ncomponents, Word_t fid, bool add) {
	_F_
	if (l1 == NULL && l2 == NULL) { ncomponents = 0; return NULL; }

	int max_result = n1 + n2;
	BaseFaceComponent *result = (BaseFaceComponent *) malloc(max_result * sizeof(BaseFaceComponent));
	BaseFaceComponent *current = result;
	BaseFaceComponent *last = NULL;

	// main loop - always output the component with smaller face_id so that we get a sorted array
	int i1 = 0, i2 = 0;
	if (add) {
		while (i1 < n1 && i2 < n2) {
			if (l1[i1].face_id < l2[i2].face_id) output_component(current, last, l1 + i1++, add);
			else output_component(current, last, l2 + i2++, add);
		}

		// finish the longer baselist
		while (i1 < n1) output_component(current, last, l1 + i1++, add);
		while (i2 < n2) output_component(current, last, l2 + i2++, add);

		ncomponents = current - result;
	}
	else {
		while (i1 < n1 && i2 < n2) {
			if (l1[i1].face_id == fid) i1++;
			else if (l2[i2].face_id == fid) i2++;
			else {
				if ((l1[i1].face_id < l2[i2].face_id) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz < l2[i2].part.horz) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz == l2[i2].part.horz && l1[i1].part.vert < l2[i2].part.vert) ||
					(l1[i1].face_id == l2[i2].face_id && l1[i1].part.horz == l2[i2].part.horz && l1[i1].part.vert == l2[i2].part.vert && l1[i1].dir < l2[i2].dir))
					output_component(current, last, l1 + i1++, add);
				else
					output_component(current, last, l2 + i2++, add);
			}
		}

		// finish the longer baselist
		while (i1 < n1) {
			if (l1[i1].face_id == fid) i1++;
			else output_component(current, last, l1 + i1++, add);
		}
		while (i2 < n2) {
			if (l2[i2].face_id == fid) i2++;
			else output_component(current, last, l2 + i2++, add);
		}

		ncomponents = current - result;
	}

	if (ncomponents < max_result) return (BaseFaceComponent *) realloc(result, ncomponents * sizeof(BaseFaceComponent));
	else return result;
}

/// @param[in] vtx1 - vertex 1
/// @param[in] vtx2 - vertex 2
///
void Space::calc_vertex_vertex_ced(Word_t vtx1, Word_t vtx2) {
	_F_
	if (type == Hcurl || type == Hdiv || type == L2) return;

	assert(vtx1 != INVALID_IDX);
	assert(vtx2 != INVALID_IDX);
	VertexData *vd[] = { vn_data[vtx1], vn_data[vtx2] };
	Word_t mid_pt = mesh->peek_midpoint(vtx1, vtx2);
	assert(mid_pt != INVALID_IDX);

	PRINTF("calc vertex/vertex #%ld\n", mid_pt);

	VertexData *vd_mid = vn_data[mid_pt];
	assert(vd_mid != NULL);

	BaseVertexComponent *bl[2], dummy_bl[2];	// base lists of vtx1 and vtx2
	int nc[2] = { 0, 0 }; // number of components of bl[0] and bl[1]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (vd[k]->ced) {
			bl[k] = vd[k]->baselist;
			nc[k] = vd[k]->ncomponents;
		}
		else {	// make up an artificial baselist
			dummy_bl[k].dof = vd[k]->dof;
			dummy_bl[k].coef = (vd[k]->dof >= 0) ? 1.0 : vd[k]->bc_proj;
			bl[k] = &dummy_bl[k];
			nc[k] = 1;
		}
	}

	assert(vd_mid->ced == 1);
	::free(vd_mid->baselist);
	int ncomp = 0;
	vd_mid->baselist = merge_baselist(bl[0], nc[0], bl[1], nc[1], ncomp, false);
	vd_mid->ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd_mid->baselist[i].dof, vd_mid->baselist[i].coef);
	}
}

/// @param[in] vtx1 - vertex 1
/// @param[in] vtx2 - vertex 2
///
void Space::calc_mid_vertex_vertex_ced(Word_t mid, Word_t vtx1, Word_t vtx2, Word_t vtx3, Word_t vtx4) {
	_F_
	if (type == Hcurl || type == Hdiv || type == L2) return;

	assert(vtx1 != INVALID_IDX);
	assert(vtx2 != INVALID_IDX);
	assert(vtx3 != INVALID_IDX);
	assert(vtx4 != INVALID_IDX);
	VertexData *vd[] = { vn_data[vtx1], vn_data[vtx2], vn_data[vtx3], vn_data[vtx4] };

	PRINTF("calc mid vertex/vertex #%ld (%ld, %ld, %ld, %ld)\n", mid, vtx1, vtx2, vtx3, vtx4);

	VertexData *vd_mid = vn_data[mid];
	assert(vd_mid != NULL);

    BaseVertexComponent *bl[4], dummy_bl[4];	// base lists of vtx1-4
	int nc[4] = { 0, 0, 0, 0 }; // number of components of bl[0-3]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 4; k++) {
		if (vd[k]->ced) {
			bl[k] = vd[k]->baselist;
			nc[k] = vd[k]->ncomponents;
		}
		else {	// make up an artificial baselist
			dummy_bl[k].dof = vd[k]->dof;
			dummy_bl[k].coef = (vd[k]->dof >= 0) ? 1.0 : vd[k]->bc_proj;
			bl[k] = &dummy_bl[k];
			nc[k] = 1;
		}
	}

	int tmp_nc[] = { 0, 0 };
	BaseVertexComponent *tmp_bl[2];
	tmp_bl[0] = merge_baselist(bl[0], nc[0], bl[2], nc[2], tmp_nc[0], false);
	tmp_bl[1] = merge_baselist(bl[1], nc[1], bl[3], nc[3], tmp_nc[1], false);

	::free(vd_mid->baselist);
	int ncomp = 0;
	vd_mid->baselist = merge_baselist(tmp_bl[0], tmp_nc[0], tmp_bl[1], tmp_nc[1], ncomp, false);
	vd_mid->ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++)
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd_mid->baselist[i].dof, vd_mid->baselist[i].coef);

	::free(tmp_bl[0]);
	::free(tmp_bl[1]);
}

void Space::calc_vertex_edge_ced(Word_t vtx, Word_t eid, int ori, int part) {
	_F_
	if (type == Hcurl || type == Hdiv || type == L2) return;

	PRINTF("calc vertex/edge #%ld\n", vtx);

	PRINTF(" - eid = %ld, part = %d, ori = %d\n", eid, part, ori);

	assert(eid != INVALID_IDX);
	EdgeData *ed = en_data[eid];
	assert(ed != NULL);

	assert(vtx != INVALID_IDX);
	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

	double lo, hi;
	int ncomp = 0;
	if (ed->ced) {
		// count the number of components to merge
		int nc = 0;
		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = ed->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
			nc += cng_enode->n;
		}
		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = ed->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node
			nc += cng_fnode->n;
		}
		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(nc * sizeof(BaseVertexComponent));

		int nci = 0;
		// update the edge part
		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent *ecomp = ed->edge_baselist + i;
			EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

			if (cng_enode->n > 0) {
				int *indices = shapeset->get_edge_indices(0, ecomp->ori, cng_enode->order);
				for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, nci++) {
					order1_t order = shapeset->get_order(indices[j]).get_edge_order(0);
					int idx = shapeset->get_constrained_edge_index(0, ecomp->ori, order, ecomp->part);
					baselist[nci].dof = dof;
					baselist[nci].coef = (ecomp->coef) * shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0);
					if (cng_enode->dof == H3D_DIRICHLET_DOF) baselist[nci].coef *= cng_enode->bc_proj[j];
					else dof += stride;
				}
			}
		}
		// update the face part
		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent *fcomp = ed->face_baselist + i;
			FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node

			if (cng_fnode->n > 0) {
				int *indices = shapeset->get_face_indices(2, fcomp->ori, cng_fnode->order);
				for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, nci++) {
					// FIXME: Hex-specific
					order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
					int idx = shapeset->get_constrained_edge_face_index(0, fcomp->ori, order, fcomp->part, fcomp->dir);

					baselist[nci].dof = dof;
					baselist[nci].coef = (fcomp->coef) * shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0);
					if (cng_fnode->dof == H3D_DIRICHLET_DOF) baselist[nci].coef *= cng_fnode->bc_proj[j];
					else dof += stride;
				}
			}
		}

		BaseVertexComponent *tmp = vd->baselist;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, nc, ncomp, true);
		vd->ncomponents = ncomp;

		::free(tmp);
		::free(baselist);
	}
	else {
		get_interval_part(part, lo, hi);
		double mid = (lo + hi) * 0.5;

		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(ed->n * sizeof(BaseVertexComponent));
		if (ed->n > 0) {
			int *indices = shapeset->get_edge_indices(0, ori, ed->order);
			for (int j = 0, dof = ed->dof; j < ed->n; j++) {
				baselist[j].dof = dof;
				baselist[j].coef = shapeset->get_fn_value(indices[j], mid, -1.0, -1.0, 0);
				if (ed->dof == H3D_DIRICHLET_DOF) baselist[j].coef *= ed->bc_proj[j];
				else dof += stride;
			}
		}

		BaseVertexComponent *tmp = vd->baselist;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, ed->n, ncomp, true);
		vd->ncomponents = ncomp;

		::free(tmp);
		::free(baselist);
	}
}

void Space::calc_mid_vertex_edge_ced(Word_t vtx, Word_t fmp, Word_t eid, int ori, int part) {
	_F_
	if (type == Hcurl || type == Hdiv || type == L2) return;

	PRINTF("calc mid vertex/edge #%ld, [%ld | %ld]\n", vtx, eid, fmp);

	assert(eid != INVALID_IDX);

	EdgeData *ed = en_data[eid];
    BaseEdgeComponent *ebl, edummy_bl;
	int enc = 0;
	BaseFaceComponent *fbl, fdummy_bl;
	int fnc = 0;

	// get baselists of edge node[0] and edge node[1]
	if (ed->ced) {
		ebl = ed->edge_baselist;
		enc = ed->edge_ncomponents;

		fbl = ed->face_baselist;
		fnc = ed->face_ncomponents;
	}
	else {	// make up an artificial baselist (we care about edge id and coef
		edummy_bl.edge_id = eid;
		edummy_bl.ori = ori;
		edummy_bl.part.part = part;
		edummy_bl.coef = 1.0;

		ebl = &edummy_bl;
		enc = 1;

		fbl = NULL;
		fnc = 0;
	}

	int en_comp = enc;
	BaseEdgeComponent *edge_baselist = ebl;

	int fn_comp = fnc;
	BaseFaceComponent *face_baselist = fbl;

	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent ec = edge_baselist[i];
		PRINTF(" - [%d]: edge_id = %ld, ori = %d, part = %d, coef = %lf\n",
			i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent fc = face_baselist[i];
		PRINTF(" - [%d]: face_id = %ld, part = (%d, %d), dir = %d, ori = %d, coef = %lf\n",
			i, fc.face_id, fc.part.horz, fc.part.vert, fc.dir, fc.ori, fc.coef);
	}

	// -- /////////////////

	// count the number of components to update
	int nc = 0;
	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent *ecomp = edge_baselist + i;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node
		nc += cng_enode->n;
		PRINTF(" - enc = %d\n", cng_enode->n);
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent *fcomp = face_baselist + i;
		FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining face node
		nc += cng_fnode->n;
		PRINTF(" - fnc = %d\n", cng_fnode->n);
	}

	BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(nc * sizeof(BaseVertexComponent));
	PRINTF(" - nc = %d\n", nc);

	assert(vtx != INVALID_IDX);
	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

	int nci = 0;
	for (int i = 0; i < en_comp; i++) {
		BaseEdgeComponent *ecomp = edge_baselist + i;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

		PRINTF("   - ECOMP: edge_id = %ld, part = %d, coef = %lf, order = %d\n",
			ecomp->edge_id, ecomp->part.part, ecomp->coef, cng_enode->order);

		if (cng_enode->n > 0) {
			int *indices = shapeset->get_edge_indices(0, ecomp->ori, cng_enode->order);		// iedge bude 0 (?)
			for (int j = 0, dof = cng_enode->dof; j < cng_enode->n; j++, nci++) {
				int order = shapeset->get_order(indices[j]).get_edge_order(0);
				int idx = shapeset->get_constrained_edge_index(0, ecomp->ori, order, ecomp->part);
				baselist[nci].dof = dof;
				baselist[nci].coef = ecomp->coef * shapeset->get_fn_value(idx, 0, -1.0, -1.0, 0);
				if (cng_enode->dof == H3D_DIRICHLET_DOF) baselist[nci].coef *= cng_enode->bc_proj[j];
				else dof += stride;
			}
		}
	}
	for (int i = 0; i < fn_comp; i++) {
		BaseFaceComponent *fcomp = face_baselist + i;
		FaceData *cng_fnode = fn_data[fcomp->face_id]; 						// constraining edge node

		PRINTF("   - FCOMP = %ld, part = (%d, %d), dir = %d, ori = %d, coef = %lf\n",
			fcomp->face_id, fcomp->part.horz, fcomp->part.vert, fcomp->dir, fcomp->ori, fcomp->coef);

		if (cng_fnode->n > 0) {
			int *indices = shapeset->get_face_indices(2, fcomp->ori, cng_fnode->order);
			for (int j = 0, dof = cng_fnode->dof; j < cng_fnode->n; j++, nci++) {
				// FIXME: Hex-specific
				order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
				int idx = shapeset->get_constrained_edge_face_index(0, fcomp->ori, order, fcomp->part, fcomp->dir);

				baselist[nci].dof = dof;
				baselist[nci].coef = fcomp->coef * shapeset->get_fn_value(idx, 0.0, -1.0, -1.0, 0);
				if (cng_fnode->dof == H3D_DIRICHLET_DOF) baselist[nci].coef *= cng_fnode->bc_proj[j];
				else dof += stride;
			}
		}
	}

	BaseVertexComponent *tmp = vd->baselist;
	int ncomp = 0;
	vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, nc, ncomp, true);
	vd->ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++)
		PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);

	::free(tmp);

	// ----
	assert(fmp != INVALID_IDX);
	VertexData *fmp_vd = vn_data[fmp];
	assert(fmp_vd != NULL);

	//
	for (int i = 0; i < nc; i++)
		baselist[i].coef *= 0.5;

	BaseVertexComponent *bl, dummy_bl;	// base lists of vtx1-4
	fnc = 0; // number of components of bl[0-3]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	if (fmp_vd->ced) {
		bl = fmp_vd->baselist;
		fnc = fmp_vd->ncomponents;
	}
	else {	// make up an artificial baselist
		dummy_bl.dof = fmp_vd->dof;
		dummy_bl.coef = (fmp_vd->dof >= 0) ? 1.0 : fmp_vd->bc_proj;
		bl = &dummy_bl;
		fnc = 1;
	}

	tmp = fmp_vd->baselist;
	fmp_vd->baselist = merge_baselist(bl, fnc, baselist, nc, ncomp, true);
	fmp_vd->ncomponents = ncomp;

	::free(tmp);
	::free(baselist);
}


/// @param vtx - id of the vertex node for which we calculate the constrain
/// @param fid - id of the constraining facet
/// @param ori - the orientation of the constraining facet
void Space::calc_vertex_face_ced(Word_t vtx, Word_t fid, int ori, int iface, int hpart, int vpart) {
	_F_
	if (type == Hcurl || type == Hdiv || type == L2) return;

	PRINTF("calc vertex/face #%ld\n", vtx);

	FaceData *fd = fn_data[fid];
	assert(fd != NULL);

	VertexData *vd = vn_data[vtx];
	assert(vd != NULL);

	double h_lo, h_hi, v_lo, v_hi;
	get_interval_part(hpart, h_lo, h_hi);
	get_interval_part(vpart, v_lo, v_hi);

	PRINTF(" - fid = %ld, ori = %d, iface = %d, part = (%d, %d), comp = %d\n", fid, ori, iface, hpart, vpart, vd->ncomponents);

	if (fd->ced) {
		EXIT("Unusual vertex/face CED situation, please report.");
	}
	else {
		BaseVertexComponent *baselist = (BaseVertexComponent *) malloc(fd->n * sizeof(BaseVertexComponent));

		if (fd->n > 0) {
			int *indices = shapeset->get_face_indices(2, ori, fd->order);
			for (int j = 0, dof = fd->dof; j < fd->n; j++) {
				// FIXME: hex-specific
				order2_t order = shapeset->get_order(indices[j]).get_face_order(2);
				Part part;
				part.horz = hpart;
				part.vert = vpart;
				int idx = shapeset->get_constrained_face_index(2, ori, order, part, shapeset->get_face_fn_variant(indices[j]));

				baselist[j].dof = dof;
				baselist[j].coef = shapeset->get_fn_value(idx, 0.0, -1.0, 0.0,  0);
				if (fd->dof == H3D_DIRICHLET_DOF) baselist[j].coef *= fd->bc_proj[j];
				else dof += stride;
				PRINTF(" - [%d]: dof = %d, coef = %lf\n", j, baselist[j].dof, baselist[j].coef);
			}
		}

		BaseVertexComponent *tmp = vd->baselist;
		int ncomp = 0;
		vd->baselist = merge_baselist(vd->baselist, vd->ncomponents, baselist, fd->n, ncomp, true);
		vd->ncomponents = ncomp;

		PRINTF("--\n");
		for (int i = 0; i < vd->ncomponents; i++)
			PRINTF(" - [%d]: dof = %d, coef = %lf\n", i, vd->baselist[i].dof, vd->baselist[i].coef);

		::free(tmp);
		::free(baselist);
	}
}

/// @param[in] seid - small edge id
/// @param[in] eid - big edge id
/// @param[in] ori - prt orientation
/// @param[in] epart - part of the edge
/// @param[in] part - part of the edge
void Space::calc_edge_edge_ced(Word_t seid, Word_t eid, int ori, int epart, int part) {
	_F_
	if (type == Hdiv || type == L2) return;

	PRINTF("calc edge/edge #%ld, #%ld\n", seid, eid);

	assert(eid != INVALID_IDX);
	EdgeData *cng_ed = en_data[eid]; // constraining edge
	assert(cng_ed != NULL);

	assert(seid != INVALID_IDX);
	EdgeData *ed = en_data[seid];
	assert(ed != NULL);

	PRINTF(" *** part = %d, cng_eid = %ld\n", part, eid);

	if (cng_ed->ced) {
		PRINTF(" - CED version\n");
		// use the baselist from the "parent" edge and update the part and ori (?)
		int ncomp = cng_ed->edge_ncomponents;
		BaseEdgeComponent *edge_bl = (BaseEdgeComponent *) malloc(ncomp * sizeof(BaseEdgeComponent));
		for (int i = 0; i < ncomp; i++) {
			edge_bl[i] = cng_ed->edge_baselist[i];
			edge_bl[i].part.part = combine_face_part(edge_bl[i].part.part, epart);
		}
		::free(ed->edge_baselist);
		ed->edge_baselist = edge_bl;
		ed->edge_ncomponents = ncomp;

		// face components
		ncomp = cng_ed->face_ncomponents;
		BaseFaceComponent *face_bl = (BaseFaceComponent *) malloc(ncomp * sizeof(BaseFaceComponent));
		for (int i = 0; i < ncomp; i++) {
			face_bl[i] = cng_ed->face_baselist[i];

			if (face_bl[i].dir == PART_ORI_VERT) face_bl[i].part.vert = combine_face_part(face_bl[i].part.vert, epart);
			else face_bl[i].part.horz = combine_face_part(face_bl[i].part.horz, epart);
		}
		::free(ed->face_baselist);
		ed->face_baselist = face_bl;
		ed->face_ncomponents = ncomp;

		for (int i = 0; i < ed->edge_ncomponents; i++) {
			BaseEdgeComponent ec = ed->edge_baselist[i];
			PRINTF(" - [%d]: edge_id = %ld, ori = %d, part = %d, coef = %lf\n",
				i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
		}

		for (int i = 0; i < ed->face_ncomponents; i++) {
			BaseFaceComponent fc = ed->face_baselist[i];
			PRINTF(" - [%d]: face_id = %ld, ori = %d, part = (%d, %d), dir = %d, coef = %lf\n",
				i, fc.face_id, fc.ori, fc.part.horz, fc.part.vert, fc.dir, fc.coef);
		}

	}
	else {
		int nc = 1;
		BaseEdgeComponent *baselist = (BaseEdgeComponent *) malloc(nc * sizeof(BaseEdgeComponent));
		baselist[0].edge_id = eid;
		baselist[0].ori = ori;
		baselist[0].part.part = part;
		baselist[0].coef = 1.0;

		assert(ed->ced == 1);
		BaseEdgeComponent *tmp = ed->edge_baselist;
		int ncomp = 0;
		ed->edge_baselist = merge_baselist(ed->edge_baselist, ed->edge_ncomponents, baselist, nc, ncomp, false);
		ed->edge_ncomponents = ncomp;

		for (int i = 0; i < ncomp; i++) {
			BaseEdgeComponent ec = ed->edge_baselist[i];
			PRINTF(" - [%d]: edge_id = %ld, ori = %d, part = %d, coef = %lf\n",
				i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
		}

		::free(tmp);
		::free(baselist);
	}
}

void Space::calc_mid_edge_edge_ced(Word_t meid, Word_t eid[], int ori[], int epart, int part) {
	_F_
	if (type == Hdiv || type == L2) return;

	PRINTF("calc mid edge/edge #%ld\n", meid);

	assert(eid[0] != INVALID_IDX);
	assert(eid[1] != INVALID_IDX);

	assert(meid != INVALID_IDX);
	EdgeData *mid_ed = en_data[meid];
	assert(mid_ed != NULL);

	EdgeData *ed[] = { en_data[eid[0]], en_data[eid[1]] };
    BaseEdgeComponent *bl[2], dummy_bl[2];		// base lists of eid1 and eid2
	int nc[2] = { 0, 0 };						// number of components of bl[0] and bl[1]
	bool free_bl[2];							// tru if the bl[i] should be freed

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (ed[k]->ced) {
			PRINTF(" - CED version\n");
			// use the baselist from the "parent" edge and update the part and ori (?)
			int ncomp = ed[k]->edge_ncomponents;
			BaseEdgeComponent *edge_bl = (BaseEdgeComponent *) malloc(ncomp * sizeof(BaseEdgeComponent));
			for (int i = 0; i < ncomp; i++) {
				edge_bl[i] = ed[k]->edge_baselist[i];
				edge_bl[i].part.part = combine_face_part(edge_bl[i].part.part, epart);
			}

			bl[k] = edge_bl;
			nc[k] = ncomp;
			free_bl[k] = true;
		}
		else {	// make up an artificial baselist
			dummy_bl[k].edge_id = eid[k];
			dummy_bl[k].ori = ori[k];
			dummy_bl[k].part.part = part;
			dummy_bl[k].coef = 1.0;

			bl[k] = &dummy_bl[k];
			nc[k] = 1;
			free_bl[k] = false;
		}
	}

	int ncomp = 0;
	mid_ed->edge_baselist = merge_baselist(bl[0], nc[0], bl[1], nc[1], ncomp, false);
	mid_ed->edge_ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		BaseEdgeComponent ec = mid_ed->edge_baselist[i];
		PRINTF(" - [%d]: edge_id = %ld, ori = %d, part = %d, coef = %lf\n",
			i, ec.edge_id, ec.ori, ec.part.part, ec.coef);
	}

	if (free_bl[0]) ::free(bl[0]);
	if (free_bl[1]) ::free(bl[1]);
}


/// @param[in] eid - edge id
/// @param[in] fid - constraining face id
/// @param[in] ori - orientation of the constraining face
/// @param[in] part_ori - part ori
/// @param[in] fpart - face part
/// @param[in] epart - edge part
void Space::calc_edge_face_ced(Word_t mid_eid, Word_t eid[], Word_t fid, int ori, int iface, int part_ori, int fpart, int epart) {
	_F_
	if (type == Hdiv || type == L2) return;

	PRINTF("calc edge/face #%ld\n", mid_eid);

	assert(fid != INVALID_IDX);
	FaceData *cng_fnode = fn_data[fid]; // constraining face
	assert(cng_fnode != NULL);

	assert(mid_eid != INVALID_IDX);
	EdgeData *mid_ed = en_data[mid_eid];
	assert(mid_ed != NULL);

	PRINTF(" - n = %d | face_id = %ld, ori = %d, iface = %d, part (%d, %d), dir = %d\n",
		cng_fnode->n, fid, ori, iface, fpart, epart, part_ori);

	EdgeData *ed[] = { en_data[eid[0]], en_data[eid[1]] };
    BaseFaceComponent *bl[2];					// base lists of eid1 and eid2
	int nc[2] = { 0, 0 };						// number of components of bl[0] and bl[1]

	// get baselists of vn[0] and vn[1] - pretend we have them even if they are unconstrained
	for (int k = 0; k < 2; k++) {
		if (ed[k]->ced) {
			bl[k] = ed[k]->face_baselist;
			nc[k] = ed[k]->face_ncomponents;
		}
		else {
			bl[k] = NULL;
			nc[k] = 0;
		}
	}

	int mnc = 0;
	BaseFaceComponent *mbl = merge_baselist(bl[0], nc[0], bl[1], nc[1], mnc, fid, false);

	BaseFaceComponent dummy_bl;
	dummy_bl.face_id = fid;
	dummy_bl.ori = ori;
	dummy_bl.iface = iface;
	dummy_bl.part.horz = fpart;
	dummy_bl.part.vert = epart;
	dummy_bl.dir = part_ori;
	dummy_bl.coef = 1.0;

	::free(mid_ed->face_baselist);
	int ncomp = 0;
	mid_ed->face_baselist = merge_baselist(&dummy_bl, 1, mbl, mnc, ncomp, fid, true);
	mid_ed->face_ncomponents = ncomp;

	for (int i = 0; i < ncomp; i++) {
		BaseFaceComponent fc = mid_ed->face_baselist[i];
		PRINTF(" - [%d]: face_id = %ld, ori = %d, iface = %d, part = (%d, %d), dir = %d, coef = %lf\n",
			i, fc.face_id, fc.ori, fc.iface, fc.part.horz, fc.part.vert, fc.dir, fc.coef);
	}

	::free(mbl);
}

/// @param[in] sfid - small facet id (constrained)
/// @param[in] fid - constraining facet id
/// @param[in] ori - orientation of the constraining face
/// @param[in] hpart - horizontal part
/// @param[in] vpart - vertical part
void Space::calc_face_face_ced(Word_t sfid, Word_t fid, int ori, int hpart, int vpart) {
	_F_
	if (type == L2) return;

	PRINTF("calc face/face #%ld\n", sfid);

	FaceData *fd = fn_data[sfid];
	assert(fd != NULL);
	fd->facet_id = fid;
	fd->ori = ori;
	fd->part.horz = hpart;
	fd->part.vert = vpart;

	PRINTF(" - part = (%d, %d)\n", hpart, vpart);
}


void Space::uc_face(Word_t eid, int iface) {
	_F_
	Element *elem = mesh->elements[eid];

	Word_t fid = mesh->get_facet_id(elem, iface);
	assert(fid != INVALID_IDX);
	if (!fi_data.exists(fid)) return;

	FaceInfo *fi = fi_data[fid];
	assert(fi != NULL);

	Facet *facet = mesh->facets[fid];
	assert(facet != NULL);

	// vertices
	int nv = elem->get_num_face_vertices(iface);
	Word_t vtcs[nv];
	elem->get_face_vertices(iface, vtcs);
	// face edges
	const int *face_edge = elem->get_face_edges(iface);

	assert(fi->elem_id != INVALID_IDX);
	Element *big_elem = mesh->elements[fi->elem_id];

	int cng_face_id = mesh->get_facet_id(big_elem, fi->face);
	int cng_face_ori = big_elem->get_face_orientation(fi->face);

	Word_t emp[4], fmp;		// four edge mid-points, one face mid-point
	Word_t cng_edge_id;		// constraining edge id
	int cng_edge_ori;		// orientation of the constraining edge
	Word_t edge_id[4];		// ID of two edges (left-right | upper-lower)
	int edge_ori[4];		// orientation of two edges

	Word_t sub_fid[4], mid_edge_id;
	FaceInfo *sfi, *sub_fi[4];
	int part_ori;			// orientation of edge/face constraint

	switch (facet->ref_mask) {
		case H3D_REFT_QUAD_HORZ:
			PRINTF("HORZ\n");

			emp[0] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[1] = mesh->peek_midpoint(vtcs[3], vtcs[0]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], vtcs[1], emp[0], emp[1]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = fi->h_part;
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[1], emp[0], vtcs[2], vtcs[3]);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = fi->h_part;
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[1]] = sfi;

			// --- ////

			// vertex
			calc_vertex_vertex_ced(vtcs[1], vtcs[2]);
			calc_vertex_vertex_ced(vtcs[3], vtcs[0]);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);

			calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->v_part);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[1], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->v_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);

			calc_vertex_edge_ced(emp[1], cng_edge_id, cng_edge_ori, fi->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[3]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->v_part);

			// mid edge
			mid_edge_id = mesh->get_edge_id(emp[0], emp[1]);
			edge_id[0] = mesh->get_edge_id(vtcs[0], vtcs[1]);
			edge_id[1] = mesh->get_edge_id(vtcs[2], vtcs[3]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[0]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[2]);

			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 0, fi->h_part);

			// edge by face
			// face orientations 4-7 have switched vertical and horizontal orientation (refer to the Pavel Solin's grey book)
			part_ori = PART_ORI_HORZ;
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, fi->v_part);

			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			break;

		case H3D_REFT_QUAD_VERT:
			PRINTF("VERT\n");

			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[2], vtcs[3]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], emp[0], emp[1], vtcs[3]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = fi->v_part;
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[0], vtcs[1], vtcs[2], emp[1]);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = fi->v_part;
			fi_data[sub_fid[1]] = sfi;

			// --- ////

			// vertex
			calc_vertex_vertex_ced(vtcs[0], vtcs[1]);
			calc_vertex_vertex_ced(vtcs[2], vtcs[3]);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);

			calc_vertex_edge_ced(emp[0], cng_edge_id, cng_edge_ori, fi->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[1]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);

			calc_vertex_edge_ced(emp[1], cng_edge_id, cng_edge_ori, fi->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			// mid edge
			mid_edge_id = mesh->get_edge_id(emp[0], emp[1]);
			edge_id[0] = mesh->get_edge_id(elem, face_edge[3]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[1]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[3]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[1]);

			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 0, fi->v_part);

			// edge by face
			// face orientations 4-7 have switched vertical and horizontal orientation (refer to the Pavel Solin's grey book)
			part_ori = PART_ORI_VERT;
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, fi->v_part);

			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			break;

		case H3D_REFT_QUAD_BOTH:
			PRINTF("BOTH\n");

			emp[0] = mesh->peek_midpoint(vtcs[0], vtcs[1]);
			emp[1] = mesh->peek_midpoint(vtcs[1], vtcs[2]);
			emp[2] = mesh->peek_midpoint(vtcs[2], vtcs[3]);
			emp[3] = mesh->peek_midpoint(vtcs[3], vtcs[0]);
			fmp = mesh->peek_midpoint(emp[0], emp[2]);

			// faces
			sub_fid[0] = mesh->get_facet_id(4, vtcs[0], emp[0], fmp, emp[3]);
			sub_fi[0] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[0]] = sfi;

			sub_fid[1] = mesh->get_facet_id(4, emp[0], vtcs[1], emp[1], fmp);
			sub_fi[1] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = get_lower_part(fi->v_part);
			fi_data[sub_fid[1]] = sfi;

			sub_fid[2] = mesh->get_facet_id(4, fmp, emp[1], vtcs[2], emp[2]);
			sub_fi[2] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_higher_part(fi->h_part);
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[2]] = sfi;

			sub_fid[3] = mesh->get_facet_id(4, emp[3], fmp, emp[2], vtcs[3]);
			sub_fi[3] = sfi = new FaceInfo(MODE_QUAD, fi->elem_id, fi->face);
			MEM_CHECK(sfi);
			sfi->h_part = get_lower_part(fi->h_part);
			sfi->v_part = get_higher_part(fi->v_part);
			fi_data[sub_fid[3]] = sfi;

			// --- ////

			// vertex by vertex
			calc_vertex_vertex_ced(vtcs[0], vtcs[1]);
			calc_vertex_vertex_ced(vtcs[1], vtcs[2]);
			calc_vertex_vertex_ced(vtcs[2], vtcs[3]);
			calc_vertex_vertex_ced(vtcs[3], vtcs[0]);
			calc_mid_vertex_vertex_ced(fmp, emp[0], emp[1], emp[2], emp[3]);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);
			calc_mid_vertex_edge_ced(emp[0], fmp, cng_edge_id, cng_edge_ori, fi->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);
			calc_mid_vertex_edge_ced(emp[1], fmp, cng_edge_id, cng_edge_ori, fi->v_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);
			calc_mid_vertex_edge_ced(emp[2], fmp, cng_edge_id, cng_edge_ori, fi->h_part);

			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);
			calc_mid_vertex_edge_ced(emp[3], fmp, cng_edge_id, cng_edge_ori, fi->v_part);

			// vertex by face
			calc_vertex_face_ced(fmp, cng_face_id, cng_face_ori, iface, fi->h_part, fi->v_part);

			// edge by edge
			cng_edge_id = mesh->get_edge_id(elem, face_edge[0]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[0]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[0], emp[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[0], vtcs[1]), cng_edge_id, cng_edge_ori, 2, sub_fi[1]->h_part);

			// egde by edge (right)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[1]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[1]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[1], emp[1]), cng_edge_id, cng_edge_ori, 1, sub_fi[1]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[1], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[2]->v_part);

			// edge by edge (upper)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[2]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[2]);

			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[2]), cng_edge_id, cng_edge_ori, 1, sub_fi[3]->h_part);
			calc_edge_edge_ced(mesh->get_edge_id(emp[2], vtcs[2]), cng_edge_id, cng_edge_ori, 2, sub_fi[2]->h_part);

			// edge by edge (left)
			cng_edge_id = mesh->get_edge_id(elem, face_edge[3]);
			cng_edge_ori = elem->get_edge_orientation(face_edge[3]);

			calc_edge_edge_ced(mesh->get_edge_id(emp[3], vtcs[0]), cng_edge_id, cng_edge_ori, 1, sub_fi[0]->v_part);
			calc_edge_edge_ced(mesh->get_edge_id(vtcs[3], emp[3]), cng_edge_id, cng_edge_ori, 2, sub_fi[3]->v_part);

			// mid egdes (vert)
			edge_id[0] = mesh->get_edge_id(elem, face_edge[3]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[1]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[3]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[1]);

			mid_edge_id = mesh->get_edge_id(emp[0], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 1, sub_fi[0]->v_part);

			mid_edge_id = mesh->get_edge_id(emp[2], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 2, sub_fi[3]->v_part);

			// mid edges (horz)
			edge_id[0] = mesh->get_edge_id(elem, face_edge[0]);
			edge_id[1] = mesh->get_edge_id(elem, face_edge[2]);
			edge_ori[0] = elem->get_edge_orientation(face_edge[0]);
			edge_ori[1] = elem->get_edge_orientation(face_edge[2]);
			//

			mid_edge_id = mesh->get_edge_id(emp[3], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 1, sub_fi[0]->h_part);
			//
			mid_edge_id = mesh->get_edge_id(emp[1], fmp);
			calc_mid_edge_edge_ced(mid_edge_id, edge_id, edge_ori, 2, sub_fi[1]->h_part);

			// edge by face
			part_ori = PART_ORI_VERT;

			edge_id[0] = mesh->get_edge_id(vtcs[0], emp[3]);
			edge_id[1] = mesh->get_edge_id(vtcs[1], emp[1]);
			mid_edge_id = mesh->get_edge_id(emp[0], fmp);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, sub_fi[0]->v_part);

			edge_id[0] = mesh->get_edge_id(emp[3], vtcs[3]);
			edge_id[1] = mesh->get_edge_id(emp[1], vtcs[2]);
			mid_edge_id = mesh->get_edge_id(emp[2], fmp);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, fi->h_part, sub_fi[3]->v_part);

			part_ori = PART_ORI_HORZ;

			edge_id[0] = mesh->get_edge_id(vtcs[0], emp[0]);
			edge_id[1] = mesh->get_edge_id(vtcs[3], emp[2]);
			mid_edge_id = mesh->get_edge_id(emp[3], fmp);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[0]->h_part, fi->v_part);

			edge_id[0] = mesh->get_edge_id(emp[0], vtcs[1]);
			edge_id[1] = mesh->get_edge_id(emp[2], vtcs[2]);
			mid_edge_id = mesh->get_edge_id(fmp, emp[1]);
			calc_edge_face_ced(mid_edge_id, edge_id, cng_face_id, cng_face_ori, iface, part_ori, sub_fi[1]->h_part, fi->v_part);

			// face by face
			calc_face_face_ced(sub_fid[0], cng_face_id, cng_face_ori, sub_fi[0]->h_part, sub_fi[0]->v_part);
			calc_face_face_ced(sub_fid[1], cng_face_id, cng_face_ori, sub_fi[1]->h_part, sub_fi[1]->v_part);
			calc_face_face_ced(sub_fid[2], cng_face_id, cng_face_ori, sub_fi[2]->h_part, sub_fi[2]->v_part);
			calc_face_face_ced(sub_fid[3], cng_face_id, cng_face_ori, sub_fi[3]->h_part, sub_fi[3]->v_part);

			break;
	}
}

void Space::uc_element(Word_t idx) {
	_F_
	if (idx == INVALID_IDX) return;

	Element *e = mesh->elements[idx];


	for (int iface = 0; iface < e->get_num_faces(); iface++) {
		Word_t fid = mesh->get_facet_id(e, iface);
		Facet *facet = mesh->facets[fid];

		const int *edge = e->get_face_edges(iface);
		for (int iedge = 0; iedge < e->get_num_face_edges(iface); iedge++) {
			Word_t edge_id = mesh->get_edge_id(e, edge[iedge]);
			Edge edg = mesh->edges[edge_id];

			if (edg.is_active())
				calc_edge_boundary_projection(e, edge[iedge]);
		}

		if (facet->ractive && facet->lactive && mesh->facets[fid]->type == Facet::OUTER)
			calc_face_boundary_projection(e, iface);

		if (face_ced.is_set(fid)) {
			if (!fi_data.exists(fid)) {
				switch (facet->mode) {
					case MODE_QUAD:
						fi_data[fid] = new FaceInfo(MODE_QUAD, idx, iface);
						MEM_CHECK(fi_data[fid]);
						break;

					case MODE_TRIANGLE:
						EXIT(H3D_ERR_NOT_IMPLEMENTED);
						break;

					default:
						EXIT(H3D_ERR_UNKNOWN_MODE);
						break;
				}
			}

			// face
			uc_face(idx, iface);
		}
	}
}

//

int Space::assign_dofs(int first_dof, int stride) {
	_F_
	this->first_dof = next_dof = first_dof;
	this->stride = stride;

	// free data
	FOR_ALL_VERTEX_NODES(i)
		delete vn_data[i];
	vn_data.remove_all();

	FOR_ALL_EDGE_NODES(i)
		delete en_data[i];
	en_data.remove_all();

	FOR_ALL_FACE_NODES(i)
		delete fn_data[i];
	fn_data.remove_all();

	for (Word_t i = fi_data.first(); i != INVALID_IDX; i = fi_data.next(i))
		delete fi_data[i];
	fi_data.remove_all();

	// find constraints
	find_constraints();

	enforce_minimum_rule();
	set_bc_information();

	assign_dofs_internal();
	update_constraints();

	mesh_seq = mesh->get_seq();
	was_assigned = true;
	seq++;

	return get_dof_count();
}


void Space::uc_dep(Word_t eid)
{
	_F_
	// find all direct dependencies and include them into deps array
	// dependencies already solved are not included (flags are kept in uc_deps array)
#define H3D_MAX_DEP				1000
	Word_t deps[H3D_MAX_DEP];
	int idep = 0;

	Element *e = mesh->elements[eid];
	for (int iface = 0; iface < e->get_num_faces(); iface++) {
		Word_t fid = mesh->get_facet_id(e, iface);
		Facet *facet = mesh->facets[fid];

		if (facet->type == Facet::OUTER) {
			Word_t parent_id = facet->parent;
			if (parent_id != INVALID_IDX) {
				Facet *parent = mesh->facets[parent_id];
				if (!uc_deps.is_set(parent->left) && parent->left != INVALID_IDX) {
					deps[idep++] = parent->left;
					uc_deps.set(parent->left);
				}
			}
		}
		else {
			Word_t parent_id = facet->parent;
			if (parent_id != INVALID_IDX) {
				Facet *parent = mesh->facets[parent_id];
				if (parent->type == Facet::INNER && (parent->left == INVALID_IDX || parent->right == INVALID_IDX)) {
					// CED -> take the value from the other side (i.e. constraining one)
					if (facet->left == eid) {
						if (!uc_deps.is_set(parent->right) && parent->right != INVALID_IDX) {
							deps[idep++] = parent->right;
							uc_deps.set(parent->right);
						}
					}
					else {
						if (!uc_deps.is_set(parent->left) && parent->left != INVALID_IDX) {
							deps[idep++] = parent->left;
							uc_deps.set(parent->left);
						}
					}
				}
				else {
					// no CED, take tha parent element
					if (facet->left == eid) {
						if (!uc_deps.is_set(parent->left) && parent->left != INVALID_IDX) {
							deps[idep++] = parent->left;
							uc_deps.set(parent->left);
						}
					}
					else {
						if (!uc_deps.is_set(parent->right) && parent->right != INVALID_IDX) {
							deps[idep++] = parent->right;
							uc_deps.set(parent->right);
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < idep; i++)
		uc_dep(deps[i]);

	uc_element(eid);
	uc_deps.set(eid);
}

void Space::update_constraints()
{
	_F_
	uc_deps.free();
	// first calc BC projs in all vertices
	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::OUTER) {
				// mark the vertices on the boundary
				const int *vtx = e->get_face_vertices(iface);
				for (int iv = 0; iv < e->get_num_face_vertices(iface); iv++)
					calc_vertex_boundary_projection(e, vtx[iv]);

				// mark the edges on the boundary
				const int *edge = e->get_face_edges(iface);
				for (int ie = 0; ie < e->get_num_face_edges(iface); ie++) {
					Word_t edge_id = mesh->get_edge_id(e, edge[ie]);
					if (mesh->edges[edge_id].bnd == 0)
						EXIT("Edge #%ld should be a boundary edge.\n", edge_id);
				}
			}
		}
	}

	// update constrains recursively
	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		uc_dep(eid);
	}
}

//// BC stuff /////////////////////////////////////////////////////////////////////////////////////

static BCType default_bc_type(int marker) {
	return BC_ESSENTIAL;
}

static scalar default_bc_value_by_coord(int ess_bdy_marker, double x, double y, double z) {
	return 0.0;
}

static scalar3 &default_bc_vec_value_by_coord(int ess_bdy_marker, double x, double y, double z) {
	static scalar3 val = { 0.0, 0.0, 0.0 };
	return val;
}

void Space::set_bc_types(BCType(*bc_type_callback)(int)) {
	_F_
	if (bc_type_callback == NULL) bc_type_callback = default_bc_type;
	this->bc_type_callback = bc_type_callback;
	seq++;
}

void Space::set_essential_bc_values(scalar(*bc_value_callback_by_coord)(int, double, double, double)) {
	_F_
	if (bc_value_callback_by_coord == NULL) bc_value_callback_by_coord = default_bc_value_by_coord;
	this->bc_value_callback_by_coord = bc_value_callback_by_coord;
	seq++;
}

void Space::set_essential_bc_values(scalar3 &(*bc_vec_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z)) {
	_F_
	if (bc_vec_value_callback_by_coord == NULL) bc_vec_value_callback_by_coord = default_bc_vec_value_by_coord;
	this->bc_vec_value_callback_by_coord = bc_vec_value_callback_by_coord;
	seq++;
}

void Space::copy_callbacks(const Space *space) {
	_F_
	bc_type_callback = space->bc_type_callback;
	bc_value_callback_by_coord = space->bc_value_callback_by_coord;
	bc_vec_value_callback_by_coord = space->bc_vec_value_callback_by_coord;
}

void Space::calc_boundary_projections() {
	_F_
	FOR_ALL_ACTIVE_ELEMENTS(elm_idx, mesh) {
		Element *e = mesh->elements[elm_idx];
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(e, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::OUTER) {
				const int *vtx = e->get_face_vertices(iface);
				for (int iv = 0; iv < e->get_num_face_vertices(iface); iv++) {
					calc_vertex_boundary_projection(e, vtx[iv]);
				}

				const int *edge = e->get_face_edges(iface);
				for (int ie = 0; ie < e->get_num_face_edges(iface); ie++)
					calc_edge_boundary_projection(e, edge[ie]);

				calc_face_boundary_projection(e, iface);
			}
		}
	}
}

void Space::dump() {
	FOR_ALL_ACTIVE_ELEMENTS(eid, mesh) {
		Element *e = mesh->elements[eid];

		Word_t vtcs[Hex::NUM_VERTICES];
		e->get_vertices(vtcs);
		for (int iv = 0; iv < Hex::NUM_VERTICES; iv++) {
			vn_data[vtcs[iv]]->dump(vtcs[iv]);
		}

		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			Word_t edge = mesh->get_edge_id(e, iedge);
			en_data[edge]->dump(edge);
		}

		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			Word_t face = mesh->get_facet_id(e, iface);
			fn_data[face]->dump(face);
		}

		elm_data[eid]->dump(eid);
	}
}
