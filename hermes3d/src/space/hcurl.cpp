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

#include "hcurl.h"
#include "../shapeset/hcurllobattohex.h"
#include "../refmap.h"
#include "../refdomain.h"
#include "../../../hermes_common/matrix.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/error.h"

HcurlSpace::HcurlSpace(Mesh* mesh, BCType (*bc_type_callback)(int), 
                 scalar (*bc_value_callback_by_coord)(int, double, double, double), Ord3 p_init, 
                 Shapeset* shapeset) 
          : Space(mesh, shapeset, bc_type_callback, bc_value_callback_by_coord, p_init)
{
  _F_
  // FIXME: this will fail if the mesh contains tetrahedra. 
  if (shapeset == NULL) this->shapeset = new HcurlShapesetLobattoHex;
  this->type = HERMES_HCURL_SPACE;

  // set uniform poly order in elements
  this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::~HcurlSpace() {
	_F_
}

Space *HcurlSpace::dup(Mesh *mesh) const {
	_F_
	  HcurlSpace *space = new HcurlSpace(mesh, NULL, NULL, Ord3(-1,-1,-1), shapeset);
	space->copy_callbacks(this);
	return space;
}

void HcurlSpace::set_shapeset(Shapeset *shapeset)
{
  if(shapeset->get_id() < 20 && shapeset->get_id() > 9)
    this->shapeset = shapeset;
  else
    error("Wrong shapeset type in HcurlSpace::set_shapeset()");
}

// ndofs ////

int HcurlSpace::get_vertex_ndofs() {
	return 0;
}

int HcurlSpace::get_edge_ndofs(Ord1 order) {
	return order + 1;
}

int HcurlSpace::get_face_ndofs(Ord2 order) {
	switch (order.type) {
		case HERMES_MODE_QUAD: return (order.x + 1) * order.y + order.x * (order.y + 1);
		case HERMES_MODE_TRIANGLE: EXIT(HERMES_ERR_NOT_IMPLEMENTED); return -1;
		default: EXIT(HERMES_ERR_UNKNOWN_MODE); return -1;
	}
}

int HcurlSpace::get_element_ndofs(Ord3 order) {
	switch (order.type) {
		case HERMES_MODE_HEX: return (order.x + 1) * order.y * order.z + order.x * (order.y + 1) * order.z + order.x * order.y * (order.z + 1);
		case HERMES_MODE_TET: EXIT(HERMES_ERR_NOT_IMPLEMENTED); return -1;
		default: EXIT(HERMES_ERR_UNKNOWN_MODE); return -1;
	}
}

//

void HcurlSpace::assign_dofs_internal() {
	_F_
	std::map<Edge::Key, bool> init_edges;
	std::map<Facet::Key, bool> init_faces;

	// edge dofs
  for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = mesh->elements[it->first];
		  // edge dofs
		  for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
			  Edge::Key eid = mesh->get_edge_id(e, iedge);
			  EdgeData *ed = en_data[eid];
			  assert(ed != NULL);
			  if (!init_edges[eid] && !ed->ced) {
				  assign_edge_dofs(eid);
				  init_edges[eid] = true;
			  }
		  }
	  }

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = mesh->elements[it->first];
		// face dofs
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Facet::Key fid = mesh->get_facet_id(e, iface);
			FaceData *fd = fn_data[fid];
			assert(fd != NULL);
			if (!init_faces[fid] && !fd->ced) {
				assign_face_dofs(fid);
				init_faces[fid] = true;
			}
		}
	}

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active)
		  assign_bubble_dofs(it->first);
}

// assembly lists ////

void HcurlSpace::get_element_assembly_list(Element *e, AsmList *al) {
	_F_
	al->clear();
	for (int i = 0; i < e->get_num_edges(); i++) get_edge_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_faces(); i++) get_face_assembly_list(e, i, al);
	get_bubble_assembly_list(e, al);
}

void HcurlSpace::get_boundary_assembly_list(Element *e, int face, AsmList *al) {
	_F_
	al->clear();
	const int *face_edges = e->get_face_edges(face);
	for (int i = 0; i < e->get_num_face_edges(face); i++) get_edge_assembly_list(e, face_edges[i], al);
	get_face_assembly_list(e, face, al);
}

// boundary projections
// We allow only zero BC (see hcurl.h), so we only check whether this is true
// and fill projection with zeros

void HcurlSpace::calc_vertex_boundary_projection(Element *elem, int ivertex) {
	_F_
	unsigned int vtx = elem->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	Vertex *v = mesh->vertices[vtx];
	if (vnode->bc_type == H3D_BC_ESSENTIAL) {
		// FIXME: use bc_vec_value_callback_by_coord
		if ((vnode->bc_proj = bc_value_callback_by_coord(vnode->marker, v->x, v->y, v->z)) != 0.0)
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);  //projection of nonzero bc not implemented, see comment in .h
	}
}

void HcurlSpace::calc_edge_boundary_projection(Element *elem, int iedge) {
	_F_
  Edge::Key edge = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge];
	if (enode->bc_type != H3D_BC_ESSENTIAL) return;			// process only Dirichlet BC
	if (enode->bc_proj != NULL) return;					// projection already calculated

	int num_fns;
	if (enode->ced) {
		assert(enode->edge_ncomponents > 0);
    Edge::Key edge_id = enode->edge_baselist[0].edge_id;
		num_fns = en_data[edge_id]->n;
	}
	else
		num_fns = enode->n;
	if (num_fns <= 0) return;

	scalar *proj_rhs = new scalar[num_fns];
	MEM_CHECK(proj_rhs);
	for (int i = 0; i < num_fns; i++) proj_rhs[i] = 0.0;

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	Ord1 order_rhs = quad->get_edge_max_order(iedge);
	int np = quad->get_edge_num_points(iedge, order_rhs);
	QuadPt3D *pt = quad->get_edge_points(iedge, order_rhs);

	double *edge_phys_x = ref_map.get_phys_x(np, pt);
	double *edge_phys_y = ref_map.get_phys_y(np, pt);
	double *edge_phys_z = ref_map.get_phys_z(np, pt);

	for (int k = 0; k < np; k++) {
		// FIXME: use bc_vec_value_callback_by_coord
		if (bc_value_callback_by_coord(enode->marker, edge_phys_x[k], edge_phys_y[k], edge_phys_z[k]) != 0.)
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);  // projection of nonzero BC not implemented, see comment in .h
	}

	delete [] edge_phys_x;
	delete [] edge_phys_y;
	delete [] edge_phys_z;

	// save vector of zeros as a projection
	enode->bc_proj = proj_rhs;
}

void HcurlSpace::calc_face_boundary_projection(Element *elem, int iface) {
	_F_
	Facet::Key facet_idx = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[facet_idx];

	if (fnode->bc_type != H3D_BC_ESSENTIAL) return;
	if (fnode->bc_proj != NULL) return;
	if (fnode->n <= 0) return;

	scalar *proj_rhs = new scalar[fnode->n];
	MEM_CHECK(proj_rhs);
	for (int i = 0; i < fnode->n; i++) proj_rhs[i] = 0.0;

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	Ord2 order_rhs = quad->get_face_max_order(iface);
	int np = quad->get_face_num_points(iface, order_rhs);
	QuadPt3D *pt = quad->get_face_points(iface, order_rhs);

	double *face_phys_x = ref_map.get_phys_x(np, pt);
	double *face_phys_y = ref_map.get_phys_y(np, pt);
	double *face_phys_z = ref_map.get_phys_z(np, pt);

	for (int k = 0; k < quad->get_face_num_points(iface, order_rhs); k++) {
		// FIXME: use bc_vec_value_callback_by_coord
		if (bc_value_callback_by_coord(fnode->marker, face_phys_x[k], face_phys_y[k], face_phys_z[k]) != 0.)
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);  // projection of nonzero BC not implemented, see comment in .h
	}

	delete [] face_phys_x;
	delete [] face_phys_y;
	delete [] face_phys_z;

	// save vector of zeros as a projection
	fnode->bc_proj = proj_rhs;
}
