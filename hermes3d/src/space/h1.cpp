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

#include "h1.h"
#include "../refmap.h"
#include "../shapeset/h1lobattohex.h"
#include "../shapeset/h1lobattotetra.h"
#include "../../../hermes_common/matrix.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/error.h"


H1Space::H1Space(Mesh* mesh, BCType (*bc_type_callback)(int), 
                 scalar (*bc_value_callback_by_coord)(int, double, double, double), Ord3 p_init, 
                 Shapeset* shapeset)
       : Space(mesh, shapeset, bc_type_callback, bc_value_callback_by_coord, p_init)
{
  _F_ 
  if (shapeset == NULL) {
    switch (p_init.type) {
      case HERMES_MODE_TET: this->shapeset = new H1ShapesetLobattoTetra; break;
      case HERMES_MODE_HEX:  this->shapeset = new H1ShapesetLobattoHex; break;
      //case HERMES_MODE_PRISM: this->shapeset = new H1ShapesetLobattoPrism; break;
      default: error("Unknown element type in H1Space::H1Space().");
    }
  }
  this->type = HERMES_H1_SPACE;

  // set uniform poly order in elements
  switch (p_init.type) {
    case HERMES_MODE_HEX: 
      if (p_init.x < 1 || p_init.y < 1 || p_init.z < 1) {
        error("P_INIT must be >= 1 in all directions in an H1 space on hexahedra.");
      }
      else this->set_uniform_order_internal(p_init);
      break;
    case HERMES_MODE_TET: 
      if (p_init.order < 1) error("P_INIT must be >= 1 in an H1 space on tetrahedra.");
      else this->set_uniform_order_internal(p_init);
      break;
    default: error("Unknown element type in H1Space::H1Space().");
  }

  // enumerate basis functions
  this->assign_dofs();
}

H1Space::~H1Space() {
  _F_
}

Space *H1Space::dup(Mesh *mesh_ext) const 
{
  _F_
  // FIXME; this only works for hexahedra.
  H1Space *space = new H1Space(mesh_ext, NULL, NULL, Ord3(1, 1, 1), this->shapeset);
  space->copy_callbacks(this);
  
  // enumerate basis functions
  space->assign_dofs();

  return space;
}

void H1Space::set_shapeset(Shapeset *shapeset)
{
  if(shapeset->get_id() < 10)
    this->shapeset = shapeset;
  else
    error("Wrong shapeset type in H1Space::set_shapeset()");
}

// ndofs ////

int H1Space::get_vertex_ndofs() 
{
  return 1;
}

int H1Space::get_edge_ndofs(Ord1 order) 
{
  return order - 1;
}

int H1Space::get_face_ndofs(Ord2 order) 
{
  switch (order.type) {
    case HERMES_MODE_TRIANGLE: return (order.order - 1) * (order.order - 2) / 2;
    case HERMES_MODE_QUAD: return (order.x - 1) * (order.y - 1);
    default: EXIT(HERMES_ERR_UNKNOWN_MODE); return -1;
  }
}

int H1Space::get_element_ndofs(Ord3 order) {
  switch (order.type) {
    case HERMES_MODE_TET: return (order.order - 1) * (order.order - 2) * (order.order - 3) / 6;
    case HERMES_MODE_HEX: return (order.x - 1) * (order.y - 1) * (order.z - 1);
    default: EXIT(HERMES_ERR_UNKNOWN_MODE); return -1;
  }
}

void H1Space::assign_dofs_internal() {
	_F_
	std::map<unsigned int, bool> init_vertices;
	std::map<Edge::Key, bool> init_edges;
	std::map<Facet::Key, bool> init_faces;

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = mesh->elements[it->first];
		  // vertex dofs
		  for (int ivtx = 0; ivtx < e->get_num_vertices(); ivtx++) {
			  unsigned int vid = e->get_vertex(ivtx);
			  VertexData *vd = vn_data[vid];
			  assert(vd != NULL);
			  if (!init_vertices[vid] && !vd->ced) {
				  assign_vertex_dofs(vid);
				  init_vertices[vid] = true;
			  }
		  }
	  }

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

void H1Space::get_element_assembly_list(Element *e, AsmList *al) {
	_F_
	al->clear();
	for (int i = 0; i < e->get_num_vertices(); i++) get_vertex_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_edges(); i++) get_edge_assembly_list(e, i, al);
	for (int i = 0; i < e->get_num_faces(); i++) get_face_assembly_list(e, i, al);
	get_bubble_assembly_list(e, al);
}

void H1Space::get_boundary_assembly_list(Element *e, int face, AsmList *al) {
	_F_
	al->clear();
	const int *face_vtcs = e->get_face_vertices(face);
	const int *face_edges = e->get_face_edges(face);
	for (int i = 0; i < e->get_num_face_vertices(face); i++) get_vertex_assembly_list(e, face_vtcs[i], al);
	for (int i = 0; i < e->get_num_face_edges(face); i++) get_edge_assembly_list(e, face_edges[i], al);
	get_face_assembly_list(e, face, al);
}


// boundary projections ////

void H1Space::calc_vertex_boundary_projection(Element *elem, int ivertex) {
	_F_
	unsigned int vtx = elem->get_vertex(ivertex);
	VertexData *vnode = vn_data[vtx];
	Vertex *v = mesh->vertices[vtx];
	if (vnode->bc_type == H3D_BC_ESSENTIAL)
		vnode->bc_proj = bc_value_callback_by_coord(vnode->marker, v->x, v->y, v->z);
}

void H1Space::calc_edge_boundary_projection(Element *elem, int iedge) {
	_F_
	Edge::Key edge_id = mesh->get_edge_id(elem, iedge);
	EdgeData *enode = en_data[edge_id];
	if (enode->bc_type != H3D_BC_ESSENTIAL) return;			// process only Dirichlet BC
	if (enode->bc_proj != NULL) return;					// projection already calculated

	int num_fns;
	if (enode->ced) {
		if(enode->edge_ncomponents <= 0)
      num_fns = 0;
    else {
		  edge_id = enode->edge_baselist[0].edge_id;
		  num_fns = en_data[edge_id]->n;
    }
	}
	else {
		num_fns = enode->n;
	}
	if (num_fns <= 0) return;

	double **proj_mat = new_matrix<double>(num_fns, num_fns);
	MEM_CHECK(proj_mat);
	scalar *proj_rhs = new scalar[num_fns];
	MEM_CHECK(proj_rhs);
	memset(proj_rhs, 0, sizeof(scalar) * num_fns);

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());
	// local edge vertex numbers
	const int *local_edge_vtx = elem->get_edge_vertices(iedge);
	// edge vertices (global indices)
	unsigned int edge_vtx[2] = { elem->get_vertex(local_edge_vtx[0]), elem->get_vertex(local_edge_vtx[1]) };

	int vtx_fn_idx[] = { shapeset->get_vertex_index(local_edge_vtx[0]), shapeset->get_vertex_index(local_edge_vtx[1]) };
	// function values at vertices
	scalar vtx_fn_coef[] = { vn_data[edge_vtx[0]]->bc_proj, vn_data[edge_vtx[1]]->bc_proj };
	int *edge_fn_idx = new int[num_fns];
	if (enode->ced && enode->edge_ncomponents > 0) {
		BaseEdgeComponent *ecomp = enode->edge_baselist + 0;
		EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

		int *indices = shapeset->get_edge_indices(iedge, ecomp->ori, cng_enode->order);
		for (int j = 0; j < cng_enode->n; j++) {
			int order = shapeset->get_order(indices[j]).get_edge_order(iedge);
			edge_fn_idx[j] = shapeset->get_constrained_edge_index(iedge, ecomp->ori, order, ecomp->part);
		}
	}
	else {
		int ori = elem->get_edge_orientation(iedge);								// edge orientation
		int *idx = shapeset->get_edge_indices(iedge, ori, enode->order);			// indices of edge functions
		for (int m = 0; m < enode->n; m++)
			edge_fn_idx[m] = idx[m];
	}

	for (int i = 0; i < num_fns; i++) {
		int iidx = edge_fn_idx[i];
		for (int j = i; j < num_fns; j++) {
			int jidx = edge_fn_idx[j];
			Ord3 order = shapeset->get_order(iidx) + shapeset->get_order(jidx);
			int edge_order = order.get_edge_order(iedge);

			int np = quad->get_edge_num_points(iedge, edge_order);
			QuadPt3D *pt = quad->get_edge_points(iedge, edge_order);
			double * vi = new double[np];
      double * vj = new double[np];
			shapeset->get_fn_values(iidx, np, pt, 0, vi);
			shapeset->get_fn_values(jidx, np, pt, 0, vj);

			double value = 0.0;
			for (int k = 0; k < np; k++)
				value += pt[k].w * vi[k] * vj[k];
			proj_mat[i][j] += value;
      delete [] vi;
      delete [] vj;
		}

		int order_rhs = quad->get_edge_max_order(iedge);
		int np = quad->get_edge_num_points(iedge, order_rhs);
		QuadPt3D *pt = quad->get_edge_points(iedge, order_rhs);

		double *edge_phys_x = ref_map.get_phys_x(np, pt);
		double *edge_phys_y = ref_map.get_phys_y(np, pt);
		double *edge_phys_z = ref_map.get_phys_z(np, pt);

		double *v0 = new double[np];
    double *v1 = new double[np];
    double *vi = new double[np];
		shapeset->get_fn_values(vtx_fn_idx[0], np, pt, 0, v0);
		shapeset->get_fn_values(vtx_fn_idx[1], np, pt, 0, v1);
		shapeset->get_fn_values(iidx, np, pt, 0, vi);

		scalar value = 0.0;
		for (int k = 0; k < np; k++) {
			scalar g = vtx_fn_coef[0] * v0[k] + vtx_fn_coef[1] * v1[k];
			value += pt[k].w * vi[k] *
				(bc_value_callback_by_coord(enode->marker, edge_phys_x[k], edge_phys_y[k], edge_phys_z[k]) - g);
		}
		proj_rhs[i] += value;

		delete [] edge_phys_x;
		delete [] edge_phys_y;
		delete [] edge_phys_z;
    delete [] v0;
    delete [] v1;
    delete [] vi;
	}

	double *chol_p = new double[num_fns];
	choldc(proj_mat, num_fns, chol_p);
	cholsl(proj_mat, num_fns, chol_p, proj_rhs, proj_rhs);
	delete [] chol_p;

	enode->bc_proj = proj_rhs;

	delete [] proj_mat;
	delete [] edge_fn_idx;
}

void H1Space::calc_face_boundary_projection(Element *elem, int iface) {
	_F_
	Facet::Key facet_idx = mesh->get_facet_id(elem, iface);
	FaceData *fnode = fn_data[facet_idx];

	if (fnode->bc_type != H3D_BC_ESSENTIAL) return;
	if (fnode->bc_proj != NULL) return;
	if (fnode->n <= 0) return;

	double **proj_mat = new_matrix<double>(fnode->n, fnode->n);
	MEM_CHECK(proj_mat);
	scalar *proj_rhs = new scalar[fnode->n];
	MEM_CHECK(proj_rhs);
	memset(proj_rhs, 0, sizeof(double) * fnode->n);

	RefMap ref_map(mesh);
	ref_map.set_active_element(elem);

	Quad3D *quad = get_quadrature(elem->get_mode());

	const int *local_face_vertex = elem->get_face_vertices(iface);
	const int *local_face_edge = elem->get_face_edges(iface);

	// get total number of vertex + edge functions
	int num_fns = elem->get_num_face_vertices(iface);
	for (int edge = 0; edge < elem->get_num_face_edges(iface); edge++) {
    Edge::Key edge_idx = mesh->get_edge_id(elem, local_face_edge[edge]);
		EdgeData *enode = en_data[edge_idx];
		if (enode->ced && enode->edge_ncomponents > 0 && enode->edge_baselist != NULL) {
			assert(enode->edge_ncomponents > 0);
      Edge::Key eid = enode->edge_baselist[0].edge_id;
			num_fns += en_data[eid]->n;
		}
		else
			num_fns += enode->n;
	}

	scalar *coef = new scalar[num_fns];
	MEM_CHECK(coef);
	int *fn_idx = new int[num_fns];
	MEM_CHECK(fn_idx);

	int m = 0;
	// vertex projection coefficients
	for (int vtx = 0; vtx < elem->get_num_face_vertices(iface); vtx++, m++) {
		VertexData *vnode = vn_data[elem->get_vertex(local_face_vertex[vtx])];
		coef[m] = vnode->bc_proj;
		fn_idx[m] = shapeset->get_vertex_index(local_face_vertex[vtx]);
	}
	// edge projection coefficients
	for (int edge = 0; edge < elem->get_num_face_edges(iface); edge++) {
		Edge::Key edge_idx = mesh->get_edge_id(elem, local_face_edge[edge]);
		EdgeData *enode = en_data[edge_idx];

		if (enode->ced && enode->edge_ncomponents > 0 && enode->edge_baselist != NULL) {
				BaseEdgeComponent *ecomp = enode->edge_baselist + 0;
				EdgeData *cng_enode = en_data[ecomp->edge_id]; 						// constraining edge node

				if (cng_enode->n > 0) {
					int *indices = shapeset->get_edge_indices(local_face_edge[edge], ecomp->ori, cng_enode->order);
					for (int j = 0; j < cng_enode->n; j++, m++) {
						int order = shapeset->get_order(indices[j]).get_edge_order(local_face_edge[edge]);
						fn_idx[m] = shapeset->get_constrained_edge_index(local_face_edge[edge], ecomp->ori, order, ecomp->part);
						assert(cng_enode->bc_proj != NULL);
						coef[m] = cng_enode->bc_proj[j];
					}
				}
		}
		else {
			if (enode->n > 0) {
				int edge_ori = elem->get_edge_orientation(local_face_edge[edge]);
				int *edge_fn_idx = shapeset->get_edge_indices(local_face_edge[edge], edge_ori, enode->order);

				for (int i = 0; i < enode->n; i++, m++) {
					assert(enode->bc_proj != NULL);
					coef[m] = enode->bc_proj[i];
					fn_idx[m] = edge_fn_idx[i];
				}
			}
		}
	}

	// do it //
	int face_ori = elem->get_face_orientation(iface);
	int *face_fn_idx = shapeset->get_face_indices(iface, face_ori, fnode->order);

	for (int i = 0; i < fnode->n; i++) {
		int iidx = face_fn_idx[i];
		for (int j = i; j < fnode->n; j++) {
			int jidx = face_fn_idx[j];

			Ord3 order = shapeset->get_order(iidx) + shapeset->get_order(jidx);
			Ord2 face_order = order.get_face_order(iface);

			int np = quad->get_face_num_points(iface, face_order);
			QuadPt3D *pt = quad->get_face_points(iface, face_order);
			double *vi = new double[np];
			double *vj = new double[np];
			shapeset->get_fn_values(iidx, np, pt, 0, vi);
			shapeset->get_fn_values(jidx, np, pt, 0, vj);

			double value = 0.0;
			for (int k = 0; k < np; k++)
				value += pt[k].w * vi[k] * vj[k];
			proj_mat[i][j] += value;

      delete [] vi;
      delete [] vj;
		}

		Ord2 order_rhs = quad->get_face_max_order(iface);
		int np = quad->get_face_num_points(iface, order_rhs);
		QuadPt3D *pt = quad->get_face_points(iface, order_rhs);

		double *face_phys_x = ref_map.get_phys_x(np, pt);
		double *face_phys_y = ref_map.get_phys_y(np, pt);
		double *face_phys_z = ref_map.get_phys_z(np, pt);

		scalar value = 0.0;
		scalar *g = new scalar[np];		// lin. combination of vertex + edge functions
		memset(g, 0, np * sizeof(double));
		for (int l = 0; l < num_fns; l++) {
			double *fn = new double[np];
			shapeset->get_fn_values(fn_idx[l], np, pt, 0, fn);
			for (int k = 0; k < np; k++)
				g[k] += coef[l] * fn[k];
      delete [] fn;
		}

		double *vi = new double[np];
		shapeset->get_fn_values(iidx, np, pt, 0, vi);
		for (int k = 0; k < np; k++) {
			value += pt[k].w * vi[k] *
				(bc_value_callback_by_coord(fnode->marker, face_phys_x[k], face_phys_y[k], face_phys_z[k]) - g[k]);
		}

		proj_rhs[i] += value;

		delete [] face_phys_x;
		delete [] face_phys_y;
		delete [] face_phys_z;
    delete [] g;
    delete [] vi;
	}

	// solve the system using a precalculated Cholesky decomposed projection matrix
	double *chol_p = new double[fnode->n];
	choldc(proj_mat, fnode->n, chol_p);
	cholsl(proj_mat, fnode->n, chol_p, proj_rhs, proj_rhs);
	delete [] chol_p;

	fnode->bc_proj = proj_rhs;

	delete [] fn_idx;
	delete [] coef;
	delete [] proj_mat;
}
