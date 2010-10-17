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

#ifndef _SPACE_H1_H_
#define _SPACE_H1_H_

#include "space.h"

/// H1 space
///
/// @ingroup spaces
class H3D_API H1Space : public Space {
public:
	H1Space(Mesh *mesh, BCType (*bc_type_callback)(int) = NULL, 
		scalar (*bc_value_callback_by_coord)(int, double, double, double) = NULL, Ord3 p_init = Ord3(1,1,1),
                Shapeset* shapeset = NULL);
	virtual ~H1Space();

	virtual Space *dup(Mesh *mesh) const;

  virtual void set_shapeset(Shapeset* shapeset);

   	virtual void get_element_assembly_list(Element *e, AsmList *al);
	virtual void get_boundary_assembly_list(Element *e, int face, AsmList *al);

protected:
	virtual int get_vertex_ndofs();
	virtual int get_edge_ndofs(Ord1 order);
	virtual int get_face_ndofs(Ord2 order);
	virtual int get_element_ndofs(Ord3 order);

	virtual void assign_dofs_internal();

	virtual void calc_vertex_boundary_projection(Element *elem, int ivertex);
	virtual void calc_edge_boundary_projection(Element *elem, int iedge);
	virtual void calc_face_boundary_projection(Element *elem, int iface);
};

#endif
