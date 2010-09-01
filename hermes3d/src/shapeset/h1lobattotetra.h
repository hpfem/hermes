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

#ifndef _SHAPESET_H1_LOBATTO_TETRA_H_
#define _SHAPESET_H1_LOBATTO_TETRA_H_

#include "../shapeset.h"
#include "../refdomain.h"
#include "tetra.h"


/// H1 Lobatto shapeset for tetrahedra
///
/// @ingroup shapesets
class H1ShapesetLobattoTetra : public Shapeset {
public:
	H1ShapesetLobattoTetra();
	virtual ~H1ShapesetLobattoTetra();

	virtual int get_vertex_index(int vertex) const {
		CHECK_VERTEX(vertex);
		return vertex_indices[vertex];
	}

	virtual int *get_edge_indices(int edge, int ori, order1_t order) {
		CHECK_EDGE(edge);
		return edge_indices[edge][ori];
	}

	virtual int *get_face_indices(int face, int ori, order2_t order) {
		return face_indices[face][ori];
	}

	virtual int *get_bubble_indices(order3_t order) {
		return bubble_indices[order.get_idx()];
	}

	virtual int get_num_edge_fns(order1_t order) const {
		return edge_count[order];
	}

	virtual int get_num_face_fns(order2_t order) const {
		return face_count[order.get_idx()];
	}

	virtual int get_num_bubble_fns(order3_t order) const {
		return bubble_count[order.get_idx()];
	}

	virtual int get_face_orientations(int face) const { return RefTetra::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefTetra::get_edge_orientations(); }

	virtual order3_t get_order(int index) const;

	virtual order3_t get_dcmp(int index) const { return order3_t(-1); }

	virtual int get_shape_type(int index) const {
		return -1;
	}

	virtual void get_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals) {
		CHECK_COMPONENT(component);
		for (int k = 0; k < np; k++)
			vals[k] = shape_table[n][component][index](pt[k].x, pt[k].y, pt[k].z);
	}

	virtual double get_value(int n, int index, double x, double y, double z, int component) {
		CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}

protected:
	shape_fn_t **shape_table[VALUE_TYPES];

	/// Indices of vertex shape functions on reference element, indexing: []
	int   *vertex_indices;
	/// Indices of edge shape functions on reference element, indexing: [edge index][ori][]
	int ***edge_indices;
	/// Indices of face shape functions on reference element, indexing: [face index][ori][]
	int ***face_indices;
	/// Indices of bubble functions on reference element, indexing: [order][]
	int  **bubble_indices;

	/// Number of edge shape functions on reference element, indexing: [order]
	int *edge_count;
	/// Number of face shape functions on reference element, indexing: [order]
	int *face_count;
	/// Number of bubble functions on reference element, indexing: [order]
	int *bubble_count;
};

#undef CHECK_VERTEX
#undef CHECK_EDGE
#undef CHECK_FACE
#undef CHECK_FACE_MODE
#undef CHECK_FACE_ORI

#endif
