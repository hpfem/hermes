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

#ifndef _QUAD_H_
#define _QUAD_H_

#include "common.h"
#include "order.h"
#include <common/error.h>
#include <common/array.h>

/// @defgroup quadratures Numerical quadratures
///
/// There is a class called QuadXD (where X = 1, 2, 3) that represents numerical
/// quadratures in X dimensions. Each class holds quadratures points for every order.
/// The order of integration is bounded by max_order. The class provides interface
/// for getting quadratures points, numbers of points and maximal order of integration.
/// Numerical quadratures of higher dimensions (2, 3) also provides integration points
/// for surface integrals. This simplifies evaluating of these integrals on elements
/// (one do NOT have to take care about transformations, it is done by the QuadXXX
/// class).
///
/// There are prepared quadratures for standard domains (line, triangle, quad, hex,
/// tetra and prism) which are described in Pavel's Gray book (see quadstd.h). You
/// can create your own quadratures by deriving a class from QuadXD.
///
///
///
/// NOTES:
/// - Quad is short name for quadratures, but can be mismatched with Quadrilateral
/// - QuadPtXD - quadratures point (find better name)
///


/// Quadrature point in 1D
///
/// @ingroup quadratures
struct QuadPt1D {
	double x;		// x-coordinate
	double w;		// weight

	QuadPt1D() { }		// default c-tor
	QuadPt1D(double x, double w) {
		this->x = x;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else EXIT("Index out of bounds");
	}
};

/// Quadrature point in 2D
///
/// @ingroup quadratures
struct QuadPt2D {
	double x, y;		// x and y-coordinate
	double w;			// weight

	QuadPt2D() { }		// default c-tor
	QuadPt2D(double x, double y, double w) {
		this->x = x;
		this->y = y;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else if (idx == 1) return this->y;
		else EXIT("Index out of bounds");
	}
};

/// Quadrature point in 3D
///
/// @ingroup quadratures
struct QuadPt3D {
	double x, y, z;		// x, y and z-coordinate
	double w;			// weight

	QuadPt3D() { } 		// default c-tor
	QuadPt3D(double x, double y, double z, double w) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else if (idx == 1) return this->y;
		else if (idx == 2) return this->z;
		else EXIT("Index out of bounds");
	}
};


/// Numerical quadrature in 1D
///
/// @ingroup quadratures
class Quad1D {
public:
	QuadPt1D *get_points(int order) const { return tables[order]; }
	inline int get_num_points(int order) const { return np[order]; };

	int get_max_order() const { return max_order; }

protected:
	/// maximal order for integration
	int max_order;
	/// number of integration points
	/// indexing: [order]
	int *np;
	/// tables with integration points
	/// indexing: [order][point no.]
	QuadPt1D **tables;
};


///
/// 2D quadratures
///

/// Numerical quadratures in 2D
///
/// @ingroup quadratures
class Quad2D {
public:
	QuadPt2D *get_points(int order) const { return tables[order]; }
	inline int get_num_points(int order) const { return np[order]; };

	QuadPt2D *get_edge_points(int edge, int order) const { return edge_tables[edge][order]; }

	int get_max_order() const { return max_order; }

	EMode2D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_TRIANGLE, MODE_QUAD)
	EMode2D mode;
	/// maximal order for integration (interpretation depened on the mode)
	int max_order;
	/// number of integration points
	/// indexing: [order]
	int *np;
	/// tables with integration points
	/// indexing: [order][point no.]
	QuadPt2D **tables;
	/// tables with integration points for edges (?)
	/// indexing: [edge][order][point no.]
	QuadPt2D ***edge_tables;
};

//
// 3D quadratures
//

#define CHECK_MODE assert(order.type == mode)

/// Numerical quadratures in 3D
///
/// @ingroup quadratures
class Quad3D {
public:
	virtual ~Quad3D() { }

	virtual QuadPt3D *get_points(const order3_t &order) { CHECK_MODE; return tables[order.get_idx()]; }
	virtual int get_num_points(const order3_t &order) { CHECK_MODE; return np[order.get_idx()]; }

	virtual QuadPt3D *get_edge_points(int edge, const order1_t &order) { return edge_tables[edge][order]; }
	int get_edge_num_points(int edge, const order1_t &order) const { return np_edge[order]; }

	virtual QuadPt3D *get_face_points(int face, const order2_t &order) { return face_tables[face][order.get_idx()]; }
	int get_face_num_points(int face, const order2_t &order) const { return np_face[order.get_idx()]; }

	virtual QuadPt3D *get_vertex_points() { return vertex_table; }
	int get_vertex_num_points() const { return np_vertex; }

	order1_t get_edge_max_order(int edge) const { return max_edge_order; }
	order2_t get_face_max_order(int face) const { return max_face_order; }
	order3_t get_max_order() const { return max_order; }

	EMode3D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM)
	EMode3D mode;
	/// maximal order for integration (interpretation depends on the mode)
	order1_t max_edge_order;
	order2_t max_face_order;
	order3_t max_order;

	Array<QuadPt3D *> tables;
	Array<QuadPt3D *> *edge_tables;
	Array<QuadPt3D *> *face_tables;
	QuadPt3D *vertex_table;
	Array<int> np;
	Array<int> np_edge;
	Array<int> np_face;
	int np_vertex;
};


// interface for getting quadratures - library wide ////////////////////////////////////////////////

Quad1D *get_quadrature_1d();
Quad3D *get_quadrature(EMode3D mode);


#ifndef DEBUG_ORDER
	#define LIMIT_TRI_ORDER(o)
	#define LIMIT_QUAD_ORDER(o)
#else
	#define LIMIT_TRI_ORDER(o) 							o = H3D_MAX_QUAD_ORDER_TRI;
	#define LIMIT_QUAD_ORDER(o) 						o = MAKE_QUAD_ORDER(H3D_MAX_QUAD_ORDER, H3D_MAX_QUAD_ORDER);
#endif

#endif
