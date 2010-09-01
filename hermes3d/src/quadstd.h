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

#ifndef _QUAD_STD_H_
#define _QUAD_STD_H_

#include "quad.h"

/// Numeric quadrature for 1D (Gauss points)
///
/// @ingroup quadrature
class QuadStd1D : public Quad1D {
public:
	QuadStd1D();
	~QuadStd1D();
};


/// Numeric quadrature for 2D triangle
///
/// @ingroup quadrature
class QuadStdTri : public Quad2D {
public:
	QuadStdTri();
	~QuadStdTri();
};


/// Numerical quadrature for 3D hexahedron
///
/// @ingroup quadrature
class QuadStdHex : public Quad3D {
public:
	QuadStdHex();
	~QuadStdHex();

	virtual QuadPt3D *get_points(const order3_t &order) {
		CHECK_MODE;
		if (!tables.exists(order.get_idx())) calc_table(order);
		return tables[order.get_idx()];
	}

	virtual QuadPt3D *get_face_points(int face, const order2_t &order) {
		if (!face_tables[face].exists(order.get_idx())) calc_face_table(face, order);
		return face_tables[face][order.get_idx()];
	}

protected:
	void calc_table(const order3_t &order);
	void calc_face_table(int face, const order2_t &order);
	///
	order3_t lower_order_same_accuracy(const order3_t &ord);
};


/// Numerical quadrature for 3D tetrahedron
///
/// @ingroup quadrature
class QuadStdTetra : public Quad3D {
public:
	QuadStdTetra();
	~QuadStdTetra();
};


/// Numerical quadrature for 3D prisms
///
/// @ingroup quadrature
class QuadStdPrism : public Quad3D {
public:
	QuadStdPrism();
	~QuadStdPrism();
};

#endif
