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

#ifndef _QUAD_CHEB_H_
#define _QUAD_CHEB_H_

#include "quad.h"

/// QuadChebTetra is a special "quadrature" consisting of product Chebyshev
/// points on the reference brick. It is used for expressing
/// the solution on an element as a linear combination of monomials.
///
class HERMES_API QuadChebTetra : public Quad3D {
public:
	QuadChebTetra();
	~QuadChebTetra();
};


/// QuadChebHex is a special "quadrature" consisting of product Chebyshev
/// points on the reference brick. It is used for expressing
/// the solution on an element as a linear combination of monomials.
///
class HERMES_API QuadChebHex : public Quad3D {
public:
	QuadChebHex();
	~QuadChebHex();

	virtual QuadPt3D *get_points(const Ord3 &order) {
    if (tables->find(order.get_idx()) == tables->end()) 
      calc_table(order);
		return (*tables)[order.get_idx()];
	}

protected:
	void calc_table(const Ord3 &order);
};

#endif
