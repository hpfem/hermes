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

#ifndef _ORDER_H_
#define _ORDER_H_

#include "h3d_common.h"
#include "../../hermes_common/error.h"


// maximal order of quadratures for 1D
#define H3D_MAX_QUAD_ORDER		         24
// maximal order of quadratures for triangle
#define H3D_MAX_QUAD_ORDER_TRI			 20
// maximal order of quadratures for tetra
#define H3D_MAX_QUAD_ORDER_TETRA		 20


inline int limit_tri_ord(int x) { return (x > H3D_MAX_QUAD_ORDER_TRI) ? H3D_MAX_QUAD_ORDER_TRI : x; }
inline int limit_tet_ord(int x) { return (x > H3D_MAX_QUAD_ORDER_TETRA) ? H3D_MAX_QUAD_ORDER_TETRA : x; }
inline int limit_quad_ord(int x) { return (x > H3D_MAX_QUAD_ORDER) ? H3D_MAX_QUAD_ORDER : x; }

// orders

// 1D polynomial order
typedef
	int Ord1;


// 2D polynomial order
struct Ord2 {
	Ord2() { type = 3; order = 31; }
	Ord2(int order) { type = MODE_TRIANGLE; this->order = order; }
	Ord2(int x, int y) { type = MODE_QUAD; this->x = x; this->y = y; }

	unsigned type:2;		// EMode2D
	union {
		// tri
		struct {
			unsigned order:5;
		};
		// quad
		struct {
			// directional orders
			unsigned x:5;
			unsigned y:5;
		};
	};

	bool invalid() { return (type == 3); }

	// Operators

	Ord2 operator+(const Ord2 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	return Ord2(limit_tri_ord(this->order + o.order));
			case MODE_QUAD:
				return Ord2(limit_quad_ord(this->x + o.x),
				                limit_quad_ord(this->y + o.y));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord2(-1);
	}

	Ord2 operator+=(const Ord2 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	this->order = limit_tri_ord(this->order + o.order); break;
			case MODE_QUAD:
				this->x = limit_quad_ord(this->x + o.x);
				this->y = limit_quad_ord(this->y + o.y);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	Ord2 operator*(const int c) {
		switch (type) {
			case MODE_TRIANGLE: return Ord2(limit_tri_ord(c * this->order));
			case MODE_QUAD:
				return Ord2(limit_quad_ord(c * this->x),
				                limit_quad_ord(c * this->y));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord2(-1);
	}

	Ord2 operator*(const Ord2 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	return Ord2(limit_tri_ord(this->order * o.order));
			case MODE_QUAD:
				return Ord2(limit_quad_ord(this->x * o.x),
				                limit_quad_ord(this->y * o.y));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord2(-1);
	}

	Ord2 operator*=(const int c) {
		switch (type) {
			case MODE_TRIANGLE:	this->order = limit_tri_ord(this->order * c); break;
			case MODE_QUAD:
				this->x = limit_quad_ord(this->x * c);
				this->y = limit_quad_ord(this->y * c);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	Ord2 operator*=(const Ord2 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TRIANGLE:	this->order = limit_tri_ord(this->order * o.order); break;
			case MODE_QUAD:
				this->x = limit_quad_ord(this->x * o.x);
				this->y = limit_quad_ord(this->x * o.y);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	// relation operators
	bool operator==(const Ord2 &o) {
		if (this->type != o.type) return false;
		switch (this->type) {
			case MODE_TRIANGLE: return this->order == o.order;
			case MODE_QUAD: return (this->x == o.x) && (this->y == o.y);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	bool operator!=(const Ord2 &o) {
		if (this->type != o.type) return true;
		switch (this->type) {
			case MODE_TRIANGLE: return this->order != o.order;
			case MODE_QUAD: return (this->x != o.x) || (this->y != o.y);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	const char *str() {
		static char s[64];
		switch (type) {
			case MODE_TRIANGLE: sprintf(s, "(%d)", this->order); break;
			case MODE_QUAD: sprintf(s, "(%d, %d)", this->x, this->y); break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return s;
	}

	int get_idx() const {
		switch (type) {
			case MODE_TRIANGLE: return (this->type << 10) | this->order;
			case MODE_QUAD: return (((this->type << 5) | this->y) << 5) | this->x;
			default: assert(false); EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return -1;
	}

	static Ord2 from_int(int o) {
		int type = (o >> 10) & 0x03;
		switch (type) {
			case MODE_TRIANGLE: return Ord2(o & 0x1F); break;
			case MODE_QUAD: return Ord2(o & 0x1F, (o >> 5) & 0x1F); break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord2(-1);
	}
};

inline Ord2 max(Ord2 a, Ord2 b) {
	assert(a.type == b.type);
	switch (a.type) {
		case MODE_TRIANGLE: return Ord2(std::max(a.order, b.order));
		case MODE_QUAD: return Ord2(std::max(a.x, b.x), std::max(a.y, b.y));
		default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
	}
	return Ord2(-1);
}


// 3D polynomial order
//
// all 1s mean invalid (not set) - see default ctor
//
struct Ord3 {
	Ord3() { invalid(); }
	Ord3(int order) { type = MODE_TETRAHEDRON; this->order = limit_tet_ord(order); }
	Ord3(int x, int y, int z) {
		type = MODE_HEXAHEDRON;
		this->x = limit_quad_ord(x);
		this->y = limit_quad_ord(y);
		this->z = limit_quad_ord(z);
	}

	signed type:3;		// EMode3D
	union {
		// tetra
		struct {
			unsigned order: 15;
		};
		// hex
		struct {
			// directional orders
			unsigned x:5;
			unsigned y:5;
			unsigned z:5;
		};
	};

	bool is_invalid() const { return (type == -1); }
	void invalid() { type = 7; }

	int get_ord() const {
		switch (type) {
			case MODE_TETRAHEDRON:	return order;
			case MODE_HEXAHEDRON: return std::max(std::max(this->x, this->y), this->z);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return -1;
	}

	// Operators

	Ord3 operator+(const Ord3 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON: return Ord3(limit_tet_ord(this->order + o.order));
			case MODE_HEXAHEDRON:
				return Ord3(limit_quad_ord(this->x + o.x),
				                limit_quad_ord(this->y + o.y),
			                    limit_quad_ord(this->z + o.z));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord3(-1);
	}

	Ord3 operator+=(const Ord3 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON: this->order = limit_tet_ord(this->order + o.order); break;
			case MODE_HEXAHEDRON:
				this->x = limit_quad_ord(this->x + o.x);
				this->y = limit_quad_ord(this->y + o.y);
				this->z = limit_quad_ord(this->z + o.z);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	Ord3 operator*(const int c) {
		switch (type) {
			case MODE_TETRAHEDRON: return Ord3(limit_tet_ord(c * this->order));
			case MODE_HEXAHEDRON:
				return Ord3(limit_quad_ord(c * this->x),
				                limit_quad_ord(c * this->y),
				                limit_quad_ord(c * this->z));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord3(-1);
	}

	Ord3 operator*(const Ord3 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON:	return Ord3(limit_tet_ord(this->order * o.order));
			case MODE_HEXAHEDRON:
				return Ord3(limit_quad_ord(this->x * o.x),
				                limit_quad_ord(this->y * o.y),
				                limit_quad_ord(this->z * o.z));
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord3(-1);
	}

	Ord3 operator*=(const int c) {
		switch (type) {
			case MODE_TETRAHEDRON:	this->order = limit_tet_ord(this->order * c); break;
			case MODE_HEXAHEDRON:
				this->x = limit_quad_ord(this->x * c);
				this->y = limit_quad_ord(this->y * c);
				this->z = limit_quad_ord(this->z * c);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	Ord3 operator*=(const Ord3 &o) {
		assert(type == o.type);
		switch (type) {
			case MODE_TETRAHEDRON: this->order = limit_tet_ord(this->order * o.order); break;
			case MODE_HEXAHEDRON:
				this->x = limit_quad_ord(this->x * o.x);
				this->y = limit_quad_ord(this->y * o.y);
				this->z = limit_quad_ord(this->z * o.z);
				break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return *this;
	}

	// relation operators
	bool operator==(const Ord3 &o) {
		if (this->type != o.type) return false;
		switch (this->type) {
			case MODE_TETRAHEDRON: return this->order == o.order;
			case MODE_HEXAHEDRON: return (this->x == o.x) && (this->y == o.y) && (this->z == o.z);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	bool operator!=(const Ord3 &o) {
		if (this->type != o.type) return true;
		switch (this->type) {
			case MODE_TETRAHEDRON: return this->order != o.order;
			case MODE_HEXAHEDRON: return (this->x != o.x) || (this->y != o.y) || (this->z != o.z);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return false;
	}

	const char *str() const {
		static char s[64];
		switch (type) {
			case MODE_TETRAHEDRON: sprintf(s, "(%d)", this->order); break;
			case MODE_HEXAHEDRON: sprintf(s, "(%d, %d, %d)", this->x, this->y, this->z); break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return s;
	}

	int get_idx() const {
		assert(!is_invalid());
		switch (type) {
			case MODE_TETRAHEDRON: return ((this->type) << 15) | this->order;
			case MODE_HEXAHEDRON:
				return (((((this->type << 5) | this->z) << 5) | this->y) << 5) | this->x;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return -1;
	}

	static Ord3 from_int(int o) {
		int type = (o >> 15) & 0x07;
		switch (type) {
			case MODE_TETRAHEDRON: return Ord3(o & 0x7FFF);
			case MODE_HEXAHEDRON: return Ord3(o & 0x1F, (o >> 5) & 0x1F, (o >> 10) & 0x1F);
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}
		return Ord3(-1);
	}

	Ord1 get_edge_order(int edge) const {
		switch (type) {
			case MODE_TETRAHEDRON: return this->order;
			case MODE_HEXAHEDRON:
				if ((edge == 0) || (edge == 2) || (edge == 10) || (edge == 8)) return this->x;
				else if((edge == 1) || (edge == 3) || (edge == 11) || (edge == 9)) return this->y;
				else if((edge == 4) || (edge == 5) || (edge == 6) || (edge == 7)) return this->z;
				else EXIT(H3D_ERR_EDGE_INDEX_OUT_OF_RANGE);
			default:
				EXIT(HERMES_ERR_UNKNOWN_MODE);
				break;
		}
		return -1;
	}

	Ord2 get_face_order(int face) const {
		switch (type) {
			case MODE_TETRAHEDRON: return this->order;
			case MODE_HEXAHEDRON:
				if ((face == 0) || (face == 1)) return Ord2(this->y, this->z);
				else if ((face == 2) || (face == 3)) return Ord2(this->x, this->z);
				else if ((face == 4) || (face == 5)) return Ord2(this->x, this->y);
				else EXIT(H3D_ERR_FACE_INDEX_OUT_OF_RANGE);
			default:
				EXIT(HERMES_ERR_UNKNOWN_MODE);
				break;
		}
		return Ord2(-1);
	}

	void limit() {
#ifdef DEBUG_ORDER
		set_maximal();
#endif
	}

	void set_maximal() {
		switch (type) {
			case MODE_TETRAHEDRON: this->order = H3D_MAX_QUAD_ORDER_TETRA; break;

			case MODE_HEXAHEDRON:
				this->x = H3D_MAX_QUAD_ORDER;
				this->y = H3D_MAX_QUAD_ORDER;
				this->z = H3D_MAX_QUAD_ORDER;
				break;

			default:
				EXIT(HERMES_ERR_UNKNOWN_MODE);
				break;
		}
	}
};

inline Ord3 operator*(const int c, const Ord3 &a) {
	switch (a.type) {
		case MODE_TETRAHEDRON:	return Ord3(limit_tet_ord(c * a.order));
		case MODE_HEXAHEDRON:
			return Ord3(limit_quad_ord(c * a.x),
			                limit_quad_ord(c * a.y),
			                limit_quad_ord(c * a.z));
		default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
	}
	return Ord3(-1);
}

inline Ord3 max(Ord3 a, Ord3 b) {
	if (a.type == 7) return b;
	if (b.type == 7) return a;

	assert(a.type == b.type);
	switch (a.type) {
		case MODE_TETRAHEDRON: return Ord3(std::max(a.order, b.order));
		case MODE_HEXAHEDRON:
			return Ord3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
		default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
	}
	return Ord3(-1);
}

inline Ord3 turn_hex_face_order(int iface, Ord3 ord) {
	int o1 = ord.x;
	int o2 = ord.y;
	int o3 = ord.z;
	if (iface == 0 || iface == 1) std::swap(o2, o3);
	else if (iface == 2 || iface == 3) std::swap(o1, o3);
	else if (iface == 4 || iface == 5) std::swap(o1, o2);
	else assert(false);
	return Ord3(o1, o2, o3);
}

#endif
