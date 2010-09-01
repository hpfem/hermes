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

#ifndef _FORMS_H_
#define _FORMS_H_

#include "common.h"
#include "quad.h"
#include "function.h"
#include "solution.h"
#include "refmap.h"

// TODO: better name?
#define FORM_CB(a)	a<double, scalar>, a<ord_t, ord_t>

// Base type for orders of functions
//
// We defined a special arithmetics with this type to be able to analyze forms
// and determine the necessary integration order.  This works for forms, but it also
// works for user-defined functions.
class ord_t {
public:
	ord_t() { order = 0; }
	ord_t(scalar d) { order = 0; }
	ord_t(int o) { order = o; }

	int get_order() const { return order; }
	// lowest common max order is tetrahedral
	// TODO: different max order for hexes and tets
	int get_max_order() const { return H3D_MAX_QUAD_ORDER_TETRA; }

	ord_t operator+(const ord_t &o) { return ord_t(std::max(this->order, o.order)); }
	ord_t operator-(const ord_t &o) { return ord_t(std::max(this->order, o.order)); }
	ord_t operator-(double d) { return *this; }
	ord_t operator*(const ord_t &o) { return ord_t(this->order + o.order); }
	ord_t operator/(const ord_t &o) { return ord_t(this->get_max_order()); }

	ord_t operator/(double d) { return *this; }

	ord_t operator+=(const ord_t &o) { this->order = std::max(this->order, o.order); return *this; }

protected:
	int order;
};

inline ord_t operator/(const scalar &a, const ord_t &b) { return ord_t(b.get_max_order()); }
inline ord_t operator*(const scalar &a, const ord_t &b) { return b; }
inline ord_t operator+(const scalar &a, const ord_t &b) { return b; }
inline ord_t operator-(const scalar &a, const ord_t &b) { return b; }
inline ord_t operator-(const ord_t &a) { return a; }

inline ord_t pow(const ord_t &a, const double &b) { return ord_t((int) ceil(fabs(b)) * a.get_order()); }
inline ord_t sqrt(const ord_t &a) { return a; }
inline ord_t sqr(const ord_t &a) { return ord_t(2 * a.get_order()); }
inline ord_t conj(const ord_t &a) { return a; }

#ifdef H3D_COMPLEX
namespace std {
	inline ord_t conj(const ord_t &a) { return a; }
	inline ord_t abs(const ord_t &a) { return a; }
};
#endif

inline ord_t atan2(const ord_t &a, const ord_t &b) { return ord_t(a.get_max_order()); }
inline ord_t atan(const ord_t &a) { return ord_t(a.get_max_order()); }
inline ord_t sin(const ord_t &a) { return ord_t(a.get_max_order()); }
inline ord_t cos(const ord_t &a) { return ord_t(a.get_max_order()); }


// Function
template<typename T>
class fn_t {
public:
	int nc;							// number of components
	T *fn;							// function values
	T *dx, *dy, *dz;				// derivatives

	T *fn0, *fn1, *fn2;				// components of function values
	T *dx0, *dx1, *dx2;				// components of derivatives
	T *dy0, *dy1, *dy2;
	T *dz0, *dz1, *dz2;

	T *curl0, *curl1, *curl2;		// components of curl

	fn_t() {
		fn = fn0 = fn1 = fn2 = NULL;
		dx = dx0 = dx1 = dx2 = NULL;
		dy = dy0 = dy1 = dy2 = NULL;
		dz = dz0 = dz1 = dz2 = NULL;
		curl0 = curl1 = curl2 = NULL;
	}
};


typedef fn_t<double> sfn_t;				// used for transformed shape functions
typedef fn_t<scalar> mfn_t;				// used for transformed mesh functions


// Geometry of the element
template<typename T>
class geom_t {
public:
	int marker;					// element/boundary marker
	T *x, *y, *z;				// coordinates
	T *nx, *ny, *nz;			// normals
	T *tx, *ty, *tz;			// tangents

	geom_t() {
		x = y = z = NULL;
		nx = ny = nz = NULL;
		tx = ty = tz = NULL;
	}
};

/// Init element geometry for calculating the integration order
geom_t<ord_t> init_geom(int marker);
/// Init element geometry for volumetric integrals
geom_t<double> init_geom(int marker, RefMap *rm, const int np, const QuadPt3D *pt);
/// Init element geometry for surface integrals
geom_t<double> init_geom(int marker, RefMap *rm, int iface, const int np, const QuadPt3D *pt);
/// Free data related to the element geometry
void free_geom(geom_t<double> *e);

/// Init the function for calculation the integration order
fn_t<ord_t> init_fn(const order3_t &order);
/// Init the function for the evaluation of the volumetric integral
sfn_t *init_fn(ShapeFunction *fu, RefMap *rm, const int np, const QuadPt3D *pt);
/// Init the function for the evaluation of the surface integral
sfn_t *init_fn(ShapeFunction *shfn, RefMap *rm, int iface, const int np, const QuadPt3D *pt);
/// Init the mesh-function for the evaluation of the volumetric/surface integral
mfn_t *init_fn(MeshFunction *f, RefMap *rm, const int np, const QuadPt3D *pt);

void free_fn(fn_t<ord_t> *f);
void free_fn(sfn_t *f);
void free_fn(mfn_t *f);

//

/// User defined data that can go to the bilin and lin forms. It also holds arbitraty number of functions, that user can use.
/// Typically, these functions are solutions from the previous time levels.
template<typename T>
class user_data_t {
public:
	int nf;						// number of functions in 'fn' array
	fn_t<T> *ext;				// array of pointers to functions

	user_data_t() {
		nf = 0;
		ext = NULL;
	}

	~user_data_t() {
		delete [] ext;
	}
};

#endif /* _FORMS_H_ */
