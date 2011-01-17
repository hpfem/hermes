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

#include "h3d_common.h"
#include "quad.h"
#include "function.h"
#include "solution.h"
#include "refmap.h"

#define callback(a)	a<double, scalar>, a<Ord, Ord>
//#define FORM_CB(a)	a<double, scalar>, a<Ord, Ord>

// Base type for orders of functions
//
// We defined a special arithmetics with this type to be able to analyze forms
// and determine the necessary integration order.  This works for forms, but it also
// works for user-defined functions.
class HERMES_API Ord {
public:
	Ord() { order = 0; }
	Ord(scalar d) { order = 0; }
	Ord(int o) { order = o; }

	int get_order() const { return order; }
	// lowest common max order is tetrahedral
	// TODO: different max order for hexes and tets
	int get_max_order() const { return H3D_MAX_QUAD_ORDER_TETRA; }

	Ord operator+(const Ord &o) { return Ord(std::max(this->order, o.order)); }
	Ord operator-(const Ord &o) { return Ord(std::max(this->order, o.order)); }
	Ord operator-(double d) { return *this; }
	Ord operator*(const Ord &o) { return Ord(this->order + o.order); }
	Ord operator/(const Ord &o) { return Ord(this->get_max_order()); }

	Ord operator/(double d) { return *this; }

	Ord operator+=(const Ord &o) { this->order = std::max(this->order, o.order); return *this; }

protected:
	int order;
};

inline Ord operator/(const scalar &a, const Ord &b) { return Ord(b.get_max_order()); }
inline Ord operator*(const scalar &a, const Ord &b) { return b; }
inline Ord operator+(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const scalar &a, const Ord &b) { return b; }
inline Ord operator-(const Ord &a) { return a; }

inline Ord pow(const Ord &a, const double &b) { return Ord((int) ceil(fabs(b)) * a.get_order()); }
inline Ord sqrt(const Ord &a) { return a; }
inline Ord sqr(const Ord &a) { return Ord(2 * a.get_order()); }
inline Ord conj(const Ord &a) { return a; }

#ifdef H3D_COMPLEX
	inline Ord abs(const Ord &a) { return a; }
#endif

inline Ord atan2(const Ord &a, const Ord &b) { return Ord(a.get_max_order()); }
inline Ord atan(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord sin(const Ord &a) { return Ord(a.get_max_order()); }
inline Ord cos(const Ord &a) { return Ord(a.get_max_order()); }


// Function
template<typename T>
class HERMES_API Func {
public:
  int num_gip;  ///< A number of integration points used by this instance.
	int nc;							// number of components
	T *val;							// function values
	T *dx, *dy, *dz;				// derivatives

	T *val0, *val1, *val2;				// components of function values
	T *dx0, *dx1, *dx2;				// components of derivatives
	T *dy0, *dy1, *dy2;
	T *dz0, *dz1, *dz2;

	T *curl0, *curl1, *curl2;		// components of curl

	Func() {
		val = val0 = val1 = val2 = NULL;
		dx = dx0 = dx1 = dx2 = NULL;
		dy = dy0 = dy1 = dy2 = NULL;
		dz = dz0 = dz1 = dz2 = NULL;
		curl0 = curl1 = curl2 = NULL;
	};
  /// Subtract arrays stored in a given attribute from the same array in provided function.
#define H3D_SUBTRACT_IF_NOT_NULL(__ATTRIB, __OTHER_FUNC) { if (__ATTRIB != NULL) { \
  assert_msg(__OTHER_FUNC.__ATTRIB != NULL, "Unable to subtract a function expansion " #__ATTRIB " is NULL in the other function."); \
  for(int i = 0; i < num_gip; i++) __ATTRIB[i] -= __OTHER_FUNC.__ATTRIB[i]; } }

  /// Calculate this -= func for each function expations and each integration point.
  /** \param[in] func A function which is subtracted from *this. A number of integratio points and a number of component has to match. */
  void subtract(const Func<T>& func) {
    assert_msg(num_gip == func.num_gip, "Unable to subtract a function due to a different number of integration points (this: %d, other: %d)", num_gip, func.num_gip);
    assert_msg(nc == func.nc, "Unable to subtract a function due to a different number of components (this: %d, other: %d)", nc, func.nc);
    if (nc == 1)
    {
      H3D_SUBTRACT_IF_NOT_NULL(val, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz, func)
    }
    if (nc == 2) {
      H3D_SUBTRACT_IF_NOT_NULL(val0, func)
      H3D_SUBTRACT_IF_NOT_NULL(val1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz1, func)
      if(func.curl0 && this->curl0)
        H3D_SUBTRACT_IF_NOT_NULL(curl0, func)
      if(func.curl1 && this->curl1)
        H3D_SUBTRACT_IF_NOT_NULL(curl1, func)
    }
    if (nc == 3)
    {
      H3D_SUBTRACT_IF_NOT_NULL(val0, func)
      H3D_SUBTRACT_IF_NOT_NULL(val1, func)
      H3D_SUBTRACT_IF_NOT_NULL(val2, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dx2, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dy2, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz0, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz1, func)
      H3D_SUBTRACT_IF_NOT_NULL(dz2, func)
      if(func.curl0 && this->curl0)
        H3D_SUBTRACT_IF_NOT_NULL(curl0, func)
      if(func.curl1 && this->curl1)
        H3D_SUBTRACT_IF_NOT_NULL(curl1, func)
      if(func.curl2 && this->curl2)
        H3D_SUBTRACT_IF_NOT_NULL(curl2, func)
    }
  };
#undef H3D_SUBTRACT_IF_NOT_NULL
};


typedef Func<double> sFunc;				// used for transformed shape functions
typedef Func<scalar> mFunc;				// used for transformed mesh functions


// Geometry of the element
template<typename T>
class HERMES_API Geom {
public:
	int marker;					// element/boundary marker
	T *x, *y, *z;				// coordinates
	T *nx, *ny, *nz;			// normals
	T *tx, *ty, *tz;			// tangents

	Geom() {
		x = y = z = NULL;
		nx = ny = nz = NULL;
		tx = ty = tz = NULL;
	}
};

/// Init element geometry for calculating the integration order
Geom<Ord> init_geom(int marker);

/// Init element geometry for volumetric integrals
Geom<double> init_geom(int marker, RefMap *rm, const int np, const QuadPt3D *pt);

/// Init element geometry for surface integrals
Geom<double> init_geom(int marker, RefMap *rm, int iface, const int np, const QuadPt3D *pt);

/// Free data related to the element geometry
void free_geom(Geom<double> *e);

/// Init the function for calculation the integration order
Func<Ord> *init_fn_ord(const Ord3 &order);

/// Init the function for the evaluation of the volumetric integral
sFunc *init_fn(ShapeFunction *fu, RefMap *rm, const int np, const QuadPt3D *pt);

/// Init the function for the evaluation of the surface integral
sFunc *init_fn(ShapeFunction *shfn, RefMap *rm, int iface, const int np, const QuadPt3D *pt);

/// Init the mesh-function for the evaluation of the volumetric/surface integral
mFunc *init_fn(MeshFunction *f, RefMap *rm, const int np, const QuadPt3D *pt);

void free_fn(Func<Ord> *f);
void free_fn(sFunc *f);
void free_fn(mFunc *f);

//

/// User defined data that can go to the bilin and lin forms. It also holds arbitraty number of functions, that user can use.
/// Typically, these functions are solutions from the previous time levels.
template<typename T>
class HERMES_API ExtData {
public:
	int nf;					// number of functions in 'fn' array
	Func<T> **fn;				// array of pointers to functions

	ExtData() {
		nf = 0;
		fn = NULL;
	}

	~ExtData() {
		delete [] fn;
	}
};

void HERMES_API free_ext_fns_ord(ExtData<Ord> * ext);

#endif /* _FORMS_H_ */
