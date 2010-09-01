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

#ifndef _FILTER_H_
#define _FILTER_H_

#include "solution.h"

struct UniData;

/// @defgroup filters Filters
///
/// Filtering
///

/// Abstract class for creating filters
///
/// Filter is a general postprocessing class, intended for visualization.
/// The output of Filter is an arbitrary combination of up to three input functions,
/// which usually are Solutions to PDEs, but also ExactSolutions and even other Filters.
///
/// (This class cannot be instantiated.)
///
/// @ingroup filters
class Filter : public MeshFunction {
public:
	Filter(MeshFunction *sln1);
	Filter(MeshFunction *sln1, MeshFunction *sln2);
	Filter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3);
	Filter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3, MeshFunction *sln4);
	virtual ~Filter();

	virtual void set_active_element(Element *e);
	virtual void free();

	virtual void push_transform(int son);
	virtual void pop_transform();

	virtual order3_t get_order();

protected:
	int num;
	MeshFunction* sln[4];
	uint64 sln_sub[4];
	void *tables[4];

	bool unimesh;
	UniData **unidata;

	void init();
	void copy_base(Filter* flt);
};


/// Base class for predefined simple filters
///
/// SimpleFilter is a base class for predefined simple filters (MagFilter, DiffFilter...).
/// The 'simplicity' lies in the fact that only one value per input function can be
/// combined (e.g., not a value and a derivative). If this is not sufficient, a full-fledged
/// filter must be derived from the Filter class (see VonMisesFilter). SimpleFilter is also
/// intended for the user to be able to easily create custom filters only by supplying the
/// combining function.
///
/// The user specifies the combining function, the arguments ('sln1', 'sln2', 'sln3'), and
/// optionally the 'item' for each argument, which can be any of FN_VAL_0, FN_DX_0, FN_DY_0
/// etc.
///
/// SimpleFilter is vector-valued, if at least one input function is vector-valued and
/// both components are specified in 'item', e.g., item1 = FN_DX (which is FN_DX_0 | FN_DX_1).
/// Otherwise it is scalar-valued.
///
/// @ingroup filters
class SimpleFilter : public Filter {
public:
	SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *result), MeshFunction *sln1, int item1 = FN_VAL);
    SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *val2, scalar *result), MeshFunction *sln1, MeshFunction *sln2, int item1 = FN_VAL, int item2 = FN_VAL);
    SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *val2, scalar *val3, scalar *result), MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3, int item1 = FN_VAL, int item2 = FN_VAL, int item3 = FN_VAL);

	virtual void precalculate(const int np, const QuadPt3D *pt, int mask);

protected:
	int item[3];

	void (*filter_fn_1)(int n, scalar *val1, scalar *result);
	void (*filter_fn_2)(int n, scalar *val1, scalar *val2, scalar *result);
	void (*filter_fn_3)(int n, scalar *val1, scalar *val2, scalar *val3, scalar *result);

	void init_components();
};


/// Calculates the magnitude of a vector function.
///
/// MagFilter takes three functions representing the components of a vector function and
/// calculates the vector magnitude, sqrt(x^2 + y^2 + z^2).
///
/// @ingroup filters
class MagFilter : public SimpleFilter {
public:
	MagFilter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3, int item1 = FN_VAL, int item2 = FN_VAL, int item3 = FN_VAL);
	MagFilter(MeshFunction *sln1, int item1 = FN_VAL); // for vector-valued sln1
};


/// Calculates the difference of two functions.
///
/// @ingroup filters
class DiffFilter : public SimpleFilter {
public:
	DiffFilter(MeshFunction *sln1, MeshFunction *sln2, int item1 = FN_VAL, int item2 = FN_VAL);
};


/// Calculates the sum of two functions.
///
/// @ingroup filters
class SumFilter : public SimpleFilter {
public:
	SumFilter(MeshFunction* sln1, MeshFunction* sln2, int item1 = FN_VAL, int item2 = FN_VAL);
};


/// Calculates the square of a function.
///
/// @ingroup filters
class SquareFilter : public SimpleFilter {
public:
	SquareFilter(MeshFunction *sln1, int item1 = FN_VAL);
};

/// Calculates the real part of a complex function.
///
/// @ingroup filters
class RealPartFilter : public SimpleFilter {
public:
	RealPartFilter(MeshFunction *sln1, int item1 = FN_VAL);
};

/// Calculates the imaginar part of a complex function.
///
/// @ingroup filters
class ImagPartFilter : public SimpleFilter {
public:
	ImagPartFilter(MeshFunction *sln1, int item1 = FN_VAL);
};


/// Calculates the Von Mises stress.
///
/// VonMisesFilter is a postprocessing filter for visualizing elastic stresses in a body.
/// It calculates the stress tensor and applies the Von Mises equivalent stress formula
/// to obtain the resulting stress measure.
///
/// FIXME: port to 3D
///
/// @ingroup filters
class VonMisesFilter : public Filter {
public:
	// TODO: cylindrical coordinates
    VonMisesFilter(MeshFunction *sln1, MeshFunction *sln2, double lambda, double mu, int cyl = 0, int item1 = FN_VAL, int item2 = FN_VAL);

protected:
	double lambda, mu;
	int cyl, item1, item2;
};

#endif
