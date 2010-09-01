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

#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "mesh.h"
#include "quad.h"


// 3D transform
struct Trf {
	double3 m;			/// diagonal of a transformation matrix (3x3)
	double3 t;			/// translation vector (3 components - x, y, z)
};


/// Transformable is a base class for all classes that perform some kind of precalculation of
/// function values on elements. These classes (ShapeFunction, Solution, RefMap) inherit
/// from Transformable the ability to transform integration points to the sub-elements
/// of an element.
///
/// The purpose of this class to transform elements to sub-elements (multi-mesh)
///
class Transformable {
public:
	Transformable();
	virtual ~Transformable() {}

	/// Called by the assembling procedure and by other functions. In PrecalcShapeset it
	/// sets an internal variable that can be later retrieved by get_active_element().
	/// In Solution it selects the element to retrieve solution values for, etc.
	/// @param[in] e - Element associated with the function being represented by the class.
	virtual void set_active_element(Element *e) { element = e; }

	/// @return The element associated with the function being represented by the class.
	Element *get_active_element() const { return element; }

	/// Multiplies the current transformation matrix on the right by a transformation to the
	/// specified son element and pushes it on top of the matrix stack. All integration
	/// points will then be transformed to this sub-element. This process can be repeated.
	/// @param[in] son - Son element number in the range [0-25] for hexes.
	virtual void push_transform(int son);

	/// Removes the current transformation matrix from the top of the stack. The new top becomes
	/// the current transformation matrix. This returns the transform to the state before the
	/// last push_transform() was performed.
	virtual void pop_transform();

	/// Sets the current transform at once as if it was created by multiple calls to push_transform().
	/// @param[in] idx The number of the sub-element, as returned by get_transform().
	void set_transform(uint64 idx);

	/// @return The current transform index.
	uint64 get_transform() const { return sub_idx; }

	/// Empties the stack, loads identity transform.
	void reset_transform();

	/// @return The jacobian of the current transformation matrix.
	double get_transform_jacobian() const { return ctm->m[0] * ctm->m[1] * ctm->m[2]; }

	/// @return The current transformation matrix.
	Trf *get_ctm() const { return ctm; }

protected:
	static const int H3D_STACK_SIZE = 10;

	Element *element;				/// the active element

	Trf *ctm;						/// current sub-element transformation matrix
	uint64 sub_idx;					/// sub-element transformation index
	static const unsigned max_idx = 0x4000;

	Trf stack[H3D_STACK_SIZE];					/// transformation matrix stack
	int top;       					/// stack top

	static Trf hex_trf[];
	static Trf tetra_trf[];

	friend class Projection;
};

/// Transform points using a transformation
///
/// @param[in] np - the number of points to transform
/// @param[in] pt - the array points to transform
/// @param[in] trf - transformation to use
/// @param[out] tpt - the array of points to store the result
void transform_points(const int np, const QuadPt3D *pt, const Trf *trf, QuadPt3D *tpt);

#endif
