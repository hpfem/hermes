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

#ifndef _SHAPEFN_H_
#define _SHAPEFN_H_

#include "function.h"
#include "shapeset.h"

// Represents a shape function on a ref. domain
//
//
class ShapeFunction : public RealFunction {
public:
	/// Constructs a standard precalculated shapeset class.
	/// @param shapeset [in] Pointer to the shapeset to be precalculated.
	ShapeFunction(Shapeset *shapeset);

	ShapeFunction();

	/// Destructor.
	virtual ~ShapeFunction();

	void free();

	ESpaceType get_type() { assert(shapeset != NULL); return shapeset->get_type(); }

	/// Ensures subsequent calls to get_active_element() will be returning 'e'.
	/// Switches the class to the appropriate mode (triangle, quad).
	virtual void set_active_element(Element *e);

	/// Activates a shape function given by its index. The values of the shape function
	/// can then be obtained by setting the required integration rule order by calling
	/// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
	/// @param index [in] Shape index.
	void set_active_shape(int index);

	/// @return Index of the active shape (can be negative if the shape is constrained).
	int get_active_shape() const { return index; };

	/// @return Pointer to the shapeset which is being precalculated.
	Shapeset *get_shapeset() const { return shapeset; }

	void set_shapeset(Shapeset *ss);

	///
	void set_transform(ShapeFunction *shfn);

	virtual void precalculate(const int np, const QuadPt3D *pt, int mask);

protected:
	Shapeset *shapeset;
	int index;					/// index of active shape function

	/// Forces a transform without using push_transform() etc.
	/// Used by the Solution class. <b>For internal use only</b>.
	void force_transform(uint64 sub_idx, Trf *ctm) {
		this->sub_idx = sub_idx;
		this->ctm = ctm;
	}

	friend class Solution;
	friend class RefMap;
};


#endif
