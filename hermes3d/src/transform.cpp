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

//
// NOTES: make this class more robust
// - instead of static stack of transformation, make it dynamic
// - do not code the sub-elements indices bit-wise, make a stack for them
// - method set_transform(idx) is using one number as parameter, where all the
//   transformation history is stored. It should use (BYTE[] history, int length)
//   to be more robust
// - careful with sub_idx member. it is used in child classes for indexing cache
//   tables.
//

#include <common/error.h>
#include "h3dconfig.h"
#include "common.h"
#include "transform.h"
#include <common/callstack.h>

// TODO: Transformations for tetras
Trf Transformable::tetra_trf[] = {
	{ { 1, 1, 1 }, { 0, 0, 0 } }			// no transformation
};

// transformations for hexahedrons
Trf Transformable::hex_trf[] = {
	// 8 sub elements (XYZ)
	{ { 0.5, 0.5, 0.5 }, { -0.5, -0.5, -0.5 } },
	{ { 0.5, 0.5, 0.5 }, {  0.5, -0.5, -0.5 } },
	{ { 0.5, 0.5, 0.5 }, {  0.5,  0.5, -0.5 } },
	{ { 0.5, 0.5, 0.5 }, { -0.5,  0.5, -0.5 } },
	{ { 0.5, 0.5, 0.5 }, { -0.5, -0.5,  0.5 } },
	{ { 0.5, 0.5, 0.5 }, {  0.5, -0.5,  0.5 } },
	{ { 0.5, 0.5, 0.5 }, {  0.5,  0.5,  0.5 } },
	{ { 0.5, 0.5, 0.5 }, { -0.5,  0.5,  0.5 } },
	// 4 sub elements (XY)
	{ { 0.5, 0.5, 1.0 }, { -0.5, -0.5,  0.0 } },
	{ { 0.5, 0.5, 1.0 }, {  0.5, -0.5,  0.0 } },
	{ { 0.5, 0.5, 1.0 }, {  0.5,  0.5,  0.0 } },
	{ { 0.5, 0.5, 1.0 }, { -0.5,  0.5,  0.0 } },
	// 4 sub elements (XZ)
	{ { 0.5, 1.0, 0.5 }, { -0.5,  0.0, -0.5 } },
	{ { 0.5, 1.0, 0.5 }, {  0.5,  0.0, -0.5 } },
	{ { 0.5, 1.0, 0.5 }, {  0.5,  0.0,  0.5 } },
	{ { 0.5, 1.0, 0.5 }, { -0.5,  0.0,  0.5 } },
	// 4 sub elements (YZ)
	{ { 1.0, 0.5, 0.5 }, {  0.0, -0.5, -0.5 } },
	{ { 1.0, 0.5, 0.5 }, {  0.0,  0.5, -0.5 } },
	{ { 1.0, 0.5, 0.5 }, {  0.0,  0.5,  0.5 } },
	{ { 1.0, 0.5, 0.5 }, {  0.0, -0.5,  0.5 } },
	// 2 sub elements (X)
	{ { 0.5, 1.0, 1.0 }, { -0.5,  0.0,  0.0 } },
	{ { 0.5, 1.0, 1.0 }, {  0.5,  0.0,  0.0 } },
	// 2 sub elements (Y)
	{ { 1.0, 0.5, 1.0 }, {  0.0, -0.5,  0.0 } },
	{ { 1.0, 0.5, 1.0 }, {  0.0,  0.5,  0.0 } },
	// 2 sub elements (Z)
	{ { 1.0, 1.0, 0.5 }, {  0.0,  0.0, -0.5 } },
	{ { 1.0, 1.0, 0.5 }, {  0.0,  0.0,  0.5 } }
};

// TODO: transformations for prisms
static Trf prism_trf[] = {
	{ { 1, 1, 1 }, { 0, 0, 0 } }			// no transformation
};

static const int SUB_ELEMENT_BITS = 5;
static const int SUB_ELEMENT_MASK = 0x1F;

Transformable::Transformable() {
	_F_
    memset(stack, 0, sizeof(stack));
	reset_transform();
	element = NULL;
}

void Transformable::push_transform(int son) {
	_F_
    assert(element != NULL);
    if (top >= H3D_STACK_SIZE) EXIT("Too deep transform.");

	Trf *mat = stack + (++top);
	Trf *tr = NULL;
	switch (element->get_mode()) {
		case MODE_TETRAHEDRON: tr = tetra_trf + son; break;
		case MODE_HEXAHEDRON:  tr = hex_trf + son; break;
		case MODE_PRISM:       tr = prism_trf + son; break;
		default: EXIT(H3D_ERR_UNKNOWN_MODE);
	}

	// update transformation matrix
	mat->m[0] = ctm->m[0] * tr->m[0];
	mat->m[1] = ctm->m[1] * tr->m[1];
	mat->m[2] = ctm->m[2] * tr->m[2];

	// update translation vector
	mat->t[0] = ctm->m[0] * tr->t[0] + ctm->t[0];
	mat->t[1] = ctm->m[1] * tr->t[1] + ctm->t[1];
	mat->t[2] = ctm->m[2] * tr->t[2] + ctm->t[2];

	ctm = mat;
	sub_idx = (sub_idx << SUB_ELEMENT_BITS) + son + 1;
}

void Transformable::pop_transform() {
	_F_
	assert(top > 0);
	ctm = stack + (--top);
	sub_idx = (sub_idx - 1) >> SUB_ELEMENT_BITS;
}

void Transformable::reset_transform() {
	_F_
    stack[0].m[0] = stack[0].m[1] = stack[0].m[2] = 1.0;
    stack[0].t[0] = stack[0].t[1] = stack[0].t[2] = 0.0;
    ctm = stack;
    top = sub_idx = 0;
}

void Transformable::set_transform(uint64 idx) {
	_F_
	int son[25];
	int i = 0;

	while (idx > 0) {
		son[i++] = (idx - 1) & SUB_ELEMENT_MASK;
		idx = (idx - 1) >> SUB_ELEMENT_BITS;
	}
	reset_transform();
	for (int k = i - 1; k >= 0; k--)
		push_transform(son[k]);
}

void transform_points(const int np, const QuadPt3D *pt, const Trf *tr, QuadPt3D *tpt) {
	_F_
	for (int k = 0; k < np; k++) {
		tpt[k].x = tr->m[0] * pt[k].x + tr->t[0];
		tpt[k].y = tr->m[1] * pt[k].y + tr->t[1];
		tpt[k].z = tr->m[2] * pt[k].z + tr->t[2];
	}
}
