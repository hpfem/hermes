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

#include "h3dconfig.h"
#include "common.h"
#include "shapefn.h"
#include "function.cpp" // non-inline template members
#include <common/error.h>
#include <common/callstack.h>


ShapeFunction::ShapeFunction() : Function<double>() {
	_F_
	this->shapeset = NULL;
	this->num_components = 0;
}

ShapeFunction::ShapeFunction(Shapeset *shapeset) :
	Function<double>()
{
	_F_
	set_shapeset(shapeset);
}

ShapeFunction::~ShapeFunction() {
	_F_
	free();
}

void ShapeFunction::set_active_shape(int index) {
	_F_

	free_cur_node();
	this->index = index;
	this->order = shapeset->get_order(index);
}


void ShapeFunction::set_active_element(Element *e) {
	_F_
	if (e->get_mode() != shapeset->get_mode())
		EXIT("Using element with incorrect shapeset.");

	free_cur_node();
	element = e;
}

void ShapeFunction::free() {
	_F_
	free_cur_node();
}

void ShapeFunction::set_shapeset(Shapeset *ss) {
	_F_
	free_cur_node();
	this->shapeset = ss;
	this->num_components = ss->get_num_components();
	assert(this->num_components == 1 || this->num_components == 3);
}


void ShapeFunction::set_transform(ShapeFunction *fn) {
	_F_
	assert(fn != NULL);
	free_cur_node();

	sub_idx = fn->sub_idx;
	top = fn->top;
	stack[top] = *(fn->ctm);
	ctm = stack + top;
}

void ShapeFunction::precalculate(const int np, const QuadPt3D *pt, int mask) {
	_F_

	int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
	int newmask = mask | oldmask;
	Node *node = new_node(newmask, np);

	// precalculate all required tables
	for (int ic = 0; ic < num_components; ic++) {
		for (int j = 0; j < VALUE_TYPES; j++) {
			if (newmask & idx2mask[j][ic]) {
				// transform quadrature points
				QuadPt3D trans_pt[np];
				for (int k = 0; k < np; k++) {
					trans_pt[k].x = ctm->m[0] * pt[k].x + ctm->t[0];
					trans_pt[k].y = ctm->m[1] * pt[k].y + ctm->t[1];
					trans_pt[k].z = ctm->m[2] * pt[k].z + ctm->t[2];
				}
				shapeset->get_values(j, index, np, trans_pt, ic, node->values[ic][j]);
			}
		}
	}

	// remove the old node and attach the new one to the Judy array
	replace_cur_node(node);
}
