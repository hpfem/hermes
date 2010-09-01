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
#include "filter.h"
#include "traverse.h"
#include <common/callstack.h>

//// Filter ////////////////////////////////////////////////////////////////////////////////////////

Filter::Filter(MeshFunction *sln1) :
	MeshFunction(NULL)
{
	_F_
	num = 1;
	sln[0] = sln1;
	init();
}

Filter::Filter(MeshFunction *sln1, MeshFunction *sln2) :
	MeshFunction(NULL)
{
	_F_
	num = 2;
	sln[0] = sln1;
	sln[1] = sln2;
	init();
}

Filter::Filter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3) :
	MeshFunction(NULL)
{
	_F_
	num = 3;
	sln[0] = sln1;
	sln[1] = sln2;
	sln[2] = sln3;
	init();
}

Filter::Filter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3, MeshFunction *sln4) :
	MeshFunction(NULL)
{
	_F_
	num = 4;
	sln[0] = sln1;
	sln[1] = sln2;
	sln[2] = sln3;
	sln[3] = sln4;
	init();
}

Filter::~Filter() {
	_F_
	free();
	if (unimesh) {
		delete mesh;
		for (int i = 0; i < num; i++)
			delete [] unidata[i];
		delete [] unidata;
	}
}

void Filter::init() {
	_F_
	// construct the union mesh, if necessary
	Mesh *meshes[4] = {
		sln[0]->get_mesh(),
		(num >= 2) ? sln[1]->get_mesh() : NULL,
		(num >= 3) ? sln[2]->get_mesh() : NULL,
		(num >= 4) ? sln[3]->get_mesh() : NULL
	};

	mesh = meshes[0];
	unimesh = false;

	// FIXME: better detection of same meshes
	for (int i = 1; i < num; i++)
		if (meshes[0] != meshes[1]) unimesh = true;

	if (unimesh) {
		Traverse trav;
		trav.begin(num, meshes);
		mesh = new Mesh;
		MEM_CHECK(mesh);
		unidata = trav.construct_union_mesh(mesh);
		trav.finish();
	}

	refmap->set_mesh(mesh);

	// misc init
	num_components = 1;
	memset(tables, 0, sizeof(tables));
	memset(sln_sub, 0, sizeof(sln_sub));
}

void Filter::set_active_element(Element *e) {
	_F_
	MeshFunction::set_active_element(e);
	if (!unimesh) {
		for (int i = 0; i < num; i++)
			sln[i]->set_active_element(e);
		memset(sln_sub, 0, sizeof(sln_sub));
	}
	else {
		for (int i = 0; i < num; i++) {
			sln[i]->set_active_element(unidata[i][e->id].e);
			sln[i]->set_transform(unidata[i][e->id].idx);
			sln_sub[i] = sln[i]->get_transform();
		}
	}

	switch (mode) {
		case MODE_TETRAHEDRON: order = order3_t(H3D_MAX_QUAD_ORDER_TETRA); break;
		case MODE_HEXAHEDRON: order = order3_t(H3D_MAX_QUAD_ORDER, H3D_MAX_QUAD_ORDER, H3D_MAX_QUAD_ORDER); break;
		default: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
	}
}

void Filter::free() {
	_F_
}

void Filter::push_transform(int son) {
	_F_
	MeshFunction::push_transform(son);
	for (int i = 0; i < num; i++) {
		// sln_sub[i] contains the value sln[i]->sub_idx, which the Filter thinks
		// the solution has, or at least had last time. If the filter graph is
		// cyclic, it could happen that some solutions would get pushed the transform
		// more than once. This mechanism prevents it. If the filter sees that the
		// solution already has a different sub_idx than it thinks it should have,
		// it assumes someone else has already pushed the correct transform. This
		// is actually the case for cyclic filter graphs and filters used in multi-mesh
		// computation.

		if (sln[i]->get_transform() == sln_sub[i])
			sln[i]->push_transform(son);
		sln_sub[i] = sln[i]->get_transform();
	}
}

void Filter::pop_transform() {
	_F_
	MeshFunction::pop_transform();
	for (int i = 0; i < num; i++) {
		if (sln[i]->get_transform() == sln_sub[i])
			sln[i]->pop_transform();
		sln_sub[i] = sln[i]->get_transform();
	}
}

order3_t Filter::get_order()
{
	_F_
	switch (element->get_mode()) {
		case MODE_HEXAHEDRON: return order3_t(10, 10, 10);
		case MODE_TETRAHEDRON: return order3_t(10);
		default: EXIT(H3D_ERR_NOT_IMPLEMENTED); return order3_t(10);
	}
}

//// SimpleFilter //////////////////////////////////////////////////////////////////////////////////

SimpleFilter::SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *result), MeshFunction *sln1, int item1) :
	Filter(sln1)
{
	_F_
	item[0] = item1;
	filter_fn_1 = filter_fn;
	init_components();
}

SimpleFilter::SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *val2, scalar *result), MeshFunction *sln1, MeshFunction *sln2, int item1, int item2) :
	Filter(sln1, sln2)
{
	_F_
	item[0] = item1;
	item[1] = item2;
	filter_fn_2 = filter_fn;
	init_components();
}

SimpleFilter::SimpleFilter(void (*filter_fn)(int n, scalar *val1, scalar *val2, scalar *val3, scalar *result), MeshFunction *sln1, MeshFunction* sln2, MeshFunction* sln3, int item1, int item2, int item3) :
	Filter(sln1, sln2, sln3)
{
	_F_
	item[0] = item1;
	item[1] = item2;
	item[2] = item3;
	filter_fn_3 = filter_fn;
	init_components();
}

void SimpleFilter::init_components() {
	_F_
	bool vec1 = false, vec2 = false;
	for (int i = 0; i < num; i++) {
		if (sln[i]->get_num_components() > 1) vec1 = true;
		if ((item[i] & FN_COMPONENT_0) && (item[i] & FN_COMPONENT_1) && (item[i] & FN_COMPONENT_2)) vec2 = true;
		if (sln[i]->get_num_components() == 1) item[i] &= FN_COMPONENT_0;
	}
	num_components = (vec1 && vec2) ? 3 : 1;
}

void SimpleFilter::precalculate(const int np, const QuadPt3D *pt, int mask) {
	_F_
	if (mask & (FN_DX | FN_DY | FN_DZ | FN_DXX | FN_DYY | FN_DZZ | FN_DXY | FN_DXZ | FN_DYZ)) {
		warning("Filter not defined for derivatives.");
		return;
	}

	Node *node = new_node(FN_VAL, np);

	// precalculate all solutions
	for (int i = 0; i < num; i++)
		sln[i]->precalculate(np, pt, item[i]);

	for (int j = 0; j < num_components; j++) {
		// obtain corresponding tables
		scalar *tab[3];
		for (int i = 0; i < num; i++) {
			int a = 0, b = 0;
			mask_to_comp_val(item[i], a, b);
			tab[i] = sln[i]->get_values(num_components == 1 ? a : j, b);
			if (tab[i] == NULL) {
				warning("'item%d' is incorrect in filter definition.", i + 1);
				return;
			}
		}

		// apply the filter
		switch (num) {
			case 1: filter_fn_1(np, tab[0], node->values[j][0]); break;
			case 2: filter_fn_2(np, tab[0], tab[1], node->values[j][0]); break;
			case 3: filter_fn_3(np, tab[0], tab[1], tab[2], node->values[j][0]); break;
			default: assert(false);
		}
	}

	// remove the old node and attach the new one
	replace_cur_node(node);
}

//// predefined simple filters /////////////////////////////////////////////////////////////////////

static void magnitude_fn_3(int n, scalar *v1, scalar *v2, scalar *v3, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = sqrt(sqr(v1[i]) + sqr(v2[i]) + sqr(v3[i]));
}

MagFilter::MagFilter(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *sln3, int item1, int item2, int item3) :
	SimpleFilter(magnitude_fn_3, sln1, sln2, sln3, item1, item2, item3)
{
	_F_
}

MagFilter::MagFilter(MeshFunction *sln1, int item1) :
	SimpleFilter(magnitude_fn_3, sln1, sln1, sln1, item1 & FN_COMPONENT_0, item1 & FN_COMPONENT_1, item1 & FN_COMPONENT_2)
{
	_F_
	if (sln1->get_num_components() < 3)
		EXIT("The single-argument constructor is intended for vector-valued solutions.");

}

static void difference_fn_2(int n, scalar *v1, scalar *v2, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = v1[i] - v2[i];
}

DiffFilter::DiffFilter(MeshFunction *sln1, MeshFunction *sln2, int item1, int item2) :
	SimpleFilter(difference_fn_2, sln1, sln2, item1, item2)
{
	_F_
}

static void sum_fn_2(int n, scalar *v1, scalar *v2, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = v1[i] + v2[i];
}

SumFilter::SumFilter(MeshFunction *sln1, MeshFunction *sln2, int item1, int item2) :
	SimpleFilter(sum_fn_2, sln1, sln2, item1, item2)
{
	_F_
}

static void square_fn_1(int n, scalar *v1, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = sqr(v1[i]);
}

SquareFilter::SquareFilter(MeshFunction *sln1, int item1) :
	SimpleFilter(square_fn_1, sln1, item1)
{
	_F_
}

//// Filters for visualisation of complex solutions
///////////////////////////////////////////////////////////////////////////////
static void real_part_1(int n, scalar *v1, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = REAL(v1[i]);
}

RealPartFilter::RealPartFilter(MeshFunction *sln1, int item1) :
		SimpleFilter(real_part_1, sln1, item1)
{
	_F_
}

static void imag_part_1(int n, scalar *v1, scalar *result) {
	for (int i = 0; i < n; i++)
		result[i] = IMAG(v1[i]);
}

ImagPartFilter::ImagPartFilter(MeshFunction *sln1, int item1) :
		SimpleFilter(imag_part_1, sln1, item1)
{
	_F_
}



//// VonMisesFilter ///////////////////////////////////////////////////////////////////////////////

#ifndef H3D_COMPLEX
#define getval(exp) (exp)
#else
#define getval(exp) (exp.real())
#endif

VonMisesFilter::VonMisesFilter(MeshFunction *sln1, MeshFunction *sln2, double lambda, double mu, int cyl, int item1, int item2) :
	Filter(sln1, sln2)
{
	_F_
	this->mu = mu;
	this->lambda = lambda;
	this->cyl = cyl;
	this->item1 = item1;
	this->item2 = item2;
}
