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
#include "shapeset.h"
#include "refdomain.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

/// TODO: move to common/shapeset?

/// Combine parts on the face
/// @return The combined part
/// @param[in] part - the part on the face
/// @param[in] finer_part - part that is added to the 'part'
/// It can be either:
///   0 - part stays the same),
///   1 - combine with the lower part of the 'part',
///   2 - combine with the higher part of the 'part'.
int combine_face_part(int part, int finer_part) {
	_F_
	assert(finer_part == 0 || finer_part == 1 || finer_part == 2);

	if (finer_part == 0) return part;							// stay the same
	else if (finer_part == 1) return get_lower_part(part);		// lower part
	else return get_higher_part(part);							// upper part
}

/// Get the part that is opposite to another part
/// @return the opposite part to the 'part' param
/// @param[in] part
int opposite_part(int part) {
	_F_
	int n;
	int m = part;
	for (n = 1; n <= part; n <<= 1)
		part -= n;
	return 3 * (n - 1) - m;
}

/// Get the endpoints of the interval
/// @param part[in] part of the interval (the number of stripe, see PIC)
/// @param lo[out] lower bound of the part the interval
/// @param hi[out] higher bound of the part the interval
void get_interval_part(int part, double &lo, double &hi) {
	_F_
	int n;											// number of pieces of the interval
	for (n = 1; n <= part; n <<= 1)
		part -= n;

	double n2 = 2.0 / n;							// length of the part
	lo = ((double) part * n2 - 1.0);
	hi = ((double) (part + 1) * n2 - 1.0);
}

/// Get the position on the edge
/// @param part[in] ID of the position on the edge (see PIC)
/// @param x[out] position on the edge
void get_edge_part(int part, double &x) {
	_F_
	if (part == 0)
		x = -1.0;
	else if (part == 1)
		x = 1.0;
	else {
		double lo, hi;
		get_interval_part(part - 2, lo, hi);
		x = (lo + hi) / 2.0;
	}
}

/// Transform the edge part
/// @return The transformed part
/// @param[in] ori - the orientation of the edge
/// @param[in] part - the part that is being transformed
Part transform_edge_part(int ori, Part part) {
	_F_
	Part rp;
	rp.part = (ori == 0) ? part.part : opposite_part(part.part);
	return rp;
}

/// Transform the face part
/// @return The transformed part
/// @param[in] ori - the orientation of the face
/// @param[in] part - the part that is being transformed
Part transform_face_part(int ori, Part part) {
	_F_
	// refer to Pavel Solin's gray book, p. 169 (?)
	int flags[8][3] = {
		{ 1, 1, 1 }, { -1, 1, 1 }, { 1, -1, 1 }, { -1, -1, 1 }, { 1, 1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, -1 }
	};

	Part rp;
	if (flags[ori][2] == 1) {
		rp.horz = (flags[ori][0] > 0) ? part.horz : opposite_part(part.horz);
		rp.vert = (flags[ori][1] > 0) ? part.vert : opposite_part(part.vert);
	}
	else {
		// switch hpart and vpart
		rp.horz = (flags[ori][1] > 0) ? part.vert : opposite_part(part.vert);
		rp.vert = (flags[ori][0] > 0) ? part.horz : opposite_part(part.horz);
	}

	return rp;
}


// Shapeset /////

Shapeset::Shapeset(int id) : id(id)
{
	_F_
	mode = 0;
	ced_idx = -1;
	num_components = -1;

#ifdef PRELOADING
	fn_prods = NULL;
	dx_prods = NULL;
	dy_prods = NULL;
	dz_prods = NULL;
#endif
}

Shapeset::~Shapeset() {
	_F_
#ifdef PRELOADING
	delete [] fn_prods;
	delete [] dx_prods;
	delete [] dy_prods;
	delete [] dz_prods;
#endif
	free_constrained_combinations();
}

int Shapeset::get_constrained_edge_index(int edge, int ori, order1_t order, Part part) {
	_F_
	CEDKey cedkey(CED_KEY_TYPE_EDGE, edge, order, ori, part);
	int fn_idx;
	if (ced_id.lookup(cedkey, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = cedkey;
		ced_id.set(cedkey, ced_idx);
		return -1 - ced_idx;
	}
}

int Shapeset::get_constrained_edge_face_index(int edge, int ori, order2_t order, Part part, int dir, int variant) {
	_F_
	CEDKey ck(CED_KEY_TYPE_EDGE_FACE, edge, order, ori, part, dir, variant);
	int fn_idx;
	if (ced_id.lookup(ck, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = ck;
		ced_id.set(ck, ced_idx);
		return -1 - ced_idx;
	}
}

int Shapeset::get_constrained_face_index(int face, int ori, order2_t order, Part part, int variant) {
	_F_
	CEDKey cedkey(CED_KEY_TYPE_FACE, face, order, ori, part, 0, variant);
	int fn_idx;
	if (ced_id.lookup(cedkey, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = cedkey;
		ced_id.set(cedkey, ced_idx);
		return -1 - ced_idx;
	}
}

void Shapeset::free_constrained_combinations() {
	_F_
	for (Word_t i = ced_comb.first(); i != INVALID_IDX; i = ced_comb.next(i))
		delete ced_comb.get(i);
	ced_id.remove_all();
	ced_key.remove_all();
	ced_idx = -1;
}

CEDComb *Shapeset::get_ced_comb(const CEDKey &key) {
	_F_
	CEDComb *comb;
	if (ced_comb.lookup(key, comb)) {
		// ok, already calculated combination
	}
	else {
		// combination does not exist yet => calculate it
		if (key.type == CED_KEY_TYPE_EDGE)           comb = calc_constrained_edge_combination(key.ori, key.order, key.part);
		else if (key.type == CED_KEY_TYPE_EDGE_FACE) comb = calc_constrained_edge_face_combination(key.ori, order2_t::from_int(key.order), key.part, key.dir, key.variant);
		else if (key.type == CED_KEY_TYPE_FACE)      comb = calc_constrained_face_combination(key.ori, order2_t::from_int(key.order), key.part, key.variant);
		else EXIT("Unknown type of CED key.");

		ced_comb.set(key, comb);
	}

	return comb;
}

int *Shapeset::get_ced_indices(const CEDKey &key) {
	_F_
	int *idx;
	if (key.type == CED_KEY_TYPE_EDGE) {
		order1_t order = key.order;
		idx = get_edge_indices(key.edge, key.ori, order);
	}
	else if (key.type == CED_KEY_TYPE_EDGE_FACE) {
		// HEX specific
		int dir = key.dir;
		const int *eori = RefHex::get_face_edge_orientation(key.ori);
		if (key.ori >= 4) dir = (key.dir == PART_ORI_VERT) ? PART_ORI_HORZ : PART_ORI_VERT; 			// turned face
		order2_t o = order2_t::from_int(key.order);
		idx = (dir == PART_ORI_HORZ) ? get_edge_indices(key.edge, eori[0], o.x) : get_edge_indices(key.edge, eori[1], o.y);
	}
	else if (key.type == CED_KEY_TYPE_FACE) {
		order2_t order = order2_t::from_int(key.order);
		idx = get_face_indices(key.face, key.ori, order);
	}
	else
		EXIT("Unknown type of CED key.");

	return idx;
}

void Shapeset::get_constrained_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals) {
	_F_
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	CEDComb *comb = get_ced_comb(key);
	assert(comb != NULL);
	int *idx = get_ced_indices(key);
	assert(idx != NULL);

	memset(vals, 0, np * sizeof(double));
	double tmp[np];
	for (int i = 0; i < comb->n; i++) {
		get_values(n, idx[i], np, pt, component, tmp);
		for (int j = 0; j < np; j++)
			vals[j] += comb->coef[i] * tmp[j];
	}
}

double Shapeset::get_constrained_value(int n, int index, double x, double y, double z, int component) {
	_F_
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	CEDComb *comb = get_ced_comb(key);
	assert(comb != NULL);
	int *idx = get_ced_indices(key);
	assert(idx != NULL);

	double val = 0.0;
	for (int i = 0; i < comb->n; i++)
		val += comb->coef[i] * get_value(n, idx[i], x, y, z, component);

	return val;
}

order3_t Shapeset::get_ced_order(int index) const {
	_F_
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	order3_t order;
	if (key.type == CED_KEY_TYPE_EDGE || key.type == CED_KEY_TYPE_EDGE_FACE) {
		int o;
		if (key.type == CED_KEY_TYPE_EDGE_FACE) {
			int dir;
			if (key.ori < 4) dir = key.dir;
			else dir = (key.dir == PART_ORI_VERT) ? PART_ORI_HORZ : PART_ORI_VERT; 			// turned face
			order2_t fo = order2_t::from_int(key.order);
			o = (dir == PART_ORI_HORZ) ? fo.x : fo.y;
		}
		else
			o = key.order;

		switch (key.edge) {
			case  0: order = order3_t(o, 1, 1); break;
			case  1: order = order3_t(1, o, 1); break;
			case  2: order = order3_t(o, 1, 1); break;
			case  3: order = order3_t(1, o, 1); break;
			case  4: order = order3_t(1, 1, o); break;
			case  5: order = order3_t(1, 1, o); break;
			case  6: order = order3_t(1, 1, o); break;
			case  7: order = order3_t(1, 1, o); break;
			case  8: order = order3_t(o, 1, 1); break;
			case  9: order = order3_t(1, o, 1); break;
			case 10: order = order3_t(o, 1, 1); break;
			case 11: order = order3_t(1, o, 1); break;
			default: EXIT("Invalid edge number %d. Can be 0 - 11.", key.edge); break;
		}
	}
	else if (key.type == CED_KEY_TYPE_FACE) {
		order2_t o = order2_t::from_int(key.order);

		switch (key.face) {
			case 0: order = order3_t(1, o.x, o.y); break;
			case 1: order = order3_t(1, o.x, o.y); break;
			case 2: order = order3_t(o.x, 1, o.y); break;
			case 3: order = order3_t(o.x, 1, o.y); break;
			case 4: order = order3_t(o.x, o.y, 1); break;
			case 5: order = order3_t(o.x, o.y, 1); break;
			default: EXIT("Invalid face number %d. Can be 0 - 5.", key.face); break;
		}
		if (key.ori >= 4) order = turn_hex_face_order(key.face, order);		// face function is turned due to orientation
	}

	return order;
}

#ifdef PRELOADING

bool Shapeset::load_prods(const char *file_name, double *&mat) {
	_F_
	FILE *file = fopen(file_name, "r");
	if (file != NULL) {
		fread(&num_fns, sizeof(num_fns), 1, file);
		int *idx_vec = new int [num_fns];
		MEM_CHECK(idx_vec);
		fread(idx_vec, sizeof(int), num_fns, file);

		if (fnidx2idx.count() == 0) {
			for (int i = 0; i < num_fns; i++)
				fnidx2idx[idx_vec[i]] = i;
		}

		mat = new double [num_fns * num_fns];
		MEM_CHECK(mat);
		fread(mat, sizeof(double), num_fns * num_fns, file);

		fclose(file);

		delete [] idx_vec;

		return true;
	}
	else
		return false;
}

bool Shapeset::preload_products() {
	_F_
	// we do not need this for fichera
//	printf("Loading FN products.\n");
//	if (!load_prods("fn-fn", fn_prods)) return false;
	printf("Loading DX products.\n");
	if (!load_prods("dx-dx", dx_prods)) return false;
	printf("Loading DY products.\n");
	if (!load_prods("dy-dy", dy_prods)) return false;
	printf("Loading DZ products.\n");
	if (!load_prods("dz-dz", dz_prods)) return false;
	return true;
}

scalar Shapeset::get_product_val(int idx1, int idx2, double *vals) {
	_F_
	assert(fnidx2idx.count() > 0 && vals != NULL);

	int idx[] = { idx1, idx2 };
	double ci[] = { 1.0, 1.0 };

	int nfn[2];
	int *fnidx[2];
	double *coef[2];
	for (int k = 0; k < 2; k++) {
		if (idx[k] < 0) {
			// ced
			assert(ced_key.exists(-1 - idx[k]));
			CEDKey key = ced_key[-1 - idx[k]];

			CEDComb *comb = get_ced_comb(key);

			fnidx[k] = get_ced_indices(key);
			nfn[k] = comb->n;
			coef[k] = comb->coef;
		}
		else {
			fnidx[k] = idx + k;
			nfn[k] = 1;
			coef[k] = ci + k;
		}
	}

	scalar res = 0.0;
	for (int i = 0; i < nfn[0]; i++) {
		scalar rj = 0.0;
		for (int j = 0; j < nfn[1]; j++)
			rj += coef[1][j] * vals[num_fns * fnidx2idx[fnidx[0][i]] + fnidx2idx[fnidx[1][j]]];
		res += coef[0][i] * rj;
	}

	return res;
}

#endif
