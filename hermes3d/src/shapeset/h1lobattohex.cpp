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
// h1lobbatohex.cc
//

#include "../h3dconfig.h"
#include "lobatto.h"
#include "h1lobattohex.h"
#include <common/error.h>
#include <common/callstack.h>
#include "matrix.h"

#include "mesh.h"

#ifdef WITH_HEX

struct h1_hex_index_t {
	unsigned type:2;
	unsigned ef:4;
	unsigned ori:3;
	unsigned x:4;
	unsigned y:4;
	unsigned z:4;

	h1_hex_index_t(int idx) {
		this->type = (idx >> 19) & 0x03;
		this->ef =   (idx >> 15) & 0x0F;
		this->ori =  (idx >> 12) & 0x07;
		this->x = (idx >> 8) & 0x0F;
		this->y = (idx >> 4) & 0x0F;
		this->z = (idx >> 0) & 0x0F;
	}

	h1_hex_index_t(int type, int ef, int x, int y, int z, int ori = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->ori = ori;
		this->type = type;
		this->ef = ef;
	}

	operator int() { return (type << 19) | (ef << 15) | (ori << 12) | (x << 8) | (y << 4) | z; }

	void to_str() { printf("hex = [type=%d, ef=%d, ori=%d, x=%d, y=%d, z=%d]", type, ef, ori, x, y, z); }
};


static int lobatto_hex_vertex_indices[] = {
	h1_hex_index_t(SHFN_VERTEX, 0, 0, 0, 0), h1_hex_index_t(SHFN_VERTEX, 0, 1, 0, 0),
	h1_hex_index_t(SHFN_VERTEX, 0, 1, 1, 0), h1_hex_index_t(SHFN_VERTEX, 0, 0, 1, 0),
	h1_hex_index_t(SHFN_VERTEX, 0, 0, 0, 1), h1_hex_index_t(SHFN_VERTEX, 0, 1, 0, 1),
	h1_hex_index_t(SHFN_VERTEX, 0, 1, 1, 1), h1_hex_index_t(SHFN_VERTEX, 0, 0, 1, 1)
};

// -- helpers -- //

static void find_permutation(int *indices, int *permut, int &num_01) {
	_F_
	for (int i = 0; i < 3; i++)
		permut[i] = i;
	num_01 = 0;
	if (indices[0] < 2) num_01++;
	if (indices[1] < 2) {
		num_01++;
		if (num_01 == 1) std::swap(permut[0], permut[1]);
	}
	if (indices[2] < 2) {
		num_01++;
		if (num_01 == 1) {
			std::swap(permut[1], permut[2]);
			std::swap(permut[0], permut[1]);
		}
		if (num_01 == 2) std::swap(permut[1], permut[2]);
	}
}


static void decompose(h1_hex_index_t index, int indices[3], int ori[3], bool swapori = true) {
	_F_
	int permut[3];
	int num_01;

	// get order in each direction
	indices[0] = index.x;
	indices[1] = index.y;
	indices[2] = index.z;

	find_permutation(indices, permut, num_01);

	ori[0] = ori[1] = ori[2] = 0;
	if (num_01 == 2) {
		// edge function
		assert(index.ori == 0 || index.ori == 1);
		ori[permut[2]] = index.ori;
	}
	else if (num_01 == 1) {
		// face function
		assert(index.ori >= 0 && index.ori < 8);
		if (index.ori % 2 == 1) ori[permut[1]] = 1;
		if (index.ori % 4 >= 2) ori[permut[2]] = 1;
		if (index.ori >= 4) {
			std::swap(indices[permut[1]], indices[permut[2]]);
			if (swapori) std::swap(ori[permut[1]], ori[permut[2]]);
		}
	}
	else {
		// vertex or bubble function
		assert(index.ori == 0);
	}
}

// -- functions that calculate values of fn, dx, dy, dz on the fly -- //

static void calc_fn_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		val[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
	}
}


static void calc_dx_values(int index, int np, QuadPt3D *pt, int component, double *dx) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dx[k] = lobatto_der_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
		if (oris[0] == 1) dx[k] = -dx[k];
	}
}


static void calc_dy_values(int index, int np, QuadPt3D *pt, int component, double *dy) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dy[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_der_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
		if (oris[1] == 1) dy[k] = -dy[k];
	}
}


static void calc_dz_values(int index, int np, QuadPt3D *pt, int component, double *dz) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dz[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_der_tab_1d[indices[2]](z);
		if (oris[2] == 1) dz[k] = -dz[k];
	}
}

#endif

// -- -- //

H1ShapesetLobattoHex::H1ShapesetLobattoHex() : Shapeset(1)
{
	_F_
#ifdef WITH_HEX
	type = H1;
	mode = MODE_HEXAHEDRON;
	num_components = 1;

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = calc_fn_values;
	shape_table_deleg[DX]  = calc_dx_values;
	shape_table_deleg[DY]  = calc_dy_values;
	shape_table_deleg[DZ]  = calc_dz_values;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	// vertices
	vertex_indices = lobatto_hex_vertex_indices;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

H1ShapesetLobattoHex::~H1ShapesetLobattoHex() {
	_F_
#ifdef WITH_HEX
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++)
		for (int ori = 0; ori < NUM_EDGE_ORIS; ori++)
			for (Word_t idx = edge_indices[edge][ori].first(); idx != INVALID_IDX; idx = edge_indices[edge][ori].next(idx))
				delete [] edge_indices[edge][ori][idx];

	for (int face = 0; face < Hex::NUM_FACES; face++)
		for (int ori = 0; ori < NUM_FACE_ORIS; ori++)
			for (Word_t idx = face_indices[face][ori].first(); idx != INVALID_IDX; idx = face_indices[face][ori].next(idx))
				delete [] face_indices[face][ori][idx];

	for (Word_t idx = bubble_indices.first(); idx != INVALID_IDX; idx = bubble_indices.next(idx))
		delete [] bubble_indices[idx];
#endif
}

order3_t H1ShapesetLobattoHex::get_order(int index) const {
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		h1_hex_index_t idx(index);
		order3_t ord(lobatto_order_1d[idx.x], lobatto_order_1d[idx.y], lobatto_order_1d[idx.z]);
		if (idx.type == SHFN_FACE && idx.ori >= 4) ord = turn_hex_face_order(idx.ef, ord);		// face function is turned due to orientation
		return ord;
	}
	else
		return get_ced_order(index);
#else
	return order3_t(0, 0, 0);
#endif
}

order3_t H1ShapesetLobattoHex::get_dcmp(int index) const
{
	if (index >= 0) {
		h1_hex_index_t idx(index);
		order3_t ord(idx.x, idx.y, idx.z);
		return ord;
	}
	else
		return order3_t(-1);
}

int H1ShapesetLobattoHex::get_shape_type(int index) const
{
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		h1_hex_index_t idx(index);
		return idx.type;
	}
	else
		return SHFN_CONSTRAINED;
#else
	return -1;
#endif
}

void H1ShapesetLobattoHex::compute_edge_indices(int edge, int ori, order1_t order) {
	_F_
#ifdef WITH_HEX
	assert(order > 1);
	int *indices = new int[order - 1];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 0, i, 0, 0, ori); break;
		case  1: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 1, 1, i, 0, ori); break;
		case  2: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 2, i, 1, 0, ori); break;
		case  3: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 3, 0, i, 0, ori); break;
		case  4: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 4, 0, 0, i, ori); break;
		case  5: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 5, 1, 0, i, ori); break;
		case  6: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 6, 1, 1, i, ori); break;
		case  7: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 7, 0, 1, i, ori); break;
		case  8: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 8, i, 0, 1, ori); break;
		case  9: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 9, 1, i, 1, ori); break;
		case 10: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 10, i, 1, 1, ori); break;
		case 11: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 11, 0, i, 1, ori); break;
		default: EXIT("Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#endif
}

void H1ShapesetLobattoHex::compute_face_indices(int face, int ori, order2_t order) {
	_F_
#ifdef WITH_HEX
	assert(order.x > 1);
	assert(order.y > 1);
	int horder = order.x, vorder = order.y;
	int *indices = new int[(horder - 1) * (vorder - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (face) {
		case 0:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 0, 0, i, j, ori);
			break;

		case 1:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 1, 1, i, j, ori);
			break;

		case 2:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 2, i, 0, j, ori);
			break;

		case 3:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 3, i, 1, j, ori);
			break;

		case 4:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 4, i, j, 0, ori);
			break;

		case 5:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 5, i, j, 1, ori);
			break;

		default:
			EXIT("Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order.get_idx()] = indices;
#endif
}

void H1ShapesetLobattoHex::compute_bubble_indices(order3_t order) {
	_F_
#ifdef WITH_HEX
	assert(order.x > 1);
	assert(order.y > 1);
	assert(order.z > 1);
	int *indices = new int[(order.x - 1) * (order.y - 1) * (order.z - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	for (unsigned int i = 2; i <= order.x; i++)
		for (unsigned int j = 2; j <= order.y; j++)
			for (unsigned int k = 2; k <= order.z; k++)
				indices[idx++] = h1_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 0);

	bubble_indices[order.get_idx()] = indices;
#endif
}

/// --- CED specific stuff ---

//
// constraints are calculated on egde 0
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_edge_combination(int ori, const order1_t &order, Part part) {
	_F_
#ifdef WITH_HEX
	Part rp = transform_edge_part(ori, part);

	// determine the interval of the edge
	double hi, lo;
	get_interval_part(rp.part, lo, hi);

	int n = get_num_edge_fns(order);							// total number of functions on the edge
	int *fn_idx = get_edge_indices(0, 0, order);				// indices of all functions on the edge

	double f_lo = get_value(FN, fn_idx[n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
	double f_hi = get_value(FN, fn_idx[n - 1], hi, -1.0, -1.0, 0);

	double **a = new_matrix<double>(n, n);
	MEM_CHECK(a);
	double *b = new double[n];
	MEM_CHECK(b);
	for (int i = 0; i < n; i++) {
		// chebyshev point
		double p = cos((i + 1) * M_PI / order);
		double r = (p + 1.0) * 0.5;
		double s = 1.0 - r;

		// matrix row
		for (int j = 0; j < n; j++)
			a[i][j] = get_value(FN, fn_idx[j], p, -1.0, -1.0, 0);

		// rhs
		b[i] = get_value(FN, fn_idx[n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;
	}

	// solve the system
	double d;
	int *iperm = new int[n];
	MEM_CHECK(iperm);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
#else
	return NULL;
#endif
}

//
// constraints are calculated on face 5
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_edge_face_combination(int ori, const order2_t &order, Part part, int dir, int variant) {
	_F_
#ifdef WITH_HEX
	Part rp = transform_face_part(ori, part);

	if (ori >= 4) dir = (dir == PART_ORI_VERT) ? PART_ORI_HORZ : PART_ORI_VERT; 			// turned face

	double c;					// constant that the lin. combination is multiplies with
	double **a, *b;				// matrix and rhs
	int n, m;					// total number of functions on horz and vert edge
	if (dir == PART_ORI_VERT) {
		double hi, lo;
		get_interval_part(rp.vert, lo, hi);
		int epart = face_to_edge_part(rp.horz);
		double x0;
		get_edge_part(epart, x0);

		int horder = order.x;
		int vorder = order.y;

		n = get_num_edge_fns(vorder);										// total number of functions on the edge
		int *edge_fn_idx[] = {
			get_edge_indices(0, 0, horder),										// indices of all functions on the edge
			get_edge_indices(0, 0, vorder)										// indices of all functions on the edge
		};

		double f_lo = get_value(FN, edge_fn_idx[1][n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
		double f_hi = get_value(FN, edge_fn_idx[1][n - 1], hi, -1.0, -1.0, 0);

		a = new_matrix<double>(n, n); MEM_CHECK(a);
		b = new double[n]; MEM_CHECK(b);
		for (int i = 0; i < n; i++) {
			// chebyshev point
			double p = cos((i + 1) * M_PI / vorder);
			double r = (p + 1.0) * 0.5;
			double s = 1.0 - r;

			for (int j = 0; j < n; j++)
				a[i][j] = get_value(FN, edge_fn_idx[1][j], p, -1.0, -1.0, 0);
			b[i] = get_value(FN, edge_fn_idx[1][n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;	// depends on the ref. domain
		}

		m = get_num_edge_fns(horder);											// total number of functions on the edge
		c = get_value(FN, edge_fn_idx[0][m - 1], x0, -1.0, -1.0, 0);
	}
	else {
		double hi, lo;
		get_interval_part(rp.horz, lo, hi);
		int epart = face_to_edge_part(rp.vert);
		double x0;
		get_edge_part(epart, x0);

		int horder = order.x;
		int vorder = order.y;

		n = get_num_edge_fns(horder);										// total number of functions on the edge
		int *edge_fn_idx[] = {
			get_edge_indices(0, 0, horder),										// indices of all functions on the edge
			get_edge_indices(0, 0, vorder)										// indices of all functions on the edge
		};

		double f_lo = get_value(FN, edge_fn_idx[0][n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
		double f_hi = get_value(FN, edge_fn_idx[0][n - 1], hi, -1.0, -1.0, 0);

		a = new_matrix<double>(n, n); MEM_CHECK(a);
		b = new double[n]; MEM_CHECK(b);
		for (int i = 0; i < n; i++) {
			// chebyshev point
			double p = cos((i+1) * M_PI / horder);
			double r = (p + 1.0) * 0.5;
			double s = 1.0 - r;

			for (int j = 0; j < n; j++)
				a[i][j] = get_value(FN, edge_fn_idx[0][j], p, -1.0, -1.0, 0);
			b[i] = get_value(FN, edge_fn_idx[0][n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;
		}

		m = get_num_edge_fns(vorder);											// total number of functions on the edge
		c = get_value(FN, edge_fn_idx[1][m - 1], x0, -1.0, -1.0, 0);
	}

	// solve the system
	double d;
	int *iperm = new int[n]; MEM_CHECK(iperm);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	for (int i = 0; i < n; i++)
		b[i] *= c;

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);

#else
	return NULL;
#endif
}

//
// constraints are calculated on face 5
//
//
//            edge2
//  v_hi +-----------+
//       |           |
// edge3 |           | edge1
//       |           |
//  v_lo +-----------+
//     h_lo  edge0  h_hi
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_face_combination(int ori, const order2_t &order, Part part, int variant) {
	_F_
#ifdef WITH_HEX
	int n = get_num_face_fns(order);										// total number of functions on the face
	int *fn_idx = get_face_indices(5, 0, order);							// indices of all functions on the face

	Part rp = transform_face_part(ori, part);

	double h_hi, h_lo, v_hi, v_lo;
	get_interval_part(rp.horz, h_lo, h_hi);				// determine the horizontal interval of the face
	get_interval_part(rp.vert, v_lo, v_hi);				// determine the vertical interval of the face

	int horder = order.x;
	int vorder = order.y;

	get_interval_part(rp.horz, h_lo, h_hi);
	get_interval_part(rp.vert, v_lo, v_hi);

	// fn. values at vertices of the face part
	double f_lo_lo = get_value(FN, fn_idx[n - 1], h_lo, v_lo, 1.0, 0);
	double f_lo_hi = get_value(FN, fn_idx[n - 1], h_lo, v_hi, 1.0, 0);
	double f_hi_lo = get_value(FN, fn_idx[n - 1], h_hi, v_lo, 1.0, 0);
	double f_hi_hi = get_value(FN, fn_idx[n - 1], h_hi, v_hi, 1.0, 0);

	// edge parts
	int *edge_fn_idx[4];
	edge_fn_idx[0] = get_edge_indices( 9, 0, vorder);
	edge_fn_idx[1] = get_edge_indices(10, 0, horder);
	edge_fn_idx[2] = get_edge_indices(11, 0, vorder);
	edge_fn_idx[3] = get_edge_indices( 8, 0, horder);

	double f_edge[4];
	f_edge[0] = get_fn_value(edge_fn_idx[0][vorder - 2],  1.0, v_lo, 1.0, 0);
	f_edge[1] = get_fn_value(edge_fn_idx[1][horder - 2], h_hi,  1.0, 1.0, 0);
	f_edge[2] = get_fn_value(edge_fn_idx[2][vorder - 2], -1.0, v_hi, 1.0, 0);
	f_edge[3] = get_fn_value(edge_fn_idx[3][horder - 2], h_lo, -1.0, 1.0, 0);

	Part hpart;
	Part vpart;
	hpart.part = rp.horz;
	vpart.part = rp.vert;
	int ced_edge_idx[4];
	ced_edge_idx[0] = get_constrained_edge_index( 8, 0, horder, hpart);
	ced_edge_idx[1] = get_constrained_edge_index( 9, 0, vorder, vpart);
	ced_edge_idx[2] = get_constrained_edge_index(10, 0, horder, hpart);
	ced_edge_idx[3] = get_constrained_edge_index(11, 0, vorder, vpart);

	double **a = new_matrix<double>(n, n);
	MEM_CHECK(a);
	double *b = new double[n];
	MEM_CHECK(b);

	//
	for (int row = 0; row < n; row++) {
		order2_t face_order = get_order(fn_idx[row]).get_face_order(5);
		int i = face_order.x;
		int j = face_order.y;

		double hp = cos((i - 1) * M_PI / horder);
		double hr = (hp + 1.0) * 0.5;
		double hs = 1.0 - hr;

		double vp = cos((j - 1) * M_PI / vorder);
		double vr = (vp + 1.0) * 0.5;
		double vs = 1.0 - vr;

		for (int k = 0; k < n; k++)
			a[row][k] = get_value(FN, fn_idx[k], hp, vp, 1.0, 0);

		// rhs
		b[row] = get_value(FN, fn_idx[n - 1], h_lo * hs + h_hi * hr, v_lo * vs + v_hi * vr, 1.0, 0)
			// subtract linear part
			- f_lo_lo * hs * vs
			- f_lo_hi * hs * vr
			- f_hi_lo * hr * vs
			- f_hi_hi * hr * vr
			// subtract residual of edge functions
			- f_edge[0] * get_constrained_value(FN, ced_edge_idx[0], hp, -1.0, 1.0, 0) * vs
			- f_edge[1] * get_constrained_value(FN, ced_edge_idx[1],  1.0, vp, 1.0, 0) * hr
			- f_edge[2] * get_constrained_value(FN, ced_edge_idx[2], hp,  1.0, 1.0, 0) * vr
			- f_edge[3] * get_constrained_value(FN, ced_edge_idx[3], -1.0, vp, 1.0, 0) * hs;
	}

	// solve the system
	double d;
	int *iperm = new int[n];
	MEM_CHECK(iperm);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
#else
	return NULL;
#endif
}
