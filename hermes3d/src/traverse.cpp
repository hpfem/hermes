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
#include "mesh.h"
#include "transform.h"
#include "traverse.h"
#include <common/error.h>
#include <common/callstack.h>

#define PRINTF(...)
//#define PRINTF printf

const uint64 ONE = (uint64) 1 << 63;

struct Box {
	uint64 x_lo, x_hi;
	uint64 y_lo, y_hi;
	uint64 z_lo, z_hi;
};


struct State {
	bool visited;
	Element **e;
	Box cr;
	Box *er;
	int *trans;
};

static const int SON_BITS = 5;

///
inline uint64 next_idx(uint64 idx, int son) {
	return ((idx << SON_BITS) + son + 1);
}


static int get_hex_split_and_sons(Element *e, Box *cr, Box *er, int *sons) {
	_F_
	PRINTF("get_hex_split_and_sons\n");

	PRINTF(" * element # = %d, reft = %d\n", e->id, e->reft);
	PRINTF(" * cr = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n", cr->x_lo, cr->x_hi, cr->y_lo, cr->y_hi, cr->z_lo, cr->z_hi);
	PRINTF(" * er = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n", er->x_lo, er->x_hi, er->y_lo, er->y_hi, er->z_lo, er->z_hi);

	uint64 xmid = (er->x_lo + er->x_hi) >> 1;
	uint64 ymid = (er->y_lo + er->y_hi) >> 1;
	uint64 zmid = (er->z_lo + er->z_hi) >> 1;

	PRINTF(" * xmid = %llx, ymid = %llx, zmid = %llx\n", xmid, ymid, zmid);

	switch (e->reft) {
		case H3D_REFT_HEX_X:
			if (cr->x_hi <= xmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 20), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 21), H3D_SPLIT_NONE;
			else
				return (sons[0] = sons[3] = sons[4] = sons[7] = 20, sons[1] = sons[2] = sons[5] = sons[6] = 21), H3D_SPLIT_HEX_X;
			break;

		case H3D_REFT_HEX_Y:
			if (cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 22), H3D_SPLIT_NONE;
			else if (cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 23), H3D_SPLIT_NONE;
			else
				return (sons[0] = sons[1] = sons[4] = sons[5] = 22, sons[2] = sons[3] = sons[6] = sons[7] = 23), H3D_SPLIT_HEX_Y;
			break;

		case H3D_REFT_HEX_Z:
			if (cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 24), H3D_SPLIT_NONE;
			else if (cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 25), H3D_SPLIT_NONE;
			else
				return (sons[0] = sons[1] = sons[2] = sons[3] = 24, sons[4] = sons[5] = sons[6] = sons[7] = 25), H3D_SPLIT_HEX_Z;
			break;

		case H3D_H3D_REFT_HEX_XY:
			if (cr->x_hi <= xmid && cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] =  8), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] =  9), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 10), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid && cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 11), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] =  8, sons[2] = sons[3] = sons[6] = sons[7] = 11), H3D_SPLIT_HEX_X;
			else if (cr->x_lo >= xmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] =  9, sons[2] = sons[3] = sons[6] = sons[7] = 10), H3D_SPLIT_HEX_X;
			else if (cr->y_hi <= ymid)
				return (sons[0] = sons[3] = sons[4] = sons[7] =  8, sons[1] = sons[2] = sons[5] = sons[6] =  9), H3D_SPLIT_HEX_Y;
			else if (cr->y_lo >= ymid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 11, sons[1] = sons[2] = sons[5] = sons[6] = 10), H3D_SPLIT_HEX_Y;
			else
				return (sons[0] = sons[4] = 8, sons[1] = sons[5] = 9, sons[2] = sons[6] = 10, sons[3] = sons[7] = 11), H3D_H3D_SPLIT_HEX_XY;
			break;

		case H3D_H3D_REFT_HEX_XZ:
			if (cr->x_hi <= xmid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 12), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 13), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 14), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 15), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 12, sons[4] = sons[5] = sons[6] = sons[7] = 15), H3D_SPLIT_HEX_X;
			else if (cr->x_lo >= xmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 13, sons[4] = sons[5] = sons[6] = sons[7] = 14), H3D_SPLIT_HEX_X;
			else if (cr->z_hi <= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 12, sons[1] = sons[2] = sons[6] = sons[5] = 13), H3D_SPLIT_HEX_Z;
			else if (cr->z_lo >= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 15, sons[1] = sons[2] = sons[6] = sons[5] = 14), H3D_SPLIT_HEX_Z;
			else
				return (sons[0] = sons[3] = 12, sons[1] = sons[2] = 13, sons[4] = sons[7] = 15, sons[5] = sons[6] = 14), H3D_H3D_SPLIT_HEX_XZ;
			break;

		case H3D_H3D_REFT_HEX_YZ:
			if (cr->y_hi <= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 16), H3D_SPLIT_NONE;
			else if (cr->y_lo >= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 17), H3D_SPLIT_NONE;
			else if (cr->y_lo >= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 18), H3D_SPLIT_NONE;
			else if (cr->y_hi <= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 19), H3D_SPLIT_NONE;
			else if (cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 16, sons[4] = sons[5] = sons[6] = sons[7] = 19), H3D_SPLIT_HEX_Y;
			else if (cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 17, sons[4] = sons[5] = sons[6] = sons[7] = 18), H3D_SPLIT_HEX_Y;
			else if (cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 16, sons[2] = sons[3] = sons[6] = sons[7] = 17), H3D_SPLIT_HEX_Z;
			else if (cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 19, sons[2] = sons[3] = sons[6] = sons[7] = 18), H3D_SPLIT_HEX_Z;
			else
				return (sons[0] = sons[1] = 16, sons[2] = sons[3] = 17, sons[4] = sons[5] = 19, sons[6] = sons[7] = 18), H3D_H3D_SPLIT_HEX_YZ;
			break;

		case H3D_H3D_H3D_REFT_HEX_XYZ:
			// 8 sub-elements
			if (cr->x_hi <= xmid && cr->y_hi <= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 0), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_hi <= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 1), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_lo >= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 2), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid && cr->y_lo <= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 3), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid && cr->y_hi <= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 4), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_hi <= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 5), H3D_SPLIT_NONE;
			else if (cr->x_lo >= xmid && cr->y_lo >= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 6), H3D_SPLIT_NONE;
			else if (cr->x_hi <= xmid && cr->y_lo <= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = sons[4] = sons[5] = sons[6] = sons[7] = 7), H3D_SPLIT_NONE;
			// 4 sub-elements
			else if (cr->x_hi <= xmid && cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 0, sons[4] = sons[5] = sons[6] = sons[7] = 4), H3D_H3D_SPLIT_HEX_XY;
			else if (cr->x_lo >= xmid && cr->y_hi <= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 1, sons[4] = sons[5] = sons[6] = sons[7] = 5), H3D_H3D_SPLIT_HEX_XY;
			else if (cr->x_lo >= xmid && cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 2, sons[4] = sons[5] = sons[6] = sons[7] = 6), H3D_H3D_SPLIT_HEX_XY;
			else if (cr->x_hi <= xmid && cr->y_lo >= ymid)
				return (sons[0] = sons[1] = sons[2] = sons[3] = 3, sons[4] = sons[5] = sons[6] = sons[7] = 7), H3D_H3D_SPLIT_HEX_XY;
			else if (cr->x_hi <= xmid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 0, sons[2] = sons[3] = sons[6] = sons[7] = 3), H3D_H3D_SPLIT_HEX_XZ;
			else if (cr->x_lo >= xmid && cr->z_hi <= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 1, sons[2] = sons[3] = sons[6] = sons[7] = 2), H3D_H3D_SPLIT_HEX_XZ;
			else if (cr->x_lo >= xmid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 5, sons[2] = sons[3] = sons[6] = sons[7] = 6), H3D_H3D_SPLIT_HEX_XZ;
			else if (cr->x_hi <= xmid && cr->z_lo >= zmid)
				return (sons[0] = sons[1] = sons[4] = sons[5] = 4, sons[2] = sons[3] = sons[6] = sons[7] = 7), H3D_H3D_SPLIT_HEX_XZ;
			else if (cr->y_hi <= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 0, sons[1] = sons[2] = sons[5] = sons[6] = 1), H3D_H3D_SPLIT_HEX_YZ;
			else if (cr->y_lo >= ymid && cr->z_hi <= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 3, sons[1] = sons[2] = sons[5] = sons[6] = 2), H3D_H3D_SPLIT_HEX_YZ;
			else if (cr->y_lo >= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 7, sons[1] = sons[2] = sons[5] = sons[6] = 6), H3D_H3D_SPLIT_HEX_YZ;
			else if (cr->y_hi <= ymid && cr->z_lo >= zmid)
				return (sons[0] = sons[3] = sons[4] = sons[7] = 4, sons[1] = sons[2] = sons[5] = sons[6] = 5), H3D_H3D_SPLIT_HEX_YZ;
			// 2 sub-element
			else if (cr->x_hi <= xmid)
				return (sons[0] = sons[1] = 0, sons[2] = sons[3] = 3, sons[4] = sons[5] = 4, sons[6] = sons[7] = 7), H3D_SPLIT_HEX_X;
			else if (cr->x_lo >= xmid)
				return (sons[0] = sons[1] = 1, sons[2] = sons[3] = 2, sons[4] = sons[5] = 5, sons[6] = sons[7] = 6), H3D_SPLIT_HEX_X;
			else if (cr->y_hi <= ymid)
				return (sons[0] = sons[3] = 0, sons[1] = sons[2] = 1, sons[4] = sons[7] = 4, sons[5] = sons[6] = 5), H3D_SPLIT_HEX_Y;
			else if (cr->y_lo >= ymid)
				return (sons[0] = sons[3] = 3, sons[1] = sons[2] = 2, sons[4] = sons[7] = 7, sons[5] = sons[6] = 6), H3D_SPLIT_HEX_Y;
			else if (cr->z_hi <= zmid)
				return (sons[0] = sons[4] = 0, sons[1] = sons[5] = 1, sons[2] = sons[6] = 2, sons[3] = sons[7] = 3), H3D_SPLIT_HEX_Z;
			else if (cr->z_lo >= zmid)
				return (sons[0] = sons[4] = 4, sons[1] = sons[5] = 5, sons[2] = sons[6] = 6, sons[3] = sons[7] = 7), H3D_SPLIT_HEX_Z;
			// 0 sub-elements
			else
				return (sons[0] = 0, sons[1] = 1, sons[2] = 2, sons[3] = 3, sons[4] = 4, sons[5] = 5, sons[6] = 6, sons[7] = 7), H3D_H3D_H3D_SPLIT_HEX_XYZ;
			break;
	}

	EXIT(H3D_ERR_NOT_IMPLEMENTED);
	return -1;
}

// Hex specific
static void hex_move_to_son(Box *rnew, Box *rold, int son) {
	_F_
	PRINTF("hex_move_to_son: son = %d\n", son);

	uint64 xmid = (rold->x_lo + rold->x_hi) >> 1;
	uint64 ymid = (rold->y_lo + rold->y_hi) >> 1;
	uint64 zmid = (rold->z_lo + rold->z_hi) >> 1;

	memcpy(rnew, rold, sizeof(Box));

	switch (son) {
		case  0: rnew->x_hi = xmid; rnew->y_hi = ymid; rnew->z_hi = zmid; break;
		case  1: rnew->x_lo = xmid; rnew->y_hi = ymid; rnew->z_hi = zmid; break;
		case  2: rnew->x_lo = xmid; rnew->y_lo = ymid; rnew->z_hi = zmid; break;
		case  3: rnew->x_hi = xmid; rnew->y_lo = ymid; rnew->z_hi = zmid; break;
		case  4: rnew->x_hi = xmid; rnew->y_hi = ymid; rnew->z_lo = zmid; break;
		case  5: rnew->x_lo = xmid; rnew->y_hi = ymid; rnew->z_lo = zmid; break;
		case  6: rnew->x_lo = xmid; rnew->y_lo = ymid; rnew->z_lo = zmid; break;
		case  7: rnew->x_hi = xmid; rnew->y_lo = ymid; rnew->z_lo = zmid; break;
		//
		case  8: rnew->x_hi = xmid; rnew->y_hi = ymid; break;
		case  9: rnew->x_lo = xmid; rnew->y_hi = ymid; break;
		case 10: rnew->x_lo = xmid; rnew->y_lo = ymid; break;
		case 11: rnew->x_hi = xmid; rnew->y_lo = ymid; break;
		//
		case 12: rnew->x_hi = xmid; rnew->z_hi = zmid; break;
		case 13: rnew->x_lo = xmid; rnew->z_hi = zmid; break;
		case 14: rnew->x_lo = xmid; rnew->z_lo = zmid; break;
		case 15: rnew->x_hi = xmid; rnew->z_lo = zmid; break;
		//
		case 16: rnew->y_hi = ymid; rnew->z_hi = zmid; break;
		case 17: rnew->y_lo = ymid; rnew->z_hi = zmid; break;
		case 18: rnew->y_lo = ymid; rnew->z_lo = zmid; break;
		case 19: rnew->y_hi = ymid; rnew->z_lo = zmid; break;
		//
		case 20: rnew->x_hi = xmid; break;
		case 21: rnew->x_lo = xmid; break;
		//
		case 22: rnew->y_hi = ymid; break;
		case 23: rnew->y_lo = ymid; break;
		//
		case 24: rnew->z_hi = zmid; break;
		case 25: rnew->z_lo = zmid; break;
	}
}

/// @param cr[in]
/// @param r[in]
static int get_son_idx(Box *cr, Box *r) {
	_F_
	assert(cr != NULL && r != NULL);

	uint64 xmid = (r->x_lo + r->x_hi) >> 1;
	uint64 ymid = (r->y_lo + r->y_hi) >> 1;
	uint64 zmid = (r->z_lo + r->z_hi) >> 1;

	// 8 sub-elements
	if (cr->x_hi <= xmid && cr->y_hi <= ymid && cr->z_hi <= zmid) return 0;
	else if (cr->x_lo >= xmid && cr->y_hi <= ymid && cr->z_hi <= zmid) return 1;
	else if (cr->x_lo >= xmid && cr->y_lo >= ymid && cr->z_hi <= zmid) return 2;
	else if (cr->x_hi <= xmid && cr->y_lo <= ymid && cr->z_hi <= zmid) return 3;
	else if (cr->x_hi <= xmid && cr->y_hi <= ymid && cr->z_lo >= zmid) return 4;
	else if (cr->x_lo >= xmid && cr->y_hi <= ymid && cr->z_lo >= zmid) return 5;
	else if (cr->x_lo >= xmid && cr->y_lo >= ymid && cr->z_lo >= zmid) return 6;
	else if (cr->x_hi <= xmid && cr->y_lo <= ymid && cr->z_lo >= zmid) return 7;
	// 4 sub-elements
	else if (cr->x_hi <= xmid && cr->y_hi <= ymid) return 8;
	else if (cr->x_lo >= xmid && cr->y_hi <= ymid) return 9;
	else if (cr->x_lo >= xmid && cr->y_lo >= ymid) return 10;
	else if (cr->x_hi <= xmid && cr->y_lo >= ymid) return 11;
	else if (cr->x_hi <= xmid && cr->z_hi <= zmid) return 12;
	else if (cr->x_lo >= xmid && cr->z_hi <= zmid) return 13;
	else if (cr->x_lo >= xmid && cr->z_lo >= zmid) return 14;
	else if (cr->x_hi <= xmid && cr->z_lo >= zmid) return 15;
	else if (cr->y_hi <= ymid && cr->z_hi <= zmid) return 16;
	else if (cr->y_lo >= ymid && cr->z_hi <= zmid) return 17;
	else if (cr->y_lo >= ymid && cr->z_lo >= zmid) return 18;
	else if (cr->y_hi <= ymid && cr->z_lo >= zmid) return 19;
	// 2 sub-element
	else if (cr->x_hi <= xmid) return 20;
	else if (cr->x_lo >= xmid) return 21;
	else if (cr->y_hi <= ymid) return 22;
	else if (cr->y_lo >= ymid) return 23;
	else if (cr->z_hi <= zmid) return 24;
	else if (cr->z_lo >= zmid) return 25;

	EXIT("Corrupted box definition?");
	return -1;
}

static int trans_to_son_idx(int trans) {
	_F_
	if (trans >= 0 && trans < 8) return trans & 7;
	else if (trans < 12) return trans & 3;
	else if (trans < 16) return trans & 3;
	else if (trans < 20) return trans & 3;
	else if (trans < 22) return trans & 1;
	else if (trans < 24) return trans & 1;
	else if (trans < 26) return trans & 1;
	else EXIT(H3D_ERR_NOT_IMPLEMENTED);
}

static void init_transforms(Transformable *fn, Box *cr, Box *er) {
	_F_
	PRINTF("init_transforms\n");

	Box r;
	memcpy(&r, er, sizeof(Box));

	PRINTF("cr = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n r = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n",
		cr->x_lo, cr->x_hi, cr->y_lo, cr->y_hi, cr->z_lo, cr->z_hi,
		r.x_lo, r.x_hi, r.y_lo, r.y_hi, r.z_lo, r.z_hi);

	while (cr->x_lo > r.x_lo || cr->x_hi < r.x_hi ||
			cr->y_lo > r.y_lo || cr->y_hi < r.y_hi ||
			cr->z_lo > r.z_lo || cr->z_hi < r.z_hi)
	{
		int son = get_son_idx(cr, &r);
		fn->push_transform(son);
		hex_move_to_son(&r, &r, son);

		PRINTF("cr = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n r = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n",
			cr->x_lo, cr->x_hi, cr->y_lo, cr->y_hi, cr->z_lo, cr->z_hi,
			r.x_lo, r.x_hi, r.y_lo, r.y_hi, r.z_lo, r.z_hi);
	}
}


State *Traverse::push_state() {
	_F_
	if (top >= size) EXIT("Stack overflow. Increase stack size.");

	if (stack[top].e == NULL) {
		stack[top].e = new Element *[num]; MEM_CHECK(stack[top].e);
		stack[top].er = new Box[num];      MEM_CHECK(stack[top].er);
		stack[top].trans = new int[num];   MEM_CHECK(stack[top].trans);
	}

	stack[top].visited = false;
	memset(stack[top].trans, 0, num * sizeof(int));
	return stack + top++;
}


void Traverse::set_boundary_info(State *s, bool *bnd, FacePos *fp) {
	_F_
	Element *e = s->e[0];
	Mesh *m = meshes[0];

	switch (e->get_mode()) {
		case MODE_HEXAHEDRON:
			PRINTF("set_boundary_info (elem # = %d)\n", e->id);

			PRINTF(" * cr = (%llx, %llx) x (%llx, %llx) x (%llx, %llx)\n", s->cr.x_lo, s->cr.x_hi, s->cr.y_lo, s->cr.y_hi, s->cr.z_lo, s->cr.z_hi);
			PRINTF(" * fid = (%d, %d, %d, %d, %d, %d)\n", m->get_facet_id(e, 0), m->get_facet_id(e, 1),
					m->get_facet_id(e, 2), m->get_facet_id(e, 3), m->get_facet_id(e, 4), m->get_facet_id(e, 5));

			bnd[0] = (s->cr.x_lo == 0)   && m->facets[m->get_facet_id(e, 0)]->type == Facet::OUTER;
			bnd[1] = (s->cr.x_hi == ONE) && m->facets[m->get_facet_id(e, 1)]->type == Facet::OUTER;
			bnd[2] = (s->cr.y_lo == 0)   && m->facets[m->get_facet_id(e, 2)]->type == Facet::OUTER;
			bnd[3] = (s->cr.y_hi == ONE) && m->facets[m->get_facet_id(e, 3)]->type == Facet::OUTER;
			bnd[4] = (s->cr.z_lo == 0)   && m->facets[m->get_facet_id(e, 4)]->type == Facet::OUTER;
			bnd[5] = (s->cr.z_hi == ONE) && m->facets[m->get_facet_id(e, 5)]->type == Facet::OUTER;
			break;

		case MODE_TETRAHEDRON:
			// FIXME: check that the triangle is the unit triable (0, ONE) (see case MODE_HEXAHEDRON above)
			bnd[0] = m->facets[m->get_facet_id(e, 0)]->type == Facet::OUTER;
			bnd[1] = m->facets[m->get_facet_id(e, 1)]->type == Facet::OUTER;
			bnd[2] = m->facets[m->get_facet_id(e, 2)]->type == Facet::OUTER;
			bnd[3] = m->facets[m->get_facet_id(e, 3)]->type == Facet::OUTER;
			break;

		case MODE_PRISM:
			EXIT(H3D_ERR_NOT_IMPLEMENTED);
			break;

		default:
			EXIT(H3D_ERR_UNKNOWN_MODE);
			break;
	}

	for (int iface = 0; iface < e->get_num_faces(); iface++) {
		if (bnd[iface]) {
			Word_t fid = m->get_facet_id(e, iface);
			Facet *facet = m->facets[fid];
			Boundary *b = m->boundaries[facet->right];

			fp[iface].marker = b->marker;
			fp[iface].face = iface;
		}
	}
}

void Traverse::hex_push_son_states(State *s) {
	_F_
	PRINTF("Traverse::hex_push_son_states\n");

	// obtain split types and son numbers for the current Box on all elements
	int split = 0;
	for (int i = 0; i < num; i++)
		if (!s->e[i]->active)
			split |= get_hex_split_and_sons(s->e[i], &s->cr, s->er + i, sons[i]);

	PRINTF(" - split = 0x%x\n", split);

	int son0, son1;
	State *ns = NULL;
	switch (split) {
		case H3D_SPLIT_HEX_X:
		case H3D_SPLIT_HEX_Y:
		case H3D_SPLIT_HEX_Z:
		case H3D_H3D_SPLIT_HEX_XY:
		case H3D_H3D_SPLIT_HEX_XZ:
		case H3D_H3D_SPLIT_HEX_YZ:
			int sidx[4];
			switch (split) {
				case H3D_SPLIT_HEX_X:  son0 = 20; son1 = 21; sidx[0] = 0; sidx[1] = 6; break;
				case H3D_SPLIT_HEX_Y:  son0 = 22; son1 = 23; sidx[0] = 0; sidx[1] = 6; break;
				case H3D_SPLIT_HEX_Z:  son0 = 24; son1 = 25; sidx[0] = 0; sidx[1] = 6; break;
				case H3D_H3D_SPLIT_HEX_XY: son0 =  8; son1 = 11; sidx[0] = 0; sidx[1] = 1; sidx[2] = 2; sidx[3] = 3; break;
				case H3D_H3D_SPLIT_HEX_XZ: son0 = 12; son1 = 15; sidx[0] = 0; sidx[1] = 1; sidx[2] = 5; sidx[3] = 4; break;
				case H3D_H3D_SPLIT_HEX_YZ: son0 = 16; son1 = 19; sidx[0] = 0; sidx[1] = 2; sidx[2] = 6; sidx[3] = 4; break;
			}

			for (int son = son0, k = 0; son <= son1; son++, k++) {
				State *ns = push_state();
				hex_move_to_son(&ns->cr, &s->cr, son);

				int j = sidx[k];
				for (int i = 0; i < num; i++) {
					if (s->e[i]->active) {
						ns->e[i] = s->e[i];
						ns->trans[i] = son + 1;
					}
					else {
						Word_t eid = s->e[i]->get_son(trans_to_son_idx(sons[i][j]));
						ns->e[i] = meshes[i]->elements[eid];
						hex_move_to_son(ns->er + i, s->er + i, sons[i][j]);
						if (ns->e[i]->active) ns->trans[i] = -1;
					}
				}
			}
			break;

		case H3D_H3D_H3D_SPLIT_HEX_XYZ:
			for (int son = 0; son <= 7; son++) {
				State *ns = push_state();
				hex_move_to_son(&ns->cr, &s->cr, son);

				for (int i = 0; i < num; i++) {
					if (s->e[i]->active) {
						ns->e[i] = s->e[i];
						ns->trans[i] = son + 1;
					}
					else {
						Word_t eid = s->e[i]->get_son(trans_to_son_idx(sons[i][son]));
						ns->e[i] = meshes[i]->elements[eid];
						hex_move_to_son(ns->er + i, s->er + i, sons[i][son]);
						if (ns->e[i]->active) ns->trans[i] = -1;
					}
				}
			}
			break;

		default:
			ns = push_state();
			memcpy(&ns->cr, &s->cr, sizeof(Box));

			for (int i = 0; i < num; i++) {
				if (s->e[i]->active) {
					ns->e[i] = s->e[i];
				}
				else {
					Word_t eid = s->e[i]->get_son(trans_to_son_idx(sons[i][0]));
					ns->e[i] = meshes[i]->elements[eid];						// CHECK ME
					hex_move_to_son(ns->er + i, s->er + i, sons[i][0]);
					if (ns->e[i]->active) ns->trans[i] = -1;
				}
			}
			break;
	}
}

Element **Traverse::get_next_state(bool *bnd, FacePos *fp) {
	_F_
	PRINTF("Traverse::get_next_state\n");

	while (1) {
		int i;

		// if the top state was visited already, we are returning through it:
		// undo all its transformations, pop it and continue with a non-visited one
		State *s;
		while (top > 0 && (s = stack + top-1)->visited) {
			if (fn != NULL) {
				for (i = 0; i < num; i++) {
					if (s->trans[i] > 0) {
						if (fn[i]->get_transform() == subs[i])
							fn[i]->pop_transform();
						subs[i] = fn[i]->get_transform();
					}
					else if (s->trans[i] < 0)
						fn[i]->reset_transform();
				}
			}
			top--;
		}

		// the stack is empty, take next base element
		if (top <= 0) {
			// no more base elements? we're finished
			if (id > meshes[0]->get_num_base_elements())
				return NULL;

			// push the state of the new base element
			s = push_state();
			static const Box unity = { 0, ONE, 0, ONE, 0, ONE };
			s->cr = unity;
			for (int i = 0; i < num; i++) {
				s->e[i] = meshes[i]->elements[id];
				if (s->e[i]->active && fn != NULL) fn[i]->set_active_element(s->e[i]);
				s->er[i] = unity;
				subs[i] = 0;
			}

			base = s->e[0];
			id++;
		}

		// entering a new state: perform the transformations for it
		s->visited = true;
		if (fn != NULL) {
			for (int i = 0; i < num; i++) {
				if (s->trans[i] > 0) {
					if (fn[i]->get_transform() == subs[i])
						fn[i]->push_transform(s->trans[i]-1);
					subs[i] = fn[i]->get_transform();
				}
				else if (s->trans[i] < 0) {
					PRINTF("eid = %d\n", s->e[i]->id);
					fn[i]->set_active_element(s->e[i]);
					init_transforms(fn[i], &s->cr, s->er + i);
					subs[i] = fn[i]->get_transform();
				}
			}
		}

		// is this the leaf state?
		bool leaf = true;
		for (int i = 0; i < num; i++)
			if (!s->e[i]->active) { leaf = false; break; }

		// if yes, set boundary flags and return the state
		if (leaf) {
			if (bnd != NULL) set_boundary_info(s, bnd, fp);
			return s->e;
		}

		switch (base->get_mode()) {
			case MODE_HEXAHEDRON: hex_push_son_states(s); break;

			case MODE_TETRAHEDRON:
			case MODE_PRISM:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;

			default:
				EXIT(H3D_ERR_UNKNOWN_MODE);
				break;
		}
	}
}


void Traverse::begin(int n, Mesh **meshes, Transformable **fn) {
	_F_
	assert(n > 0);
	num = n;

	this->meshes = meshes;
	this->fn = fn;

	top = 0;
	size = 256;
	stack = new State[size];
	MEM_CHECK(stack);
	memset(stack, 0, size * sizeof(State));

	sons = new int[num][8];
	MEM_CHECK(sons);
	subs = new uint64[num];
	MEM_CHECK(subs);
	id = 1;				// element IDs start with one

	// TODO: check that meshes are compatible
}


static void free_state(State *state) {
	_F_
	delete [] state->e;
	delete [] state->er;
	delete [] state->trans;
	memset(state, 0, sizeof(State));
}


void Traverse::finish() {
	_F_
	if (stack == NULL) return;

	for (int i = 0; i < size; i++)
		if (stack[i].e != NULL)
			free_state(stack + i);

	delete [] stack;
	stack = NULL;

	delete [] subs;
	delete [] sons;
}



//// union mesh ////////////////////////////////////////////////////////////////////////////////////

uint64 hex_init_idx(Box *cr, Box *er) {
	_F_
	Box r;
	memcpy(&r, er, sizeof(Box));

	uint64 idx = 0;
	while (cr->x_lo > r.x_lo || cr->x_hi < r.x_hi ||
			cr->y_lo > r.y_lo || cr->y_hi < r.y_hi ||
			cr->z_lo > r.z_lo || cr->z_hi < r.z_hi)
	{
		int son = get_son_idx(cr, &r);
		hex_move_to_son(&r, &r, son);
		idx = next_idx(idx, son);
	}

	return idx;
}

void Traverse::hex_union_rec(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni) {
	_F_
	// state arrays
	Element *e_new[num];
	Box er_new[num], cr_new;
	int sons[num][8];
	uint64 idx_new[num];
	memcpy(idx_new, idx, sizeof(idx_new));

	// obtain split types and son numbers for the current Box on all elements
	int split = 0;
	for (int i = 0; i < num; i++)
		if (!e[i]->active)
			split |= get_hex_split_and_sons(e[i], cr, er + i, sons[i] + 0);

	int sidx[4];
	int son0, son1, reft;
	switch (split) {
		case H3D_SPLIT_HEX_X:
		case H3D_SPLIT_HEX_Y:
		case H3D_SPLIT_HEX_Z:
		case H3D_H3D_SPLIT_HEX_XY:
		case H3D_H3D_SPLIT_HEX_XZ:
		case H3D_H3D_SPLIT_HEX_YZ:
			switch (split) {
				case H3D_SPLIT_HEX_X:  reft = H3D_REFT_HEX_X;  son0 = 20; son1 = 21; sidx[0] = 0; sidx[1] = 1; break;
				case H3D_SPLIT_HEX_Y:  reft = H3D_REFT_HEX_Y;  son0 = 22; son1 = 23; sidx[0] = 0; sidx[1] = 2; break;
				case H3D_SPLIT_HEX_Z:  reft = H3D_REFT_HEX_Z;  son0 = 24; son1 = 25; sidx[0] = 0; sidx[1] = 4; break;
				case H3D_H3D_SPLIT_HEX_XY: reft = H3D_H3D_REFT_HEX_XY; son0 =  8; son1 = 11; sidx[0] = 0; sidx[1] = 1; sidx[2] = 2; sidx[3] = 3; break;
				case H3D_H3D_SPLIT_HEX_XZ: reft = H3D_H3D_REFT_HEX_XZ; son0 = 12; son1 = 15; sidx[0] = 0; sidx[1] = 1; sidx[2] = 4; sidx[3] = 5; break;
				case H3D_H3D_SPLIT_HEX_YZ: reft = H3D_H3D_REFT_HEX_YZ; son0 = 16; son1 = 19; sidx[0] = 0; sidx[1] = 2; sidx[2] = 4; sidx[3] = 6; break;
			}

			unimesh->refine_element(uni->id, reft);

			for (int son = son0, k = 0; son <= son1; son++, k++) {
				hex_move_to_son(&cr_new, cr, son);
				int j = sidx[k];
				for (int i = 0; i < num; i++) {
					if (e[i]->active) {
						e_new[i] = e[i];
						idx_new[i] = next_idx(idx[i], son);
					}
					else {
						Word_t eid = e[i]->get_son(trans_to_son_idx(sons[i][j]));
						e_new[i] = meshes[i]->elements[eid];
						hex_move_to_son(er_new + i, er + i, sons[i][j]);
						if (e_new[i]->active) idx_new[i] = hex_init_idx(&cr_new, er_new + i);
					}
				}
				Word_t eid = uni->get_son(trans_to_son_idx(son));
				union_recurrent(&cr_new, e_new, er_new, idx_new, unimesh->elements[eid]);
			}
			break;

		// 8 elements
		case H3D_H3D_H3D_SPLIT_HEX_XYZ:
			unimesh->refine_element(uni->id, H3D_H3D_H3D_REFT_HEX_XYZ);
			for (int son = 0; son < 8; son++) {
				hex_move_to_son(&cr_new, cr, son);
				for (int i = 0; i < num; i++) {
					if (e[i]->active) {
						e_new[i] = e[i];
						idx_new[i] = next_idx(idx[i], son);
					}
					else {
						Word_t eid = e[i]->get_son(trans_to_son_idx(sons[i][son]));
						e_new[i] = meshes[i]->elements[eid];
						hex_move_to_son(er_new + i, er + i, sons[i][son]);
						if (e_new[i]->active) idx_new[i] = hex_init_idx(&cr_new, er_new + i);
					}
				}

				Word_t eid = uni->get_son(trans_to_son_idx(son));
				union_recurrent(&cr_new, e_new, er_new, idx_new, unimesh->elements[eid]);
			}
			break;

		default:
			// no split
			memcpy(&cr_new, cr, sizeof(Box));
			for (int i = 0; i < num; i++) {
				if (e[i]->active)
					e_new[i] = e[i];
				else {
					Word_t eid = e[i]->get_son(trans_to_son_idx(sons[i][0]));
					e_new[i] = meshes[i]->elements[eid];
					hex_move_to_son(er_new + i, er + i, sons[i][0]);
					if (e_new[i]->active) idx_new[i] = hex_init_idx(&cr_new, er_new + i);
				}
			}
			union_recurrent(&cr_new, e_new, er_new, idx_new, uni);
			break;
	}
}

void Traverse::union_recurrent(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni) {
	_F_
	// are we at the bottom?
	bool leaf = true;
	for (int i = 0; i < num; i++)
		if (!e[i]->active) { leaf = false; break; }

	if (leaf) {
		// store the element transformation indices
		if (udsize <= uni->id) {
			if (!udsize) udsize = 1024;
			while (udsize <= uni->id) udsize *= 2;
			for (int i = 0; i < num; i++)
				unidata[i] = (UniData *) realloc(unidata[i], udsize * sizeof(UniData));
		}
		for (int i = 0; i < num; i++) {
			unidata[i][uni->id].e = e[i];
			unidata[i][uni->id].idx = idx[i];
		}
	}
	else {
		switch (base->get_mode()) {
			case MODE_HEXAHEDRON:  hex_union_rec(cr, e, er, idx, uni); break;
			case MODE_TETRAHEDRON: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
 			case MODE_PRISM:       EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
			default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
		}
	}
}


UniData **Traverse::construct_union_mesh(Mesh *unimesh) {
	_F_
	int i;
	Element *e[num];
	Box er[num], cr;

	this->unimesh = unimesh;
	unimesh->copy_base(*meshes[0]);

	udsize = 0;
	unidata = new UniData *[num];
	MEM_CHECK(unidata);
	memset(unidata, 0, sizeof(UniData *) * num);

	uint64 idx[num];
	memset(idx, 0, sizeof(idx));

	for (id = 1; id <= meshes[0]->get_num_base_elements(); id++) {
		for (i = 0; i < num; i++) {
			e[i] = meshes[i]->elements[id];
			static const Box unity = { 0, ONE, 0, ONE, 0, ONE };
			cr = er[i] = unity;
		}
		base = e[0];
		union_recurrent(&cr, e, er, idx, unimesh->elements[id]);
	}

	return unidata;
}

