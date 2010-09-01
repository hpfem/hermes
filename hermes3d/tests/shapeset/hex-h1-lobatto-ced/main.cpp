// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

/*
 * Test of the constrained feature of the H1 Lobatto shapeset
 * main.cc
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define EPS						10e-10

// Everything is tested on the following geometry
//
//
//      7             6             11
//        +-----------+-----------+
//       /|          /|          /|
//      / |         / |         / |
//     /  |      5 /  |     10 /  |
//  4 +-----------+-----------+   |
//    |   |       |   |       |   |
//    |   +-------|---+-------|---+
//    |  / 3      |  / 2      |  / 9
//    | /         | /         | /
//    |/          |/          |/
//    +-----------+-----------+
//   0            1            8
//
//  ^ z
//  |
//  | / y
//  |/
//  +---> x

// vertices
Point3D vtcs[] = {
	Point3D(-2, -1, -1),
	Point3D( 0, -1, -1),
	Point3D( 0,  1, -1),
	Point3D(-2,  1, -1),
	Point3D(-2, -1,  1),
	Point3D( 0, -1,  1),
	Point3D( 0,  1,  1),
	Point3D(-2,  1,  1),
	Point3D( 2, -1, -1),
	Point3D( 2,  1, -1),
	Point3D( 2, -1,  1),
	Point3D( 2,  1,  1),
};

/*
// CAUTION: do not change the following unless you know what you are doing
// the order of hexes is crucial (function get_face_perm_ori depends on the ordering)
int hex[2][24][8] = {
	// hex 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 }
	},
	// hex 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 }
	}
};
*/
// CAUTION: do not change the following unless you know what you are doing
// the order of hexes is crucial (function get_face_perm_ori depends on the ordering)
Word_t hex[48][8] = {
	// hex 1
	{  0,  1,  2,  3,  4,  5,  6,  7 },
	{  1,  2,  3,  0,  5,  6,  7,  4 },
	{  2,  3,  0,  1,  6,  7,  4,  5 },
	{  3,  0,  1,  2,  7,  4,  5,  6 },

	{  0,  3,  2,  1,  4,  7,  6,  5 },
	{  1,  0,  3,  2,  5,  4,  7,  6 },
	{  2,  1,  0,  3,  6,  5,  4,  7 },
	{  3,  2,  1,  0,  7,  6,  5,  4 },

	{  1,  5,  6,  2,  0,  4,  7,  3 },
	{  5,  6,  2,  1,  4,  7,  3,  0 },
	{  6,  2,  1,  5,  7,  3,  0,  4 },
	{  2,  1,  5,  6,  3,  0,  4,  7 },

	{  1,  2,  6,  5,  0,  3,  7,  4 },
	{  5,  1,  2,  6,  4,  0,  3,  7 },
	{  6,  5,  1,  2,  7,  4,  0,  3 },
	{  2,  6,  5,  1,  3,  7,  4,  0 },

	{  1,  0,  4,  5,  2,  3,  7,  6 },
	{  0,  4,  5,  1,  3,  7,  6,  2 },
	{  4,  5,  1,  0,  7,  6,  2,  3 },
	{  5,  1,  0,  4,  6,  2,  3,  7 },

	{  1,  5,  4,  0,  2,  6,  7,  3 },
	{  0,  1,  5,  4,  3,  2,  6,  7 },
	{  4,  0,  1,  5,  7,  3,  2,  6 },
	{  5,  4,  0,  1,  6,  7,  3,  2 },

	{  0,  3,  7,  4,  1,  2,  6,  5 },
	{  4,  0,  3,  7,  5,  1,  2,  6 },
	{  7,  4,  0,  3,  6,  5,  1,  2 },
	{  3,  7,  4,  0,  2,  6,  5,  1 },

	{  0,  4,  7,  3,  1,  5,  6,  2 },
	{  4,  7,  3,  0,  5,  6,  2,  1 },
	{  7,  3,  0,  4,  6,  2,  1,  5 },
	{  3,  0,  4,  7,  2,  1,  5,  6 },

	{  7,  6,  5,  4,  3,  2,  1,  0 },
	{  4,  7,  6,  5,  0,  3,  2,  1 },
	{  5,  4,  7,  6,  1,  0,  3,  2 },
	{  6,  5,  4,  7,  2,  1,  0,  3 },

	{  7,  4,  5,  6,  3,  0,  1,  2 },
	{  4,  5,  6,  7,  0,  1,  2,  3 },
	{  5,  6,  7,  4,  1,  2,  3,  0 },
	{  6,  7,  4,  5,  2,  3,  0,  1 },

	{  2,  6,  7,  3,  1,  5,  4,  0 },
	{  6,  7,  3,  2,  5,  4,  0,  1 },
	{  7,  3,  2,  6,  4,  0,  1,  5 },
	{  3,  2,  6,  7,  0,  1,  5,  4 },

	{  2,  3,  7,  6,  1,  0,  4,  5 },
	{  6,  2,  3,  7,  5,  1,  0,  4 },
	{  7,  6,  2,  3,  4,  5,  1,  0 },
	{  3,  7,  6,  2,  0,  4,  5,  1 }
};

//
/*
// commented out: defined in shapeset.h

/// Get the endpoints of the interval
/// @param part[in] part of the interval (the number of stripe, see PIC)
/// @param lo[out] lower bound of the part the interval
/// @param hi[out] higher bound of the part the interval
static void get_interval_part(int part, double &lo, double &hi) {
	int n;											// number of pieces of the interval
	for (n = 1; n <= part; n <<= 1)
		part -= n;

	double n2 = 2.0 / n;							// length of the part
	lo = ((double) part * n2 - 1.0);
	hi = ((double) (part + 1) * n2 - 1.0);
}

///
/// Get the position on the edge
/// @param part[in] ID of the position on the edge (see PIC)
/// @param x[out] position on the edge
static void get_edge_part(int part, double &x) {
	double lo, hi;
	get_interval_part(part, lo, hi);
	x = (lo + hi) / 2.0;
}
*/

// edge functions

bool test_edge_ced(Mesh *mesh, Shapeset *shapeset) {
//	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
	int iedge = 0; {
//		for (int ori = 0; ori < RefHex::get_edge_orientations(); ori++) {
		int ori = 0; {
			// indices of constraining functions
			int *cng_idx = shapeset->get_edge_indices(iedge, ori, H3D_MAX_ELEMENT_ORDER);

//			for (int order = 2; order <= H3D_MAX_ELEMENT_ORDER; order++) {		// edge functions on hex start with order 2
			int order = 2; {
//				for (int p = 0; p <= 14; p++) {								// 3 levels
				int p = 2; {
					Part part;
					part.part = p;
					int ced_idx = shapeset->get_constrained_edge_index(iedge, ori, order, part);

					double lo, hi;
					get_interval_part(part.part, lo, hi);

					double vfn[2] = {
						shapeset->get_fn_value(cng_idx[order - 2], lo, -1.0, -1.0, 0),
						shapeset->get_fn_value(cng_idx[order - 2], hi, -1.0, -1.0, 0)
					};

					printf("edge = %d, ori = %d, order = %02d, part = %02d (% lf - % lf), (cng_idx = %03d)\n",
						iedge, ori, order, part, lo, hi, cng_idx[order - 2]);

					int np = 10;								// number of points to check values in
					double x = lo;								// position on the part
					double h = 2.0 / np;
					for (int i = 0; i <= np; i++) {
						double ref_x = (h * i) - 1.0;			// position on the interval [-1, 1]

						// compose the edge function (linear part + ced edge function)
						double ced_val =
							vfn[0] * shapeset->get_fn_value(shapeset->get_vertex_index(0), ref_x, -1.0, -1.0, 0) +
							vfn[1] * shapeset->get_fn_value(shapeset->get_vertex_index(1), ref_x, -1.0, -1.0, 0) +
							shapeset->get_fn_value(ced_idx, ref_x, -1.0, -1.0, 0);

						double cng_val = shapeset->get_fn_value(cng_idx[order - 2], x, -1.0, -1.0, 0);

						printf(" %lf: ced = % lf | cng = % lf | (%lf) | % lf\n", x,
							ced_val, cng_val, fabs(ced_val - cng_val), shapeset->get_fn_value(ced_idx, ref_x, -1.0, -1.0, 0));

						if (fabs(ced_val - cng_val) > EPS) {
							printf("ERROR: constrained value does not fit constraning value.\n");
							return false;
						}


/*						printf(" %lf: ced = % lf * %lf + %lf * %lf + % lf = % lf | % lf | (%lf)\n", x,
							vfn[0], shapeset->get_fn_value(shapeset->get_vertex_index(0), x, -1.0, -1.0, 0),
							vfn[1], shapeset->get_fn_value(shapeset->get_vertex_index(1), x, -1.0, -1.0, 0),
							val,
							ced_val,
							cng_val, fabs(ced_val - cng_val));
*/
//						printf(" % lf: ced = % lf\n", x, ced_val);

						x += (hi - lo) / np;
					}
				}
			}
		}
	}

	return false;
}

// edge functions constrained by face functions

bool test_edge_face_ced(Mesh *mesh, Shapeset *shapeset) {
//	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
	int iface = 5; {
//		for (int ori = 0; ori < RefHex::get_face_orientations(); ori++) {
		for (int ori = 0; ori < 1; ori++) {
			printf("face %d, ori = %d\n", iface, ori);
			// indices of constraining functions
			int *cng_idx = shapeset->get_face_indices(iface, ori, MAKE_QUAD_ORDER(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER));

			int fn_idx = 0;
			for (int horder = 2; horder <= 2; horder++) {					// face functions on hex start with order 2 (in each direction)
				for (int vorder = 2; vorder <= 2; vorder++, fn_idx++) {
//			for (int horder = 2; horder <= H3D_MAX_ELEMENT_ORDER; horder++) {			// face functions on hex start with order 2 (in each direction)
//				for (int vorder = 2; vorder <= H3D_MAX_ELEMENT_ORDER; vorder++) {
					int order = MAKE_QUAD_ORDER(horder, vorder);
					printf("order = %d\n", order);

					int fpart = 2; {
						int epart = 2; {
//					for (int vpart = 0; vpart <= 0; vpart++) {
//						for (int hpart = 0; hpart <= 0; hpart++) {
							Part part;
							part.fpart = fpart;
							part.epart = epart;
							part.ori = 0;

							/// ???
							int iedge = 8;
							int ced_idx = shapeset->get_constrained_edge_face_index(iedge, ori, order, part);
//							printf("ced_idx = %d\n", ced_idx);

							double lo, hi;
							get_interval_part(part.fpart, lo, hi);
							double x0;
							get_edge_part(part.epart, x0);

							printf(" region = (% lf, % lf), x0 = % lf\n", lo, hi, x0);

							double vfn[2] = {					// fn. values at vertices
								shapeset->get_fn_value(cng_idx[fn_idx], lo, x0, 1.0, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], hi, x0, 1.0, 0),
							};

//							printf("vfn: ");
//							for (int t = 0; t < 4; t++)
//								printf("% lf ", vfn[t]);
//							printf("\n");
//
//							printf("efn: ");
//							for (int t = 0; t < 4; t++)
//								printf("% lf ", efn[t]);
//							printf("\n");

							int np = 8;
							double x = lo;								// position on the part

//							printf("\n");
							double h = 2.0 / np;
							for (int i = 0; i <= np; i++) {
								double ref_x = (h * i) - 1.0;			// position on the interval [-1, 1]

								double ced_val =
									vfn[0] * shapeset->get_fn_value(shapeset->get_vertex_index(4), ref_x, -1.0, 1.0, 0) +
									vfn[1] * shapeset->get_fn_value(shapeset->get_vertex_index(5), ref_x, -1.0, 1.0, 0) +
									shapeset->get_fn_value(ced_idx, ref_x, -1.0, 1.0, 0);

//								double val = shapeset->get_fn_value(ced_idx, -1.0, ref_x, ref_y, 0);
								printf("% lf, ", ced_val);

//								printf(" (% lf) ", shapeset->get_fn_value(ced_edge_idx[1], -1.0, ref_x, ref_y, 0));

//								printf("\n");
							}

							printf("\ncing\n");
/*
							np = 8;
							double y = -1.0;
							for (int i = 0; i <= np; i++) {
								x = lo;								// position on the part
								printf("% lf: ", y);
								for (int j = 0; j <= np; j++) {
									double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], x, y, 1.0, 0);
									printf("% lf, ", cng_val);
									x += (hi - lo) / np;
								}
								printf("\n");
								y += 2.0 / np;
							}
							printf("\n");
*/
							x = lo;
							for (int i = 0; i <= np; i++) {
								double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], x, x0, 1.0, 0);
								printf("% lf, ", cng_val);
								x += (hi - lo) / np;
							}
							printf("\n");

						}
					}
				}
			}
		}
	}

	return false;
}

bool test_face_ced(Mesh *mesh, Shapeset *shapeset) {
	int edge_ori[8][4] = {
		{ 0, 0, 0, 0 }, { 1, 0, 1, 0 }, { 0, 1, 0, 1 }, { 1, 1, 1, 1 },
		{ 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 1, 0, 1, 0 }, { 1, 1, 1, 1 }
	};

//	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
	int iface = 5; {
//		for (int ori = 0; ori < RefHex::get_face_orientations(iface); ori++) {
		int ori = 7; {
//			printf("face %d, ori = %d\n", iface, ori);
			// indices of constraining functions
//			int *cng_idx = shapeset->get_face_indices(iface, ori, MAKE_QUAD_ORDER(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER));
//			printf("------\n");

			int fn_idx = 0;
//			for (int horder = 2; horder <= 10; horder++) {					// face functions on hex start with order 2 (in each direction)
			int horder = 2; {
//				for (int vorder = 2; vorder <= 10; vorder++, fn_idx++) {
				int vorder = 3; {
					int order = MAKE_QUAD_ORDER(horder, vorder);
//					printf("order = %d\n", order);
//					printf("face %d, ori = %d, order = (%d %d)\n", iface, ori, horder, vorder);
					printf("ori = %d, order = (%d %d)\n", ori, horder, vorder);

					int *cng_idx = shapeset->get_face_indices(iface, ori, order);
					fn_idx = shapeset->get_num_face_fns(order) - 1;

					int vpart = 0; {
						int hpart = 2; {
//					for (int vpart = 0; vpart <= 510; vpart++) {				// 8 levels
//						for (int hpart = 0; hpart <= 510; hpart++) {			// 8 levels
//					for (int vpart = 0; vpart <= 16; vpart++) {					// 4 levels
//						for (int hpart = 0; hpart <= 16; hpart++) {				// 4 levels
							Part part;
							part.horz = hpart;
							part.vert = vpart;

							int ced_idx = shapeset->get_constrained_face_index(iface, ori, order, part);
//							printf("ced_idx = %d\n", ced_idx);

							double h_lo, h_hi;
							get_interval_part(part.horz, h_lo, h_hi);
							double v_lo, v_hi;
							get_interval_part(part.vert, v_lo, v_hi);

//							printf(" region = (% lf, % lf) x (% lf, % lf)\n", h_lo, h_hi, v_lo, v_hi);

							double vfn[4] = {					// fn. values at vertices
								shapeset->get_fn_value(cng_idx[fn_idx], h_lo, v_lo, 1.0, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], h_hi, v_lo, 1.0, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], h_hi, v_hi, 1.0, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], h_lo, v_hi, 1.0, 0)
							};

							int *edge_fn_idx[4];
							if (ori >= 0 && ori <= 3) {
								edge_fn_idx[0] = shapeset->get_edge_indices( 9, edge_ori[ori][1], vorder);
								edge_fn_idx[1] = shapeset->get_edge_indices(10, edge_ori[ori][2], horder);
								edge_fn_idx[2] = shapeset->get_edge_indices(11, edge_ori[ori][3], vorder);
								edge_fn_idx[3] = shapeset->get_edge_indices( 8, edge_ori[ori][0], horder);
							}
							else {
								edge_fn_idx[0] = shapeset->get_edge_indices( 9, edge_ori[ori][1], horder);
								edge_fn_idx[1] = shapeset->get_edge_indices(10, edge_ori[ori][2], vorder);
								edge_fn_idx[2] = shapeset->get_edge_indices(11, edge_ori[ori][3], horder);
								edge_fn_idx[3] = shapeset->get_edge_indices( 8, edge_ori[ori][0], vorder);
							}

							Part epart[4];
							epart[0].part = hpart;
							epart[1].part = vpart;
							epart[2].part = hpart;
							epart[3].part = vpart;

							int ced_edge_idx[4];
							double efn[4];
							if (ori >= 0 && ori <= 3) {
								ced_edge_idx[0] = shapeset->get_constrained_edge_index( 8, edge_ori[ori][0], horder, epart[0]);
								ced_edge_idx[1] = shapeset->get_constrained_edge_index( 9, edge_ori[ori][1], vorder, epart[1]);
								ced_edge_idx[2] = shapeset->get_constrained_edge_index(10, edge_ori[ori][2], horder, epart[2]);
								ced_edge_idx[3] = shapeset->get_constrained_edge_index(11, edge_ori[ori][3], vorder, epart[3]);

								efn[0] = shapeset->get_fn_value(edge_fn_idx[0][vorder - 2],  1.0, v_lo, 1.0, 0);
								efn[1] = shapeset->get_fn_value(edge_fn_idx[1][horder - 2], h_hi,  1.0, 1.0, 0);
								efn[2] = shapeset->get_fn_value(edge_fn_idx[2][vorder - 2], -1.0, v_hi, 1.0, 0);
								efn[3] = shapeset->get_fn_value(edge_fn_idx[3][horder - 2], h_lo, -1.0, 1.0, 0);
							}
							else {
								ced_edge_idx[0] = shapeset->get_constrained_edge_index( 8, edge_ori[ori][0], vorder, epart[0]);
								ced_edge_idx[1] = shapeset->get_constrained_edge_index( 9, edge_ori[ori][1], horder, epart[1]);
								ced_edge_idx[2] = shapeset->get_constrained_edge_index(10, edge_ori[ori][2], vorder, epart[2]);
								ced_edge_idx[3] = shapeset->get_constrained_edge_index(11, edge_ori[ori][3], horder, epart[3]);

								efn[0] = shapeset->get_fn_value(edge_fn_idx[0][horder - 2],  1.0, v_lo, 1.0, 0);
								efn[1] = shapeset->get_fn_value(edge_fn_idx[1][vorder - 2], h_hi,  1.0, 1.0, 0);
								efn[2] = shapeset->get_fn_value(edge_fn_idx[2][horder - 2], -1.0, v_hi, 1.0, 0);
								efn[3] = shapeset->get_fn_value(edge_fn_idx[3][vorder - 2], h_lo, -1.0, 1.0, 0);
							}
#ifdef DEBUG
//							printf("vfn: ");
//							for (int t = 0; t < 4; t++)
//								printf("% lf ", vfn[t]);
//							printf("\n");

/*							printf("efn: ");
							for (int t = 0; t < 4; t++)
								printf("% lf ", efn[t]);
							printf("\n");
*/
#endif
							int np = 5;								// number of points in one direction (matrix of points np x np)
							double x, y;								// position on the part

//							printf("\n");
							y = v_lo;
							double h = 2.0 / np;
							for (int i = 0; i <= np; i++) {
								double ref_y = (h * i) - 1.0;				// position on the interval [-1, 1]
								x = h_lo;								// position on the part
								for (int j = 0; j <= np; j++) {
									double ref_x = (h * j) - 1.0;			// position on the interval [-1, 1]

									// constrained value
									double ced_val =
										vfn[0] * shapeset->get_fn_value(shapeset->get_vertex_index(4), ref_x, ref_y, 1.0, 0) +
										vfn[1] * shapeset->get_fn_value(shapeset->get_vertex_index(5), ref_x, ref_y, 1.0, 0) +
										vfn[2] * shapeset->get_fn_value(shapeset->get_vertex_index(6), ref_x, ref_y, 1.0, 0) +
										vfn[3] * shapeset->get_fn_value(shapeset->get_vertex_index(7), ref_x, ref_y, 1.0, 0) +
										efn[0] * shapeset->get_fn_value(ced_edge_idx[0], ref_x, ref_y, 1.0, 0) +
										efn[1] * shapeset->get_fn_value(ced_edge_idx[1], ref_x, ref_y, 1.0, 0) +
										efn[2] * shapeset->get_fn_value(ced_edge_idx[2], ref_x, ref_y, 1.0, 0) +
										efn[3] * shapeset->get_fn_value(ced_edge_idx[3], ref_x, ref_y, 1.0, 0) +
										shapeset->get_fn_value(ced_idx, ref_x, ref_y, 1.0, 0);

									// constraining value
									double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], x, y, 1.0, 0);
#ifdef DEBUG
									printf("% lf, ", ced_val);
//									printf("% lf == % lf (% lf), ", cng_val, ced_val, fabs(ced_val - cng_val));
#endif
									// compare constraining and constrained values
									if (fabs(ced_val - cng_val) > EPS) {
										printf("\nERROR: constrained value does not fit constraning value (face = %d, ori = %d, order = (%d, %d), part = (%d, %d)).\n", iface, ori, horder, vorder, hpart, vpart);
										return false;
									}

//									printf("\n");

									x += (h_hi - h_lo) / np;
								}
								printf("\n");
								y += (v_hi - v_lo) / np;
							}

#ifdef DEBUG
							printf("\ncing\n");
							y = v_lo;
							for (int i = 0; i <= np; i++) {
								x = h_lo;								// position on the part
								for (int j = 0; j <= np; j++) {
									double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], x, y, 1.0, 0);
									printf("% lf, ", cng_val);
									x += (h_hi - h_lo) / np;
								}
								printf("\n");
								y += (v_hi - v_lo) / np;
							}
#endif

//							printf(".");
						}
					}
//					printf("\n");
				}
			}
		}
	}

	return false;
}

bool test_face_ced_0(Mesh *mesh, Shapeset *shapeset) {
//	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
	int iface = 0; {
//		for (int ori = 0; ori < RefHex::get_face_orientations(); ori++) {
//		for (int ori = 0; ori < 2; ori++) {
		int ori = 0; {
			printf("face %d, ori = %d\n", iface, ori);
			// indices of constraining functions
//			int *cng_idx = shapeset->get_face_indices(iface, ori, MAKE_QUAD_ORDER(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER));

//			int fn_idx = 0;
//			for (int horder = 2; horder <= 2; horder++) {					// face functions on hex start with order 2 (in each direction)
			int horder = 2; {
//				for (int vorder = 2; vorder <= 2; vorder++, fn_idx++) {
				int vorder = 2; {

					int *cng_idx = shapeset->get_face_indices(iface, ori, MAKE_QUAD_ORDER(horder, vorder));
					int fn_idx = shapeset->get_num_face_fns(MAKE_QUAD_ORDER(horder, vorder)) - 1;
//			for (int horder = 2; horder <= H3D_MAX_ELEMENT_ORDER; horder++) {			// face functions on hex start with order 2 (in each direction)
//				for (int vorder = 2; vorder <= H3D_MAX_ELEMENT_ORDER; vorder++) {
					int order = MAKE_QUAD_ORDER(horder, vorder);
					printf("order = %d\n", order);

					int vpart = 1; {
						int hpart = 1; {
//					for (int vpart = 0; vpart <= 0; vpart++) {
//						for (int hpart = 0; hpart <= 0; hpart++) {
							Part part;
							part.horz = hpart;
							part.vert = vpart;

							int ced_idx = shapeset->get_constrained_face_index(iface, ori, order, part);
//							printf("ced_idx = %d\n", ced_idx);

							double h_lo, h_hi;
							get_interval_part(part.horz, h_lo, h_hi);
							double v_lo, v_hi;
							get_interval_part(part.vert, v_lo, v_hi);

							printf(" region = (% lf, % lf) x (% lf, % lf)\n", h_lo, h_hi, v_lo, v_hi);

							double vfn[4] = {					// fn. values at vertices
								shapeset->get_fn_value(cng_idx[fn_idx], -1.0, h_lo, v_lo, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], -1.0, h_hi, v_lo, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], -1.0, h_hi, v_hi, 0),
								shapeset->get_fn_value(cng_idx[fn_idx], -1.0, h_lo, v_hi, 0)
							};

							int *edge_fn_idx[] = {
/*								shapeset->get_edge_indices(8, 0, horder),
								shapeset->get_edge_indices(8, 0, vorder),
*/								shapeset->get_edge_indices( 7, 0, vorder),
								shapeset->get_edge_indices(11, 0, horder),
								shapeset->get_edge_indices( 4, 0, vorder),
								shapeset->get_edge_indices( 3, 0, horder)
							};

							Part epart[4];
							epart[0].part = hpart;
							epart[1].part = vpart;
							epart[2].part = hpart;
							epart[3].part = vpart;

							int ced_edge_idx[4] = {
								shapeset->get_constrained_edge_index( 3, 0, horder, epart[0]),
								shapeset->get_constrained_edge_index( 7, 0, vorder, epart[1]),
								shapeset->get_constrained_edge_index(11, 0, horder, epart[2]),
								shapeset->get_constrained_edge_index( 4, 0, vorder, epart[3])
							};

							double efn[4] = {
								shapeset->get_fn_value(edge_fn_idx[0][vorder - 2], -1.0,  1.0, v_lo, 0),
								shapeset->get_fn_value(edge_fn_idx[1][horder - 2], -1.0, h_hi,  1.0, 0),
								shapeset->get_fn_value(edge_fn_idx[2][vorder - 2], -1.0, -1.0, v_hi, 0),
								shapeset->get_fn_value(edge_fn_idx[3][horder - 2], -1.0, h_lo, -1.0, 0)
							};

//							printf("vfn: ");
//							for (int t = 0; t < 4; t++)
//								printf("% lf ", vfn[t]);
//							printf("\n");
//
//							printf("efn: ");
//							for (int t = 0; t < 4; t++)
//								printf("% lf ", efn[t]);
//							printf("\n");

							printf("\n");

							int np = 5;
							double h = 2.0 / np;
//							double ref_y = v_lo;
							for (int i = 0; i <= np; i++) {
								double ref_y = (h * i) - 1.0;				// position on the interval [-1, 1]
//								double ref_x = h_lo;								// position on the part
								for (int j = 0; j <= np; j++) {
									double ref_x = (h * j) - 1.0;			// position on the interval [-1, 1]

									double ced_val =
										vfn[0] * shapeset->get_fn_value(shapeset->get_vertex_index(0), -1.0, ref_x, ref_y, 0) +
										vfn[1] * shapeset->get_fn_value(shapeset->get_vertex_index(3), -1.0, ref_x, ref_y, 0) +
										vfn[2] * shapeset->get_fn_value(shapeset->get_vertex_index(7), -1.0, ref_x, ref_y, 0) +
										vfn[3] * shapeset->get_fn_value(shapeset->get_vertex_index(4), -1.0, ref_x, ref_y, 0) +
										efn[0] * shapeset->get_fn_value(ced_edge_idx[0], -1.0, ref_x, ref_y, 0) +
										efn[1] * shapeset->get_fn_value(ced_edge_idx[1], -1.0, ref_x, ref_y, 0) +
										efn[2] * shapeset->get_fn_value(ced_edge_idx[2], -1.0, ref_x, ref_y, 0) +
										efn[3] * shapeset->get_fn_value(ced_edge_idx[3], -1.0, ref_x, ref_y, 0) +
										shapeset->get_fn_value(ced_idx, -1.0, ref_x, ref_y, 0);

//									double val = shapeset->get_fn_value(ced_idx, -1.0, ref_x, ref_y, 0);
									printf("% lf, ", ced_val);

//									printf(" (% lf) ", shapeset->get_fn_value(ced_edge_idx[1], -1.0, ref_x, ref_y, 0));
//									printf("\n");

//									x += (h_hi - h_lo) / np;
//									ref_x += (h_hi - h_lo) / np;
								}
								printf("\n");
//								y += (v_hi - v_lo) / np;
//								ref_y += (v_hi - v_lo) / np;
							}

							printf("\ncing\n");

							double y = v_lo;
							for (int i = 0; i <= np; i++) {
								double x = h_lo;								// position on the part
								for (int j = 0; j <= np; j++) {
									double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], -1.0, x, y, 0);
									printf("% lf, ", cng_val);
									x += (h_hi - h_lo) / np;
								}
								printf("\n");
								y += (v_hi - v_lo) / np;
							}
/*
							printf("---\n");
							y = -1.0;
							for (int i = 0; i <= np; i++) {
								x = -1.0;								// position on the part
								for (int j = 0; j <= 2 * np; j++) {
									double cng_val = shapeset->get_fn_value(cng_idx[fn_idx], -1.0, x, y, 0);
									printf("% lf, ", cng_val);
									x += 2.0 / (2 * np);
								}
								printf("\n");
								y += 2.0 / (np);
							}
*/
						}
					}
				}
			}
		}
	}

	return false;
}

bool test_match(Shapeset *shapeset) {
	int edge_ori[8][4] = {
		{ 0, 0, 0, 0 }, { 1, 0, 1, 0 }, { 0, 1, 0, 1 }, { 1, 1, 1, 1 },
		{ 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 1, 0, 1, 0 }, { 1, 1, 1, 1 }
	};

//		for (int ori = 0; ori < RefHex::get_face_orientations(iface); ori++) {
		int ori = 7;

//			int fn_idx = 0;
//			for (int horder = 2; horder <= 10; horder++) {					// face functions on hex start with order 2 (in each direction)
			int horder = 3; {
//				for (int vorder = 2; vorder <= 10; vorder++, fn_idx++) {
				int vorder = 4; {
					int order = MAKE_QUAD_ORDER(horder, vorder);
//					printf("order = %d\n", order);
//					printf("face %d, ori = %d, order = (%d %d)\n", iface, ori, horder, vorder);
					printf("ori = %d, order = (%d %d)\n", ori, horder, vorder);

					int *face_fn[] = {
						shapeset->get_face_indices(0, ori, order),
						shapeset->get_face_indices(5, ori, order)
					};
					int fn_idx = shapeset->get_num_face_fns(order) - 1;

							int np = 5;								// number of points in one direction (matrix of points np x np)
							double x, y;								// position on the part

//							printf("\n");
							double h = 2.0 / np;
							for (int i = 0; i <= np; i++) {
								double y = (h * i) - 1.0;				// position on the interval [-1, 1]
								for (int j = 0; j <= np; j++) {
									double x = (h * j) - 1.0;			// position on the interval [-1, 1]

									// constraining value
									double val[] = {
										shapeset->get_fn_value(face_fn[0][fn_idx], -1.0, x, y, 0),
										shapeset->get_fn_value(face_fn[1][fn_idx], x, y, 1.0, 0)
									};
//									printf("% lf, ", ced_val);
//									printf("% lf == % lf (% lf), ", val[0], val[1], fabs(val[0] - val[1]));
									printf("% lf == % lf, ", val[0], val[1]);
								}
								printf("\n");
							}

//							printf(".");
//					printf("\n");
				}
			}

	return false;
}

int main(int argc, char *argv[]) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	H1ShapesetLobattoHex shapeset;

	Mesh mesh;

	// combine all possible situations how two hexahedrons can face each other
//	for (int i = 0; i < 48; i++) {
//		for (int j = 0; j < 48; j++) {

	int i = 0;
	int j = 0;
			// build the mesh
			Array<Vertex *> vertices;
			Array<Element *> elements;
			Array<Boundary *> boundaries;
			MapOrd<Facet *> facets;

			for (int k = 0; k < countof(vtcs); k++)
				vertices.add(new Vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z));

			// build the mesh
			Hex *h1 = mesh.add_hex(hex[0] + i);
			Hex *h2 = mesh.add_hex(hex[1] + j);

			mesh.ugh();

			// test
			printf("pos = %d, %d\n", i, j);

//			res = test_edge_ced(&mesh, &shapeset) ? ERR_SUCCESS : ERR_FAILURE;
//			res = test_edge_face_ced(&mesh, &shapeset) ? ERR_SUCCESS : ERR_FAILURE;
			res = test_face_ced(&mesh, &shapeset) ? ERR_SUCCESS : ERR_FAILURE;
//			res = test_face_ced_0(&mesh, &shapeset) ? ERR_SUCCESS : ERR_FAILURE;
//			res = test_match(&shapeset);
//		}
//	}


#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return res;
}
