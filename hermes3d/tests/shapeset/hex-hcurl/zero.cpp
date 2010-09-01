// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
// - Pavel Kus
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
// testing zero function values on edges and faces
//

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>


// check edge functions
bool test_zero_values_of_edge_fns(Shapeset *shapeset) {
	_F_
	const int num_edges = 11;
	// indexing[edge] => { edges to check where the function is zero (local indices) }
	int edges[][num_edges] = {
		{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10 },
	};
	const int num_faces = 4;
	// indexing[edge] => { faces to check where the function is zero (local indices) }
	int faces[][num_faces] = {
		{ 0, 1, 3, 5 },
		{ 0, 2, 3, 5 },
		{ 0, 1, 2, 5 },
		{ 1, 2, 3, 5 },
		{ 1, 3, 4, 5 },
		{ 0, 3, 4, 5 },
		{ 0, 2, 4, 5 },
		{ 1, 2, 4, 5 },
		{ 0, 1, 3, 4 },
		{ 0, 2, 3, 4 },
		{ 0, 1, 2, 4 },
		{ 1, 2, 3, 4 },
	};

	Quad3D *quad = get_quadrature(MODE);
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		for (int ori = 0; ori < 2; ori++) {
			order1_t order = H3D_MAX_ELEMENT_ORDER;
			int n_fns = shapeset->get_num_edge_fns(order);
			int *edge_fn = shapeset->get_edge_indices(iedge, ori, order);

			for (int fn = 0; fn < n_fns; fn++) {
				printf("  * Edge fn #%d (edge = %d, ori = %d)", edge_fn[fn], iedge, ori);

				// edges
				for (int i = 0; i < num_edges; i++) {
					int max_order = quad->get_edge_max_order(edges[iedge][i]);
					QuadPt3D *pts = quad->get_edge_points(edges[iedge][i], max_order);
					for (int j = 0; j < quad->get_edge_num_points(iedge, max_order); j++) {
						int comp = RefHex::get_edge_tangent_direction(edges[iedge][i]);
						if (shapeset->get_fn_value(edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
							warning("Edge fn #%d is not zero at (% lf, %lf, %lf), edge %d, component %d.\n", edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, edges[iedge][i], comp);
							return false;
						}
					}
				}

				// faces
				for (int i = 0; i < num_faces; i++) {
					order2_t max_order = quad->get_face_max_order(faces[iedge][i]);
					QuadPt3D *pts = quad->get_face_points(faces[iedge][i], max_order);
					for (int j = 0; j < quad->get_face_num_points(faces[iedge][i], max_order); j++) {
						for (int icomp = 0; icomp < 2; icomp++) {
							int comp = RefHex::get_face_tangent_direction(faces[iedge][i], icomp);
							if (shapeset->get_fn_value(edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
								warning("Edge fn #%d is not zero at (% lf, %lf, %lf), face %d.", edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, faces[iedge][i]);
								return false;
							}
						}
					}
				}

				printf("... ok\n");
			}
		}
	}

	return true;
}

// check face functions
bool test_zero_values_of_face_fns(Shapeset *shapeset) {
	_F_
	const int num_edges = 12;
	// indexing[face] => { edges to check where the function is zero (local indices) }
	int edges[][num_edges] = {
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
	};
	const int num_faces = 5;
	// indexing[face] => { faces to check where the function is zero (local indices) }
	int faces[][num_faces] = {
		{ 1, 2, 3, 4, 5 },
		{ 0, 2, 3, 4, 5 },
		{ 0, 1, 3, 4, 5 },
		{ 0, 1, 2, 4, 5 },
		{ 0, 1, 2, 3, 5 },
		{ 0, 1, 2, 3, 4 }
	};

	Quad3D *quad = get_quadrature(MODE);
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		for (int ori = 0; ori < 8; ori++) {
			order2_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);

			int n_fns = shapeset->get_num_face_fns(order);
			int *face_fn = shapeset->get_face_indices(iface, ori, order);

			for (int fn = 0; fn < n_fns; fn++) {
				printf("  * Face fn #%d (face = %d, ori = %d) ", face_fn[fn], iface, ori);

				// edges
				for (int i = 0; i < num_edges; i++) {
					int max_order = quad->get_edge_max_order(edges[iface][i]);
					QuadPt3D *pts = quad->get_edge_points(edges[iface][i], max_order);
					for (int j = 0; j < quad->get_edge_num_points(edges[iface][i], max_order); j++) {
						int comp = RefHex::get_edge_tangent_direction(edges[iface][i]);
						if (shapeset->get_fn_value(face_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
							warning("Face fn #%d is not zero at (% lf, %lf, %lf), edge %d.", face_fn[fn], pts[j].x, pts[j].y, pts[j].z, edges[iface][i]);
							return false;
						}
					}
				}

				// faces
				for (int i = 0; i < num_faces; i++) {
					order2_t max_order = quad->get_face_max_order(faces[iface][i]);
					QuadPt3D *pts = quad->get_face_points(faces[iface][i], max_order);
					for (int j = 0; j < quad->get_face_num_points(faces[iface][i], max_order); j++) {
						for (int icomp = 0; icomp < 2; icomp++) {
							int comp = RefHex::get_face_tangent_direction(faces[iface][i], icomp);
							if (shapeset->get_fn_value(face_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
								warning("Face fn #%d is not zero at (% lf, %lf, %lf), face %d.", face_fn[fn], pts[j].x, pts[j].y, pts[j].z, faces[iface][i]);
								return false;
							}
						}
					}
				}

				printf("... ok\n");
			}
		}
	}

	return true;
}

// check bubble functions
bool test_zero_values_of_bubble_fns(Shapeset *shapeset) {
	_F_
	Quad3D *quad = get_quadrature(MODE);

	order3_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);

	int n_fns = shapeset->get_num_bubble_fns(order);
	int *bubble_fn = shapeset->get_bubble_indices(order);

	for (int fn = 0; fn < n_fns; fn++) {
		printf("  * Bubble fn #%d ", bubble_fn[fn]);

		// edges
		for (int i = 0; i < Hex::NUM_EDGES; i++) {
			int max_order = quad->get_edge_max_order(i);
			QuadPt3D *pts = quad->get_edge_points(i, max_order);
			for (int j = 0; j < quad->get_edge_num_points(i, max_order); j++) {
				int comp = RefHex::get_edge_tangent_direction(i);
				if (shapeset->get_fn_value(bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
					warning("Bubble fn #%d is not zero at (% lf, %lf, %lf), edge %d.", bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, i);
					return false;
				}
			}
		}

		// faces
		for (int i = 0; i < Hex::NUM_FACES; i++) {
			order2_t max_order = quad->get_face_max_order(i);
			QuadPt3D *pts = quad->get_face_points(i, max_order);
			for (int j = 0; j < quad->get_face_num_points(i, max_order); j++) {
				for(int i_comp = 0; i_comp < 2; i_comp++){
					int comp = RefHex::get_face_tangent_direction(i, i_comp);
					if (shapeset->get_fn_value(bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, comp) > EPS) {
						warning("Bubble fn #%d is not zero at (% lf, %lf, %lf), face %d.", bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, i);
						return false;
					}
				}
			}
		}

		printf("... ok\n");
	}

	return true;
}

bool test_zero_values(Shapeset *shapeset) {
	_F_
	printf("II. function values\n");

	if (!test_zero_values_of_edge_fns(shapeset))
		return false;

	if (!test_zero_values_of_face_fns(shapeset))
		return false;

	if (!test_zero_values_of_bubble_fns(shapeset))
		return false;

	return true;
}
