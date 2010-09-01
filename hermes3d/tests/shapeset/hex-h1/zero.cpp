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
 * zero.cc
 *
 * testing zero functions values at vertices, on edges and faces
 *
 */

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// check vertex functions
bool test_zero_values_of_vertex_fns(Shapeset *shapeset) {
	_F_
	const int num_vertices = 7;
	// indexing[vertex] => { vertices to check where the function is zero (local indices) }
	int vertices[][num_vertices] = {
		{ 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 7 },
		{ 0, 1, 2, 3, 4, 5, 6 }
	};
	const int num_edges = 9;
	// indexing[vertex] => { edges to check where the function is zero (local indices) }
	int edges[][num_edges] = {
		{ 1, 2, 5, 6, 7, 8, 9, 10, 11 },
		{ 2, 3, 4, 6, 7, 8, 9, 10, 11 },
		{ 0, 3, 4, 5, 7, 8, 9, 10, 11 },
		{ 0, 1, 4, 5, 6, 8, 9, 10, 11 },
		{ 0, 1, 2, 3, 5, 6, 7,  9, 10 },
		{ 0, 1, 2, 3, 4, 6, 7, 10, 11 },
		{ 0, 1, 2, 3, 4, 5, 7,  8, 11 },
		{ 0, 1, 2, 3, 4, 5, 6,  8,  9 },
	};
	const int num_faces = 3;
	// indexing[vertex] => { faces to check where the function is zero (local indices) }
	int faces[][num_faces] = {
		{ 5, 1, 3 },
		{ 5, 3, 0 },
		{ 5, 0, 2 },
		{ 5, 2, 1 },
		{ 4, 1, 3 },
		{ 3, 0, 4 },
		{ 0, 2, 4 },
		{ 1, 4, 2 }
	};

	Quad3D *quad = get_quadrature(MODE);
	for (int vtx = 0; vtx < Hex::NUM_VERTICES; vtx++) {
		int fn_idx = shapeset->get_vertex_index(vtx);
		printf("  * Vertex fn #%d (%d) ", vtx, fn_idx);

		// vertices
		const Point3D *vtx_pt = REF_DOMAIN::get_vertices();
		int *idx = vertices[vtx];
		for (int i = 0; i < num_vertices; i++) {
			if (shapeset->get_fn_value(fn_idx, vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, 0) > EPS) {
				warning("Vertex fn #%d (%d) is not zero at (% lf, %lf, %lf), vertex #%d.", vtx, fn_idx, vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, idx[i]);
				return false;
			}
		}

		// edges
		for (int i = 0; i < num_edges; i++) {
			order1_t max_order = quad->get_edge_max_order(edges[vtx][i]);
			int np = quad->get_edge_num_points(edges[vtx][i], max_order);
			QuadPt3D *pts = quad->get_edge_points(edges[vtx][i], max_order);

			double vals[np];
			shapeset->get_fn_values(fn_idx, np, pts, 0, vals);
			for (int j = 0; j < np; j++) {
				if (vals[j] > EPS) {
					warning("Vertex fn #%d (%d) is not zero at (% lf, %lf, %lf), edge %d.", vtx, fn_idx, pts[j].x, pts[j].y, pts[j].z, edges[vtx][i]);
					return false;
				}
			}
		}

		// faces
		for (int i = 0; i < num_faces; i++) {
			order2_t max_order = quad->get_face_max_order(faces[vtx][i]);
			int np = quad->get_face_num_points(faces[vtx][i], max_order);
			QuadPt3D *pts = quad->get_face_points(faces[vtx][i], max_order);
			double vals[np];
			shapeset->get_fn_values(fn_idx, np, pts, 0, vals);
			for (int j = 0; j < np; j++) {
				if (vals[j] > EPS) {
					warning("Vertex fn #%d (%d) is not zero at (% lf, %lf, %lf), face %d.", vtx, fn_idx, pts[j].x, pts[j].y, pts[j].z, faces[vtx][i]);
					return false;
				}
			}
		}

		printf("... ok\n");
	}

	return true;
}

// check edge functions
bool test_zero_values_of_edge_fns(Shapeset *shapeset) {
	_F_
	const int num_vertices = 8;
	// indexing[edge] => { vertices to check where the function is zero (local indices) }
	int vertices[][num_vertices] = {
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 }
	};
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
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++) {
		for (int ori = 0; ori < 2; ori++) {
			int order = H3D_MAX_ELEMENT_ORDER;
			int n_fns = shapeset->get_num_edge_fns(order);
			int *edge_fn = shapeset->get_edge_indices(edge, ori, order);

			for (int fn = 0; fn < n_fns; fn++) {
				printf("  * Edge fn #%d (edge = %d, ori = %d) ", edge_fn[fn], edge, ori);

				// vertices
				const Point3D *vtx_pt = REF_DOMAIN::get_vertices();
				int *idx = vertices[edge];
				for (int i = 0; i < num_vertices; i++) {
					if (shapeset->get_fn_value(edge_fn[fn], vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, 0) > EPS) {
						warning("Edge fn #%d is not zero at (% lf, %lf, %lf), vertex #%d.", edge_fn[fn], vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, idx[i]);
						return false;
					}
				}

				// edges
				for (int i = 0; i < num_edges; i++) {
					order1_t max_order = quad->get_edge_max_order(edges[edge][i]);
					int np = quad->get_edge_num_points(edges[edge][i], max_order);
					QuadPt3D *pts = quad->get_edge_points(edges[edge][i], max_order);
					double vals[np];
					shapeset->get_fn_values(edge_fn[fn], np, pts, 0, vals);
					for (int j = 0; j < np; j++) {
						if (vals[j] > EPS) {
							warning("Edge fn #%d is not zero at (% lf, %lf, %lf), edge %d.", edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, edges[edge][i]);
							return false;
						}
					}
				}

				// faces
				for (int i = 0; i < num_faces; i++) {
					order2_t max_order = quad->get_face_max_order(faces[edge][i]);
					int np = quad->get_face_num_points(faces[edge][i], max_order);
					QuadPt3D *pts = quad->get_face_points(faces[edge][i], max_order);
					double vals[np];
					shapeset->get_fn_values(edge_fn[fn], np, pts, 0, vals);
					for (int j = 0; j < np; j++) {
						if (vals[j] > EPS) {
							warning("Edge fn #%d is not zero at (% lf, %lf, %lf), face %d.", edge_fn[fn], pts[j].x, pts[j].y, pts[j].z, faces[edge][i]);
							return false;
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
	const int num_vertices = 8;
	// indexing[face] => { vertices to check where the function is zero (local indices) }
	int vertices[][num_vertices] = {
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 },
		{ 0, 1, 2, 3, 4, 5, 6, 7 }
	};
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
	for (int face = 0; face < Hex::NUM_FACES; face++) {
		for (int ori = 0; ori < 8; ori++) {
			order2_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);

			int n_fns = shapeset->get_num_face_fns(order);
			int *face_fn = shapeset->get_face_indices(face, ori, order);

			for (int fn = 0; fn < n_fns; fn++) {
				printf("  * Face fn #%d (face = %d, ori = %d) ", face_fn[fn], face, ori);

				// vertices
				const Point3D *vtx_pt = REF_DOMAIN::get_vertices();
				int *idx = vertices[face];
				for (int i = 0; i < num_vertices; i++) {
					if (shapeset->get_fn_value(face_fn[fn], vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, 0) > EPS) {
						warning("Face fn #%d is not zero at (% lf, %lf, %lf), vertex #%d.", face_fn[fn], vtx_pt[idx[i]].x, vtx_pt[idx[i]].y, vtx_pt[idx[i]].z, idx[i]);
						return false;
					}
				}

				// edges
				for (int i = 0; i < num_edges; i++) {
					order1_t max_order = quad->get_edge_max_order(edges[face][i]);
					int np = quad->get_edge_num_points(edges[face][i], max_order);
					QuadPt3D *pts = quad->get_edge_points(edges[face][i], max_order);
					double vals[np];
					shapeset->get_fn_values(face_fn[fn], np, pts, 0, vals);
					for (int j = 0; j < np; j++) {
						if (vals[j] > EPS) {
							warning("Face fn #%d is not zero at (% lf, %lf, %lf), edge %d.", face_fn[fn], pts[j].x, pts[j].y, pts[j].z, edges[face][i]);
							return false;
						}
					}
				}

				// faces
				for (int i = 0; i < num_faces; i++) {
					order2_t max_order = quad->get_face_max_order(faces[face][i]);
					int np = quad->get_face_num_points(faces[face][i], max_order);
					QuadPt3D *pts = quad->get_face_points(faces[face][i], max_order);
					double vals[np];
					shapeset->get_fn_values(face_fn[fn], np, pts, 0, vals);
					for (int j = 0; j < np; j++) {
						if (vals[j] > EPS) {
							warning("Face fn #%d is not zero at (% lf, %lf, %lf), face %d.", face_fn[fn], pts[j].x, pts[j].y, pts[j].z, faces[face][i]);
							return false;
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

		// vertices
		const Point3D *vtx_pt = REF_DOMAIN::get_vertices();
		for (int i = 0; i < Hex::NUM_VERTICES; i++) {
			if (shapeset->get_fn_value(bubble_fn[fn], vtx_pt[i].x, vtx_pt[i].y, vtx_pt[i].z, 0) > EPS) {
				warning("Bubble fn #%d is not zero at (% lf, %lf, %lf), vertex #%d.", bubble_fn[fn], vtx_pt[i].x, vtx_pt[i].y, vtx_pt[i].z, i);
				return false;
			}
		}

		// edges
		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			order1_t max_order = quad->get_edge_max_order(iedge);
			int np = quad->get_edge_num_points(iedge, max_order);
			QuadPt3D *pts = quad->get_edge_points(iedge, max_order);
			double vals[np];
			shapeset->get_fn_values(bubble_fn[fn], np, pts, 0, vals);
			for (int j = 0; j < np; j++) {
				if (vals[j] > EPS) {
					warning("Bubble fn #%d is not zero at (% lf, %lf, %lf), edge %d.", bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, iedge);
					return false;
				}
			}
		}

		// faces
		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			order2_t max_order = quad->get_face_max_order(iface);
			int np = quad->get_face_num_points(iface, max_order);
			QuadPt3D *pts = quad->get_face_points(iface, max_order);
			double vals[np];
			shapeset->get_fn_values(bubble_fn[fn], np, pts, 0, vals);
			for (int j = 0; j < np; j++) {
				if (vals[j] > EPS) {
					warning("Bubble fn #%d is not zero at (% lf, %lf, %lf), face %d.", bubble_fn[fn], pts[j].x, pts[j].y, pts[j].z, iface);
					return false;
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

	if (!test_zero_values_of_vertex_fns(shapeset))
		return false;

	if (!test_zero_values_of_edge_fns(shapeset))
		return false;

	if (!test_zero_values_of_face_fns(shapeset))
		return false;

	if (!test_zero_values_of_bubble_fns(shapeset))
		return false;

	return true;
}
