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

//
// testing continuity of function values over faces
//

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

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
	{ -2, -1, -1 },
	{  0, -1, -1 },
	{  0,  1, -1 },
	{ -2,  1, -1 },
	{ -2, -1,  1 },
	{  0, -1,  1 },
	{  0,  1,  1 },
	{ -2,  1,  1 },
	{  2, -1, -1 },
	{  2,  1, -1 },
	{  2, -1,  1 },
	{  2,  1,  1 },
};

// mesh
Word_t hexs[2][48][8] = {
	// hexs 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },
		{  7,  6,  5,  4,  3,  2,  1,  0 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  3,  0,  1,  2,  7,  4,  5,  6 },
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  0,  4,  5,  1,  3,  7,  6,  2 }
	},
	// hexs 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },
		{  6, 11, 10,  5,  2,  9,  8,  1 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  2,  1,  8,  9,  6,  5, 10, 11 },
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5,  10 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{  1,  5, 10,  8,  2,  6, 11,  9 }
	}
};

Word_t bnd[10][5] = {
	{ 0,  3,  7,  4, 1 },
	{ 8,  9, 11, 10, 2 },
	{ 0,  1,  5,  4, 3 },
	{ 1,  8, 10,  5, 3 },
	{ 3,  2,  6,  7, 4 },
	{ 2,  9, 11,  6, 4 },
	{ 0,  1,  2,  3, 5 },
	{ 1,  8,  9,  2, 5 },
	{ 4,  5,  6,  7, 6 },
	{ 5, 10, 11,  6, 6 }
};

struct Pts {
	QuadPt3D *ref_pt;
	double *phys_x;
	double *phys_y;
	double *phys_z;
	int *sort;
};

static Pts *sort_pts;

int compare(const void *a, const void *b) {
	int i1 = *(int *) a;
	int i2 = *(int *) b;

	double val1 = 1000000. * sort_pts->phys_x[i1] + 1000. * sort_pts->phys_y[i1] + sort_pts->phys_z[i1];
	double val2 = 1000000. * sort_pts->phys_x[i2] + 1000. * sort_pts->phys_y[i2] + sort_pts->phys_z[i2];

	if (val1 < val2) return -1;
	else if (val1 > val2) return 1;
	else return 0;
}


bool test_cont_values_of_vertex_fns(Mesh *mesh, Word_t fid, int pos0, int pos1, Shapeset *shapeset) {
	_F_
	Facet *facet = mesh->facets[fid];
	Quad3D *quad = get_quadrature(MODE);

	Element *h[2] = {
		mesh->elements[facet->left],
		mesh->elements[facet->right]
	};

	// face vertices
	const int *face_vtx[] = {
		h[0]->get_face_vertices(facet->left_face_num),
		h[1]->get_face_vertices(facet->right_face_num)
	};

	RefMap rm(mesh);
	Pts vpts[2];

	// get coordinates of vertices on the face
	int vnp = Quad::NUM_VERTICES;
	for (int i = 0; i < 2; i++) {
		QuadPt3D *rd_vtx_pt = quad->get_vertex_points();

		QuadPt3D *pt = new QuadPt3D[Quad::NUM_VERTICES];
		for (int j = 0; j < vnp; j++)
			pt[j] = rd_vtx_pt[face_vtx[i][j]];

		rm.set_active_element(h[i]);

		vpts[i].ref_pt = pt;
		vpts[i].phys_x = rm.get_phys_x(vnp, pt);
		vpts[i].phys_y = rm.get_phys_y(vnp, pt);
		vpts[i].phys_z = rm.get_phys_z(vnp, pt);

		// sort points by physical coordinates
		vpts[i].sort = new int [vnp];
		for (int j = 0; j < vnp; j++)
			vpts[i].sort[j] = j;
		sort_pts = vpts + i;
		qsort(vpts[i].sort, vnp, sizeof(int), compare);
	}

	// we know which vertices correspond to each other, since we sorted them by their phys. coordinates
	// we use it for comparing matching vertex functions

	// obtain coordinates where vertex functions are going to be evaluated
	Pts fpts[2];
	order2_t face_ord = quad->get_face_max_order(facet->left_face_num);
	int np = quad->get_face_num_points(facet->left_face_num, face_ord);

	QuadPt3D *pt[2] = {
		quad->get_face_points(facet->left_face_num, face_ord),
		quad->get_face_points(facet->right_face_num, face_ord)
	};

	// transform the coorinates so we can compare them
	for (int i = 0; i < 2; i++) {
		rm.set_active_element(h[i]);

		fpts[i].ref_pt = pt[i];
		fpts[i].phys_x = rm.get_phys_x(np, pt[i]);
		fpts[i].phys_y = rm.get_phys_y(np, pt[i]);
		fpts[i].phys_z = rm.get_phys_z(np, pt[i]);

		// sort points by physical coordinates
		fpts[i].sort = new int [np];
		for (int j = 0; j < np; j++)
			fpts[i].sort[j] = j;
		sort_pts = fpts + i;
		qsort(fpts[i].sort, np, sizeof(int), compare);
	}

	// loop through vertex functions
	for (int k = 0; k < Quad::NUM_VERTICES; k++) {
		int i1 = vpts[0].sort[k];
		int i2 = vpts[1].sort[k];

		// get values of vertex functions
		double uval[np], vval[np];
		shapeset->get_fn_values(shapeset->get_vertex_index(face_vtx[0][i1]), np, pt[0], 0, uval);
		shapeset->get_fn_values(shapeset->get_vertex_index(face_vtx[1][i2]), np, pt[1], 0, vval);

		// compare their values
		for (int i = 0; i < np; i++) {
			int j1 = fpts[0].sort[i];
			int j2 = fpts[1].sort[i];

			if (fabs(uval[j1] - vval[j2]) > EPS) {
				printf("failed\n");
				warning("Vertex fn not continuous @(% lf, % lf, % lf), diff = % e",
					fpts[0].ref_pt[j1].x, fpts[0].ref_pt[j1].y, fpts[0].ref_pt[j1].z,
					fabs(uval[j1] - vval[j2]));
				return false;
			}
		}
	}

	// free memory
	for (int i = 0; i < 2; i++) {
		delete [] vpts[i].ref_pt;
		delete [] vpts[i].phys_x;
		delete [] vpts[i].phys_y;
		delete [] vpts[i].phys_z;
		delete [] vpts[i].sort;

		delete [] fpts[i].phys_x;
		delete [] fpts[i].phys_y;
		delete [] fpts[i].phys_z;
		delete [] fpts[i].sort;
	}

	return true;
}

bool test_cont_values_of_edge_fns(Mesh *mesh, Word_t fid, int pos0, int pos1, Shapeset *shapeset) {
	_F_
	Facet *facet = mesh->facets[fid];
	Quad3D *quad = get_quadrature(MODE);

	Element *h[2] = {
		mesh->elements[facet->left],
		mesh->elements[facet->right]
	};

	// local number of edges on the face
	const int *face_edges[] = {
		h[0]->get_face_edges(facet->left_face_num),
		h[1]->get_face_edges(facet->right_face_num)
	};

	// orientations of edges on the face
	int edge_ori[2][4];
	for (int i = 0; i < 4; i++) {
		edge_ori[0][i] = h[0]->get_edge_orientation(face_edges[0][i]);
		edge_ori[1][i] = h[1]->get_edge_orientation(face_edges[1][i]);
	}

	RefMap rm(mesh);
	Pts epts[2];

	// get coordinates of the point in the middle of an edge on the face
	// it is used to find matching edges
	int enp = Quad::NUM_EDGES;
	for (int i = 0; i < 2; i++) {
		QuadPt3D *pt = new QuadPt3D[Quad::NUM_EDGES];
		for (int j = 0; j < enp; j++) {
			const int *edge_vtx = RefHex::get_edge_vertices(face_edges[i][j]);
			QuadPt3D *rd_vtx_pt = quad->get_vertex_points();
			pt[j].x = (rd_vtx_pt[edge_vtx[0]].x + rd_vtx_pt[edge_vtx[1]].x) / 2.0;
			pt[j].y = (rd_vtx_pt[edge_vtx[0]].y + rd_vtx_pt[edge_vtx[1]].y) / 2.0;
			pt[j].z = (rd_vtx_pt[edge_vtx[0]].z + rd_vtx_pt[edge_vtx[1]].z) / 2.0;
		}

		rm.set_active_element(h[i]);

		epts[i].ref_pt = pt;
		epts[i].phys_x = rm.get_phys_x(enp, pt);
		epts[i].phys_y = rm.get_phys_y(enp, pt);
		epts[i].phys_z = rm.get_phys_z(enp, pt);

		// sort points by physical coordinates
		epts[i].sort = new int [enp];
		for (int j = 0; j < enp; j++)
			epts[i].sort[j] = j;
		sort_pts = epts + i;
		qsort(epts[i].sort, enp, sizeof(int), compare);
	}

	// we know which edge correspond to each other, since we sorted them by their phys. coordinates
	// we use it for comparing matching edge functions

	// obtain coordinates where edge functions are going to be evaluated
	Pts fpts[2];
	order2_t face_ord = quad->get_face_max_order(facet->left_face_num);
	int np = quad->get_face_num_points(facet->left_face_num, face_ord);

	QuadPt3D *pt[2] = {
		quad->get_face_points(facet->left_face_num, face_ord),
		quad->get_face_points(facet->right_face_num, face_ord)
	};

	// transform the coorinates so we can compare them
	for (int i = 0; i < 2; i++) {
		rm.set_active_element(h[i]);

		fpts[i].ref_pt = pt[i];
		fpts[i].phys_x = rm.get_phys_x(np, pt[i]);
		fpts[i].phys_y = rm.get_phys_y(np, pt[i]);
		fpts[i].phys_z = rm.get_phys_z(np, pt[i]);

		// sort points by physical coordinates
		fpts[i].sort = new int [np];
		for (int j = 0; j < np; j++)
			fpts[i].sort[j] = j;
		sort_pts = fpts + i;
		qsort(fpts[i].sort, np, sizeof(int), compare);
	}

	// loop through edges
	for (int k = 0; k < Quad::NUM_EDGES; k++) {
		int i1 = epts[0].sort[k];
		int i2 = epts[1].sort[k];

		// loop through edge functions
		// get all face functions on the face
		order1_t order(H3D_MAX_ELEMENT_ORDER);
		int *edge_fn[] = {
			shapeset->get_edge_indices(face_edges[0][i1], edge_ori[0][i1], order),
			shapeset->get_edge_indices(face_edges[1][i2], edge_ori[1][i2], order)
		};

		// loop through face functions
		int n_edge_fns = shapeset->get_num_edge_fns(order);
		for (int ifn = 0; ifn < n_edge_fns; ifn++) {
			// get values of edge functions
			double uval[np], vval[np];
			shapeset->get_fn_values(edge_fn[0][ifn], np, pt[0], 0, uval);
			shapeset->get_fn_values(edge_fn[1][ifn], np, pt[1], 0, vval);

			// compare their values
			for (int i = 0; i < np; i++) {
				int j1 = fpts[0].sort[i];
				int j2 = fpts[1].sort[i];

				if (fabs(uval[j1] - vval[j2]) > EPS) {
					printf("failed\n");
					warning("Edge fn not continuous @(% lf, % lf, % lf), diff = % e",
						fpts[0].ref_pt[j1].x, fpts[0].ref_pt[j1].y, fpts[0].ref_pt[j1].z,
						fabs(uval[j1] - vval[j2]));
					return false;
				}
			}
		}
	}

	// free memory
	for (int i = 0; i < 2; i++) {
		delete [] epts[i].ref_pt;
		delete [] epts[i].phys_x;
		delete [] epts[i].phys_y;
		delete [] epts[i].phys_z;
		delete [] epts[i].sort;

		delete [] fpts[i].phys_x;
		delete [] fpts[i].phys_y;
		delete [] fpts[i].phys_z;
		delete [] fpts[i].sort;
	}

	return true;
}

bool test_cont_values_of_face_fns(Mesh *mesh, Word_t fid, int pos0, int pos1, Shapeset *shapeset) {
	_F_
	Quad3D *quad = get_quadrature(MODE);
	Facet *facet = mesh->facets[fid];

	Element *h[2] = {
		mesh->elements[facet->left],
		mesh->elements[facet->right]
	};

	int face_ori[2] = {
		h[0]->get_face_orientation(facet->left_face_num),
		h[1]->get_face_orientation(facet->right_face_num)
	};

	RefMap rm(mesh);

	// get points where the functions are going to be evaluated
	Pts fpts[2];
	order2_t face_ord = quad->get_face_max_order(facet->left_face_num);
	int np = quad->get_face_num_points(facet->left_face_num, face_ord);

	QuadPt3D *pt[2] = {
		quad->get_face_points(facet->left_face_num, face_ord),
		quad->get_face_points(facet->right_face_num, face_ord)
	};

	// transform points and sort them
	for (int i = 0; i < 2; i++) {
		rm.set_active_element(h[i]);

		fpts[i].ref_pt = pt[i];
		fpts[i].phys_x = rm.get_phys_x(np, pt[i]);
		fpts[i].phys_y = rm.get_phys_y(np, pt[i]);
		fpts[i].phys_z = rm.get_phys_z(np, pt[i]);

		// sort points by physical coordinates
		fpts[i].sort = new int [np];
		for (int j = 0; j < np; j++)
			fpts[i].sort[j] = j;
		sort_pts = fpts + i;
		qsort(fpts[i].sort, np, sizeof(int), compare);
	}

	// get all face functions on the face
	order2_t order(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
	int *face_fn[] = {
		shapeset->get_face_indices(facet->left_face_num, face_ori[0], order),
		shapeset->get_face_indices(facet->right_face_num, face_ori[1], order)
	};

	// loop through face functions
	int n_face_fns = shapeset->get_num_face_fns(order);
	for (int ifn = 0; ifn < n_face_fns; ifn++) {
		// evaluate functions
		double uval[np], vval[np];
		shapeset->get_fn_values(face_fn[0][ifn], np, pt[0], 0, uval);
		shapeset->get_fn_values(face_fn[1][ifn], np, pt[1], 0, vval);

		// compare
		for (int i = 0; i < np; i++) {
			int j1 = fpts[0].sort[i];
			int j2 = fpts[1].sort[i];

			if (fabs(uval[j1] - vval[j2]) > EPS) {
				printf("failed\n");
				warning("Face fn not continuous @(% lf, % lf, % lf), diff = % e",
					fpts[0].ref_pt[j1].x, fpts[0].ref_pt[j1].y, fpts[0].ref_pt[j1].z,
					fabs(uval[j1] - vval[j2]));
				return false;
			}
		}
	}

	for (int i = 0; i < 2; i++) {
		delete [] fpts[i].phys_x;
		delete [] fpts[i].phys_y;
		delete [] fpts[i].phys_z;
		delete [] fpts[i].sort;
	}

	return true;
}

//
// main check function
//
bool test_continuity(Shapeset *shapeset) {
	_F_
	printf("III. continuity\n");

	// combine all possible situations how two hexahedrons can face each other
	for (int i = 0; i < 48; i++) {
		for (int j = 0; j < 48; j++) {
			// build the mesh
			Mesh mesh;
			for (Word_t k = 0; k < countof(vtcs); k++)
				mesh.add_vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z);
			Word_t h1[] = {
					hexs[0][i][0] + 1, hexs[0][i][1] + 1, hexs[0][i][2] + 1, hexs[0][i][3] + 1,
					hexs[0][i][4] + 1, hexs[0][i][5] + 1, hexs[0][i][6] + 1, hexs[0][i][7] + 1 };
			mesh.add_hex(h1);
			Word_t h2[] = {
					hexs[1][j][0] + 1, hexs[1][j][1] + 1, hexs[1][j][2] + 1, hexs[1][j][3] + 1,
					hexs[1][j][4] + 1, hexs[1][j][5] + 1, hexs[1][j][6] + 1, hexs[1][j][7] + 1 };
			mesh.add_hex(h2);
			// bc
			for (Word_t k = 0; k < countof(bnd); k++) {
				Word_t facet_idxs[Quad::NUM_VERTICES] = { bnd[k][0] + 1, bnd[k][1] + 1, bnd[k][2] + 1, bnd[k][3] + 1 };
				mesh.add_quad_boundary(facet_idxs, bnd[k][4]);
			}
			mesh.ugh();

			printf("pos = %d, %d\n", i, j);

			// test continuity on inner factes
			// since we have only 2 elements, there is only one such facet
			for (Word_t fid = mesh.facets.first(); fid != INVALID_IDX; fid = mesh.facets.next(fid)) {
				Facet *facet = mesh.facets[fid];
				if (facet->type == Facet::INNER) {
					printf("  - vertex fns..."); fflush(stdout);
					if (!test_cont_values_of_vertex_fns(&mesh, fid, i, j, shapeset))
						return false;
					printf("ok\n");

					printf("  - edge fns..."); fflush(stdout);
					if (!test_cont_values_of_edge_fns(&mesh, fid, i, j, shapeset))
						return false;
					printf("ok\n");

					printf("  - face fns..."); fflush(stdout);
					if (!test_cont_values_of_face_fns(&mesh, fid, i, j, shapeset))
						return false;
					printf("ok\n");
				}
			}
		}
	}

	return true;
}
