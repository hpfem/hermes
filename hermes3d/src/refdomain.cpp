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

#include "refdomain.h"
#include "h3d_common.h"
#include "../../hermes_common/error.h"

//
// 1D domains
//

Point1D RefLine_vertices[] = {
	{ -1.0 },
	{  1.0 }
};


//
// 2D domains
//


// Triangle ///////////////////////////////////////////////////////////////////

Point2D RefTri_vertices[] = {
	{ -1.0, -1.0 },
	{  1.0, -1.0 },
	{ -1.0,  1.0 }
};

int2 RefTri_edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };


// Quad ///////////////////////////////////////////////////////////////////////

Point2D RefQuad_vertices[] = {
	{ -1.0, -1.0 },
	{  1.0, -1.0 },
	{  1.0,  1.0 },
	{ -1.0,  1.0 }
};

int2 RefQuad_edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };


//
// 3D domains
//


// Tetra //////////////////////////////////////////////////////////////////////

Point3D RefTetra_vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 }
};

int2 RefTetra_edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
int tet_face_vtcs_0[] = { 0, 1, 3 };
int tet_face_vtcs_1[] = { 1, 2, 3 };
int tet_face_vtcs_2[] = { 2, 0, 3 };
int tet_face_vtcs_3[] = { 0, 2, 1 };
int *RefTetra_face_vtcs[] = { tet_face_vtcs_0, tet_face_vtcs_1, tet_face_vtcs_2, tet_face_vtcs_3 };
int RefTetra_face_nvtcs[] = { 3, 3, 3, 3 };
int tet_face_edges_0[] = { 0, 4, 3 };
int tet_face_edges_1[] = { 1, 5, 4 };
int tet_face_edges_2[] = { 2, 3, 5 };
int tet_face_edges_3[] = { 0, 2, 1 };
int *RefTetra_face_edges[] = { tet_face_edges_0, tet_face_edges_1, tet_face_edges_2, tet_face_edges_3 };
int RefTetra_face_nedges[] = { 3, 3, 3, 3 };
EMode2D RefTetra_face_mode[] = { MODE_TRIANGLE, MODE_TRIANGLE, MODE_TRIANGLE, MODE_TRIANGLE };
int RefTetra_face_orientations[] = { 6, 6, 6, 6 };
Point3D RefTetra_face_normal[] = { { 0, -1, 0 }, { 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) }, { -1, 0, 0 }, { 0, 0, -1 } };

// Hex ////////////////////////////////////////////////////////////////////////

Point3D RefHex_vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{  1.0,  1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 },
	{  1.0, -1.0,  1.0 },
	{  1.0,  1.0,  1.0 },
	{ -1.0,  1.0,  1.0 }
};

int2 RefHex_edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 3, 2 }, { 0, 3 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 7, 6 }, { 4, 7 } };

int hex_face_vtcs_0[] = { 0, 3, 7, 4 };
int hex_face_vtcs_1[] = { 1, 2, 6, 5 };
int hex_face_vtcs_2[] = { 0, 1, 5, 4 };
int hex_face_vtcs_3[] = { 3, 2, 6, 7 };
int hex_face_vtcs_4[] = { 0, 1, 2, 3 };
int hex_face_vtcs_5[] = { 4, 5, 6, 7 };
int *RefHex_face_vtcs[] = { hex_face_vtcs_0, hex_face_vtcs_1, hex_face_vtcs_2, hex_face_vtcs_3, hex_face_vtcs_4, hex_face_vtcs_5 };
int RefHex_face_nvtcs[] = { 4, 4, 4, 4, 4, 4 };

int hex_face_edges_0[] = { 3, 7, 11, 4 };
int hex_face_edges_1[] = { 1, 6, 9, 5 };
int hex_face_edges_2[] = { 0, 5, 8, 4 };
int hex_face_edges_3[] = { 2, 6, 10, 7 };
int hex_face_edges_4[] = { 0, 1, 2, 3 };
int hex_face_edges_5[] = { 8, 9, 10, 11 };
int *RefHex_face_edges[] = { hex_face_edges_0, hex_face_edges_1, hex_face_edges_2, hex_face_edges_3, hex_face_edges_4, hex_face_edges_5 };
int RefHex_face_nedges[] = { 4, 4, 4, 4, 4, 4 };
EMode2D RefHex_face_mode[] = { MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD };
int RefHex_face_orientations[] = { 8, 8, 8, 8, 8, 8 };
int RefHex_face_edge_ori[8][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 }, { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
int RefHex_edge_tangent[12] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };
int RefHex_face_tangent[6][2] = { { 1, 2 }, { 1, 2 }, { 0, 2 }, { 0, 2 }, { 0, 1 }, { 0, 1 } };


// Prism //////////////////////////////////////////////////////////////////////

Point3D RefPrism_vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 },
	{  1.0, -1.0,  1.0 },
	{ -1.0,  1.0,  1.0 }
};

int2 RefPrism_edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3}, { 1, 4 }, { 2, 5 }, { 3, 4 }, { 4, 5 }, { 5, 3 } };
int std_3d_prism_face_vtcs_0[] = { 0, 1, 4, 3 };
int std_3d_prism_face_vtcs_1[] = { 1, 2, 5, 4 };
int std_3d_prism_face_vtcs_2[] = { 2, 0, 3, 5 };
int std_3d_prism_face_vtcs_3[] = { 0, 2, 1 };
int std_3d_prism_face_vtcs_4[] = { 3, 4, 5 };
int *RefPrism_face_vtcs[] = { std_3d_prism_face_vtcs_0, std_3d_prism_face_vtcs_1, std_3d_prism_face_vtcs_2, std_3d_prism_face_vtcs_3, std_3d_prism_face_vtcs_4 };
int RefPrism_face_nvtcs[] = { 4, 4, 4, 3, 3 };
int std_3d_prism_face_edges_0[] = { 0, 4, 6, 3 };
int std_3d_prism_face_edges_1[] = { 1, 5, 7, 4 };
int std_3d_prism_face_edges_2[] = { 2, 3, 8, 5 };
int std_3d_prism_face_edges_3[] = { 0, 2, 1 };
int std_3d_prism_face_edges_4[] = { 6, 7, 8 };
int *RefPrism_face_edges[] = { std_3d_prism_face_edges_0, std_3d_prism_face_edges_1, std_3d_prism_face_edges_2, std_3d_prism_face_edges_3, std_3d_prism_face_edges_4 };
int RefPrism_face_nedges[] = { 4, 4, 4, 3, 3 };
EMode2D RefPrism_face_mode[] = { MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_TRIANGLE, MODE_TRIANGLE };
int RefPrism_face_orientations[] = { 8, 8, 8, 6, 6 };

