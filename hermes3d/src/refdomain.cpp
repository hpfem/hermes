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
#include "refdomain.h"
#include "common.h"
#include <common/error.h>

//
// 1D domains
//

const Point1D RefLine::vertices[] = {
	{ -1.0 },
	{  1.0 }
};


//
// 2D domains
//


// Triangle ///////////////////////////////////////////////////////////////////

const Point2D RefTri::vertices[] = {
	{ -1.0, -1.0 },
	{  1.0, -1.0 },
	{ -1.0,  1.0 }
};

const int2 RefTri::edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };


// Quad ///////////////////////////////////////////////////////////////////////

const Point2D RefQuad::vertices[] = {
	{ -1.0, -1.0 },
	{  1.0, -1.0 },
	{  1.0,  1.0 },
	{ -1.0,  1.0 }
};

const int2 RefQuad::edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };


//
// 3D domains
//


// Tetra //////////////////////////////////////////////////////////////////////

const Point3D RefTetra::vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 }
};

const int2 RefTetra::edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
static int tet_face_vtcs_0[] = { 0, 1, 3 };
static int tet_face_vtcs_1[] = { 1, 2, 3 };
static int tet_face_vtcs_2[] = { 2, 0, 3 };
static int tet_face_vtcs_3[] = { 0, 2, 1 };
const int *RefTetra::face_vtcs[] = { tet_face_vtcs_0, tet_face_vtcs_1, tet_face_vtcs_2, tet_face_vtcs_3 };
const int RefTetra::face_nvtcs[] = { 3, 3, 3, 3 };
static int tet_face_edges_0[] = { 0, 4, 3 };
static int tet_face_edges_1[] = { 1, 5, 4 };
static int tet_face_edges_2[] = { 2, 3, 5 };
static int tet_face_edges_3[] = { 0, 2, 1 };
const int *RefTetra::face_edges[] = { tet_face_edges_0, tet_face_edges_1, tet_face_edges_2, tet_face_edges_3 };
const int RefTetra::face_nedges[] = { 3, 3, 3, 3 };
const EMode2D RefTetra::face_mode[] = { MODE_TRIANGLE, MODE_TRIANGLE, MODE_TRIANGLE, MODE_TRIANGLE };
const int RefTetra::face_orientations[] = { 6, 6, 6, 6 };
const Point3D RefTetra::face_normal[] = { { 0, -1, 0 }, { 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) }, { -1, 0, 0 }, { 0, 0, -1 } };

// Hex ////////////////////////////////////////////////////////////////////////

const Point3D RefHex::vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{  1.0,  1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 },
	{  1.0, -1.0,  1.0 },
	{  1.0,  1.0,  1.0 },
	{ -1.0,  1.0,  1.0 }
};

const int2 RefHex::edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 3, 2 }, { 0, 3 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 5 }, { 5, 6 }, { 7, 6 }, { 4, 7 } };

static int hex_face_vtcs_0[] = { 0, 3, 7, 4 };
static int hex_face_vtcs_1[] = { 1, 2, 6, 5 };
static int hex_face_vtcs_2[] = { 0, 1, 5, 4 };
static int hex_face_vtcs_3[] = { 3, 2, 6, 7 };
static int hex_face_vtcs_4[] = { 0, 1, 2, 3 };
static int hex_face_vtcs_5[] = { 4, 5, 6, 7 };
const int *RefHex::face_vtcs[] = { hex_face_vtcs_0, hex_face_vtcs_1, hex_face_vtcs_2, hex_face_vtcs_3, hex_face_vtcs_4, hex_face_vtcs_5 };
const int RefHex::face_nvtcs[] = { 4, 4, 4, 4, 4, 4 };

static int hex_face_edges_0[] = { 3, 7, 11, 4 };
static int hex_face_edges_1[] = { 1, 6, 9, 5 };
static int hex_face_edges_2[] = { 0, 5, 8, 4 };
static int hex_face_edges_3[] = { 2, 6, 10, 7 };
static int hex_face_edges_4[] = { 0, 1, 2, 3 };
static int hex_face_edges_5[] = { 8, 9, 10, 11 };
const int *RefHex::face_edges[] = { hex_face_edges_0, hex_face_edges_1, hex_face_edges_2, hex_face_edges_3, hex_face_edges_4, hex_face_edges_5 };
const int RefHex::face_nedges[] = { 4, 4, 4, 4, 4, 4 };
const EMode2D RefHex::face_mode[] = { MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_QUAD };
const int RefHex::face_orientations[] = { 8, 8, 8, 8, 8, 8 };
const int RefHex::face_edge_ori[8][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 }, { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
const int RefHex::edge_tangent[12] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };
const int RefHex::face_tangent[6][2] = { { 1, 2 }, { 1, 2 }, { 0, 2 }, { 0, 2 }, { 0, 1 }, { 0, 1 } };


// Prism //////////////////////////////////////////////////////////////////////

const Point3D RefPrism::vertices[] = {
	{ -1.0, -1.0, -1.0 },
	{  1.0, -1.0, -1.0 },
	{ -1.0,  1.0, -1.0 },
	{ -1.0, -1.0,  1.0 },
	{  1.0, -1.0,  1.0 },
	{ -1.0,  1.0,  1.0 }
};

const int2 RefPrism::edge_vtcs[] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3}, { 1, 4 }, { 2, 5 }, { 3, 4 }, { 4, 5 }, { 5, 3 } };
static int std_3d_prism_face_vtcs_0[] = { 0, 1, 4, 3 };
static int std_3d_prism_face_vtcs_1[] = { 1, 2, 5, 4 };
static int std_3d_prism_face_vtcs_2[] = { 2, 0, 3, 5 };
static int std_3d_prism_face_vtcs_3[] = { 0, 2, 1 };
static int std_3d_prism_face_vtcs_4[] = { 3, 4, 5 };
const int *RefPrism::face_vtcs[] = { std_3d_prism_face_vtcs_0, std_3d_prism_face_vtcs_1, std_3d_prism_face_vtcs_2, std_3d_prism_face_vtcs_3, std_3d_prism_face_vtcs_4 };
const int RefPrism::face_nvtcs[] = { 4, 4, 4, 3, 3 };
static int std_3d_prism_face_edges_0[] = { 0, 4, 6, 3 };
static int std_3d_prism_face_edges_1[] = { 1, 5, 7, 4 };
static int std_3d_prism_face_edges_2[] = { 2, 3, 8, 5 };
static int std_3d_prism_face_edges_3[] = { 0, 2, 1 };
static int std_3d_prism_face_edges_4[] = { 6, 7, 8 };
const int *RefPrism::face_edges[] = { std_3d_prism_face_edges_0, std_3d_prism_face_edges_1, std_3d_prism_face_edges_2, std_3d_prism_face_edges_3, std_3d_prism_face_edges_4 };
const int RefPrism::face_nedges[] = { 4, 4, 4, 3, 3 };
const EMode2D RefPrism::face_mode[] = { MODE_QUAD, MODE_QUAD, MODE_QUAD, MODE_TRIANGLE, MODE_TRIANGLE };
const int RefPrism::face_orientations[] = { 8, 8, 8, 6, 6 };

