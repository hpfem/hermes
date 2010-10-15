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

#ifndef _REFDOMAIN_H_
#define _REFDOMAIN_H_

#include "mesh.h"
#include "common.h"

/// @defgroup ref_domains Reference domains
///
/// Reference domains are used... ;)
///


//
// 1D domains
//

/// Reference domain for a line (1D)
///
/// FIXME: name
///
/// @ingroup ref_domains
extern HERMES_API Point1D RefLine_vertices[];
	
extern HERMES_API Point2D RefTri_vertices[];
extern HERMES_API int2 RefTri_edge_vtcs[];

extern HERMES_API Point2D RefQuad_vertices[];
extern HERMES_API int2 RefQuad_edge_vtcs[];


extern HERMES_API Point3D RefTetra_vertices[];
extern HERMES_API int2 RefTetra_edge_vtcs[];
extern HERMES_API int *RefTetra_face_vtcs[];
extern HERMES_API int *RefTetra_face_edges[];
extern HERMES_API int RefTetra_face_nvtcs[];
extern HERMES_API int RefTetra_face_nedges[];
extern HERMES_API EMode2D RefTetra_face_mode[];
extern HERMES_API int RefTetra_face_orientations[];
extern HERMES_API Point3D RefTetra_face_normal[];


extern HERMES_API Point3D RefHex_vertices[];
extern HERMES_API int2 RefHex_edge_vtcs[];
extern HERMES_API int *RefHex_face_vtcs[];
extern HERMES_API int *RefHex_face_edges[];
extern HERMES_API int RefHex_face_nvtcs[];
extern HERMES_API int RefHex_face_nedges[];
extern HERMES_API EMode2D RefHex_face_mode[];
extern HERMES_API int RefHex_face_orientations[];
extern HERMES_API int RefHex_face_edge_ori[8][2];
extern HERMES_API int RefHex_edge_tangent[];
extern HERMES_API int2 RefHex_face_tangent[];

extern HERMES_API Point3D RefPrism_vertices[];
extern HERMES_API int2 RefPrism_edge_vtcs[];
extern HERMES_API int *RefPrism_face_vtcs[];
extern HERMES_API int *RefPrism_face_edges[];
extern HERMES_API int RefPrism_face_nvtcs[];
extern HERMES_API int RefPrism_face_nedges[];
extern HERMES_API EMode2D RefPrism_face_mode[];
extern HERMES_API int RefPrism_face_orientations[];

class HERMES_API RefLine {
public:
	static const Point1D *get_vertices() { return RefLine_vertices; }
};


//
// 2D domains
//

/// Reference domain for triangle (2D)
///
/// @ingroup ref_domains
class HERMES_API RefTri {
public:
	static const Point2D *get_vertices() { return RefTri_vertices; }

	static const int *get_edge_vertices(int edge) { return RefTri_edge_vtcs[edge]; }
};


/// Reference domain for quadrilateral (2D)
///
/// @ingroup ref_domains
class HERMES_API RefQuad {
public:
	static const Point2D *get_vertices() { return RefQuad_vertices; }

	static const int *get_edge_vertices(int edge) { return RefQuad_edge_vtcs[edge]; }
};


//
// 3D domains
//

/// Reference domain for tetrahedron (3D)
///
/// @ingroup ref_domains
class HERMES_API RefTetra {
public:
	static const Point3D *get_vertices() { return RefTetra_vertices; }

	static const int *get_edge_vertices(int edge) { return RefTetra_edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return RefTetra_face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return RefTetra_face_nedges[face]; }
	static const int *get_face_vertices(int face) { return RefTetra_face_vtcs[face]; }
	static const int *get_face_edges(int face) { return RefTetra_face_edges[face]; }
	static EMode2D get_face_mode(int face) { return RefTetra_face_mode[face]; }
	static int get_face_orientations(int face) { return RefTetra_face_orientations[face]; }

	static const Point3D get_face_normal(int iface) { return RefTetra_face_normal[iface]; }
};


/// Reference domain for hexahedron (3D)
///
/// @ingroup ref_domains
class HERMES_API RefHex {
public:
	static const Point3D *get_vertices() { return RefHex_vertices; }

	static const int *get_edge_vertices(int edge) { return RefHex_edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return RefHex_face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return RefHex_face_nedges[face]; }
	static const int *get_face_vertices(int face) { return RefHex_face_vtcs[face]; }
	static const int *get_face_edges(int face) { return RefHex_face_edges[face]; }
	static EMode2D get_face_mode(int face) { return RefHex_face_mode[face]; }
	static int get_face_orientations(int face) { return RefHex_face_orientations[face]; }
	/// @param[in] ori - face orientation
	/// @return orientations of edges on a face
	/// 	all edges are oriented the same way on all faces, so we do not have to care about face number
	static const int *get_face_edge_orientation(int ori) { return RefHex_face_edge_ori[ori]; }

	static int get_edge_tangent_direction(int edge) { return RefHex_edge_tangent[edge]; }
	static int get_face_tangent_direction(int face, int which) { return RefHex_face_tangent[face][which]; }
};


/// Reference domain for prism (3D)
///
/// @ingroup ref_domains
class HERMES_API RefPrism {
public:
	static const Point3D *get_vertices() { return RefPrism_vertices; }

	static const int *get_edge_vertices(int edge) { return RefPrism_edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return RefPrism_face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return RefPrism_face_nedges[face]; }
	static const int *get_face_vertices(int face) { return RefPrism_face_vtcs[face]; }
	static const int *get_face_edges(int face) { return RefPrism_face_edges[face]; }
	static EMode2D get_face_mode(int face) { return RefPrism_face_mode[face]; }
	static int get_face_orientations(int face) { return RefPrism_face_orientations[face]; }
};

#endif
