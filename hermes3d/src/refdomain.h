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
class RefLine {
public:
	static const Point1D *get_vertices() { return vertices; }

protected:
	static const Point1D vertices[];
};


//
// 2D domains
//

/// Reference domain for triangle (2D)
///
/// @ingroup ref_domains
class RefTri {
public:
	static const Point2D *get_vertices() { return vertices; }

	static const int *get_edge_vertices(int edge) { return edge_vtcs[edge]; }

protected:
	static const Point2D vertices[];
	static const int2 edge_vtcs[];
};


/// Reference domain for quadrilateral (2D)
///
/// @ingroup ref_domains
class RefQuad {
public:
	static const Point2D *get_vertices() { return vertices; }

	static const int *get_edge_vertices(int edge) { return edge_vtcs[edge]; }

protected:
	static const Point2D vertices[];
	static const int2 edge_vtcs[];
};


//
// 3D domains
//

/// Reference domain for tetrahedron (3D)
///
/// @ingroup ref_domains
class RefTetra {
public:
	static const Point3D *get_vertices() { return vertices; }

	static const int *get_edge_vertices(int edge) { return edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return face_nedges[face]; }
	static const int *get_face_vertices(int face) { return face_vtcs[face]; }
	static const int *get_face_edges(int face) { return face_edges[face]; }
	static EMode2D get_face_mode(int face) { return face_mode[face]; }
	static int get_face_orientations(int face) { return face_orientations[face]; }

	static const Point3D get_face_normal(int iface) { return face_normal[iface]; }

protected:
	static const Point3D vertices[];
	static const int2 edge_vtcs[];
	static const int *face_vtcs[];
	static const int *face_edges[];
	static const int face_nvtcs[];
	static const int face_nedges[];
	static const EMode2D face_mode[];
	static const int face_orientations[];

	static const Point3D face_normal[];
};


/// Reference domain for hexahedron (3D)
///
/// @ingroup ref_domains
class RefHex {
public:
	static const Point3D *get_vertices() { return vertices; }

	static const int *get_edge_vertices(int edge) { return edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return face_nedges[face]; }
	static const int *get_face_vertices(int face) { return face_vtcs[face]; }
	static const int *get_face_edges(int face) { return face_edges[face]; }
	static EMode2D get_face_mode(int face) { return face_mode[face]; }
	static int get_face_orientations(int face) { return face_orientations[face]; }
	/// @param[in] ori - face orientation
	/// @return orientations of edges on a face
	/// 	all edges are oriented the same way on all faces, so we do not have to care about face number
	static const int *get_face_edge_orientation(int ori) { return face_edge_ori[ori]; }

	static int get_edge_tangent_direction(int edge) { return edge_tangent[edge]; }
	static int get_face_tangent_direction(int face, int which) { return face_tangent[face][which]; }

protected:
	static const Point3D vertices[];
	static const int2 edge_vtcs[];
	static const int *face_vtcs[];
	static const int *face_edges[];
	static const int face_nvtcs[];
	static const int face_nedges[];
	static const EMode2D face_mode[];
	static const int face_orientations[];
	static const int face_edge_ori[8][2];
	static const int edge_tangent[];
	static const int2 face_tangent[];
};


/// Reference domain for prism (3D)
///
/// @ingroup ref_domains
class RefPrism {
public:
	static const Point3D *get_vertices() { return vertices; }

	static const int *get_edge_vertices(int edge) { return edge_vtcs[edge]; }
	static int get_edge_orientations() { return 2; }		// two orientations of an edge

	static int get_num_face_vertices(int face) { return face_nvtcs[face]; }
	static int get_num_face_edges(int face) { return face_nedges[face]; }
	static const int *get_face_vertices(int face) { return face_vtcs[face]; }
	static const int *get_face_edges(int face) { return face_edges[face]; }
	static EMode2D get_face_mode(int face) { return face_mode[face]; }
	static int get_face_orientations(int face) { return face_orientations[face]; }

protected:
	static const Point3D vertices[];
	static const int2 edge_vtcs[];
	static const int *face_vtcs[];
	static const int *face_edges[];
	static const int face_nvtcs[];
	static const int face_nedges[];
	static const EMode2D face_mode[];
	static const int face_orientations[];
};

#endif
