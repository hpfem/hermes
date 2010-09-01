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

#ifndef _ADAPT_PROJIPOL_H_
#define _ADAPT_PROJIPOL_H_

#include "proj.h"

/// @defgroup hp-adapt hp-Adaptivity

/// Abstract class for projecting reference solution onto the coarse mesh using projection based
/// interpolation
///
/// NOTE: hex-specific
///
/// @ingroup hp-adapt
class ProjectionIpol : public Projection {
public:
	ProjectionIpol(Solution *afn, Element *e, Shapeset *ss);
	virtual ~ProjectionIpol();

	virtual double get_error(int split, int son, const order3_t &order) = 0;

protected:
	struct ProjItem {
		scalar coef;				// coef of the projection
		int idx;					// index of the shape function
	};

	ProjItem *vertex_proj;			// vertex projection
	ProjItem *edge_proj[12];		// edge projection
	ProjItem *face_proj[6];			// face projection
	ProjItem *bubble_proj;			// bubble projection

	ProjItem **proj;				// projection
	int proj_fns;					// number of funcions

	void free_proj();
	void calc_projection(int split, int son, const order3_t &order);
	virtual void calc_vertex_proj(int split, int son) = 0;
	virtual void calc_edge_proj(int edge, int split, int son, const order3_t &order) = 0;
	virtual void calc_face_proj(int face, int split, int son, const order3_t &order) = 0;
	virtual void calc_bubble_proj(int split, int son, const order3_t &order) = 0;
};

#endif
