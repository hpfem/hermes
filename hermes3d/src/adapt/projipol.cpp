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

#include "../h3dconfig.h"
#include "../function.h"
#include "../solution.h"
#include "projipol.h"
#include <common/callstack.h>

ProjectionIpol::ProjectionIpol(Solution *afn, Element *e, Shapeset *ss) : Projection(afn, e, ss)
{
	_F_
	// null
	vertex_proj = NULL;
	for (int i = 0; i < Hex::NUM_EDGES; i++) edge_proj[i] = NULL;
	for (int i = 0; i < Hex::NUM_FACES; i++) face_proj[i] = NULL;
	bubble_proj = NULL;
	proj = NULL;
	proj_fns = 0;
}

ProjectionIpol::~ProjectionIpol()
{
	_F_
	delete fu;
	delete fv;

	free_proj();
}

void ProjectionIpol::free_proj()
{
	_F_
	delete [] vertex_proj;
	for (int i = 0; i < Hex::NUM_EDGES; i++) delete [] edge_proj[i];
	for (int i = 0; i < Hex::NUM_FACES; i++) delete [] face_proj[i];
	delete [] bubble_proj;

	delete [] proj;
}

void ProjectionIpol::calc_projection(int split, int son, const order3_t &order)
{
	_F_
	free_proj();
	calc_vertex_proj(split, son + 1);
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
		calc_edge_proj(iedge, split, son + 1, order);
	for (int iface = 0; iface < Hex::NUM_FACES; iface++)
		calc_face_proj(iface, split, son + 1, order);
	calc_bubble_proj(split, son + 1, order);

	proj_fns = (order.x + 1) * (order.y + 1) * (order.z + 1);
	proj = new ProjItem *[proj_fns];

	int m = 0;
	for (int i = 0; i < Hex::NUM_VERTICES; i++, m++)
		proj[m] = vertex_proj + i;
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		order1_t edge_order = order.get_edge_order(iedge);
		int edge_fns = edge_order - 1;
		for (int i = 0; i < edge_fns; i++, m++)
			proj[m] = edge_proj[iedge] + i;
	}
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		order2_t face_order = order.get_face_order(iface);
		int face_fns = (face_order.x - 1) * (face_order.y - 1);
		for (int i = 0; i < face_fns; i++, m++)
			proj[m] = face_proj[iface] + i;
	}
	int bubble_fns = (order.x - 1) * (order.y - 1) * (order.z - 1);
	for (int i = 0; i < bubble_fns; i++, m++)
		proj[m] = bubble_proj + i;
}
