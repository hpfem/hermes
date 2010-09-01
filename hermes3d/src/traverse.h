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

#ifndef _TRAVERSE_H_
#define _TRAVERSE_H_

struct FacePos {
	int marker;			///< face marker
	int face;			///< element face number (local)

	// for internal use
	Element *base;
	Space *space, *space_u, *space_v;
};

struct Mesh;
class  Transformable;
struct State;
struct Box;


struct UniData {
	Element *e;
	uint64 idx;
};


/// Traverse is a multi-mesh traversal utility class. Given N meshes sharing the
/// same base mesh it walks through all (pseudo-)elements of the union of all
/// the N meshes.
///
class Traverse {
public:
	void begin(int n, Mesh **meshes, Transformable **fn = NULL);
	void finish();

	Element **get_next_state(bool *bnd, FacePos *ep);
	Element  *get_base() const { return base; }

	UniData **construct_union_mesh(Mesh *unimesh);

private:
	int num;
	Mesh **meshes;
	Transformable **fn;

	State *stack;
	int top, size;

	Word_t id;
	Element *base;
	int (*sons)[8];
	uint64 *subs;

	UniData **unidata;
	Word_t udsize;

	State *push_state();
	void set_boundary_info(State *s, bool *bnd, FacePos *ep);
	void union_recurrent(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni);

	void hex_union_rec(Box *cr, Element **e, Box *er, uint64 *idx, Element *uni);
	void hex_push_son_states(State *s);

	Mesh *unimesh;
};


#endif
