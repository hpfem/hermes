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

#include "../../hermes_common/trace.h"
#include "../../hermes_common/error.h"

#include "mesh.h"
#include "meshloader.h"
#include "refdomain.h"
#include "refmap.h"
#include "determinant.h"

const int TOP_LEVEL_REF = -1;

// forward declarations
extern bool verbose;

// Vertex /////////////////////////////////////////////////////////////////////

Vertex::Vertex() {
	x = y = z = 0.0;
}

Vertex::Vertex(double _x, double _y, double _z) {
	x = _x;
	y = _y;
	z = _z;
}

Vertex::Vertex(const Vertex &o) {
	x = o.x;
	y = o.y;
	z = o.z;
}

Vertex::~Vertex() {
}

Vertex *Vertex::copy() {
	return new Vertex(*this);
}

void Vertex::dump() {
	printf("(x = %lf, y = %lf, z = %lf)\n", x, y, z);
}

Edge::Key Edge::invalid_key = Edge::Key(NULL, 0);

Edge::Edge() {
	bnd = 0;
	ref = 0;
}

// for debugging
void Edge::dump() {
	printf("bnd = %d, ref = %d\n", bnd, ref);
}

// Facet //////////////////////////////////////////////////////////////////////

Facet::Key Facet::invalid_key = Facet::Key(NULL, 0);

Facet::Facet(ElementMode2D mode) {
	_F_
	this->mode = mode;
	this->type = Facet::INNER;
	this->left = INVALID_IDX;
	this->right = INVALID_IDX;
	this->left_face_num = -1;
	this->right_face_num = -1;
	this->lactive = false;
	this->ractive = false;
	this->ref_mask = H3D_REFT_FACE_NONE;

  this->parent = this->invalid_key;
	for (int i = 0; i < MAX_SONS; i++)
		this->sons[i] = this->invalid_key;
}

Facet::Facet(const Facet &o) {
	_F_
	mode = o.mode;
	lactive = o.lactive;
	ractive = o.ractive;
	type = o.type;
	left = o.left;
	right = o.right;
	left_face_num = o.left_face_num;
	right_face_num = o.right_face_num;
	ref_mask = o.ref_mask;

	parent = o.parent;
	for (int i = 0; i < MAX_SONS; i++)
		sons[i] = o.sons[i];
}

Facet::~Facet() {
	_F_
}

Facet *Facet::copy() {
	_F_
	return new Facet(*this);
}

Facet *Facet::copy_base() {
	_F_
	Facet *copy = new Facet(mode);
	MEM_CHECK(copy);

	copy->lactive = true;
	copy->ractive = true;
	copy->type = type;
	copy->left = left;
	copy->right = right;
	copy->left_face_num = left_face_num;
	copy->right_face_num = right_face_num;

	return copy;
}

bool Facet::ced(unsigned int idx, int iface)
{
	return type == Facet::INNER &&
		((lactive && !ractive) ||
		 (!lactive && ractive));
}

// for debugging
void Facet::dump() {
	_F_
	const char *s_type[] = { "INNER", "OUTER" };
	const char *s_mode[] = { "TRI", "QUAD" };

	printf("type = %s (%s), [%d, %d], left (elem = %d, face = %d), ", s_type[type], s_mode[mode], lactive, ractive, left, left_face_num);
	if (type == INNER) printf(" right (elem = %d, face = %d)", right, right_face_num);
	else printf(" right (bdr = %u)", right);
	if (parent != invalid_key) printf("parent");
	else printf("no parent");
	printf("\n");
}

// Element ////////////////////////////////////////////////////////////////////

Element::Element() {
	_F_
	id = INVALID_IDX;
	marker = -1;
	iro_cache = -1;

	active = true;
	used = 1;
	reft = 0;
}

Element::Element(const Element &o) {
	_F_
	id = o.id;
	marker = o.marker;

	iro_cache = o.iro_cache;
	active = o.active;
	used = o.used;
	reft = o.reft;
}

Element::~Element() {
	_F_
}

void Element::dump() {
	_F_
	printf("id = %u\n", id);
}

// Hex ////////////////////////////////////////////////////////////////////////

Hex::Hex() {
	_F_
#ifdef WITH_HEX
	for (int i = 0; i < NUM_SONS; i++)
		sons[i] = INVALID_IDX;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

Hex::Hex(unsigned int v[]) {
	_F_
#ifdef WITH_HEX
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];

	for (int i = 0; i < NUM_SONS; i++)
		sons[i] = INVALID_IDX;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

Hex::Hex(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, unsigned int v5, unsigned int v6, unsigned int v7, unsigned int v8) {
	_F_
#ifdef WITH_HEX
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
	vtcs[4] = v5;
	vtcs[5] = v6;
	vtcs[6] = v7;
	vtcs[7] = v8;

	for (int i = 0; i < NUM_SONS; i++)
		sons[i] = INVALID_IDX;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

Hex::Hex(const Hex &o) : Element(o) {
	_F_
#ifdef WITH_HEX
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = o.vtcs[i];
	for (int i = 0; i < NUM_SONS; i++)
		sons[i] = o.sons[i];
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

Hex::~Hex() {
	_F_
}

int Hex::get_edge_vertices(int edge_num, unsigned int *vtcs) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	const int *local_vetrex = RefHex::get_edge_vertices(edge_num);
	vtcs[0] = this->vtcs[local_vetrex[0]];
	vtcs[1] = this->vtcs[local_vetrex[1]];
	return Edge::NUM_VERTICES;
}

const int *Hex::get_edge_vertices(int edge_num) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	return RefHex::get_edge_vertices(edge_num);
}

int Hex::get_face_vertices(int face_num, unsigned int *vtcs) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	const int *local_numbers = RefHex::get_face_vertices(face_num);
	for (int i = 0; i < Quad::NUM_VERTICES; i++)
		vtcs[i] = this->vtcs[local_numbers[i]];
	return Quad::NUM_VERTICES;
}

const int *Hex::get_face_vertices(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	return RefHex::get_face_vertices(face_num);
}

const int *Hex::get_face_edges(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	return RefHex::get_face_edges(face_num);
}

int Hex::get_edge_orientation(int edge_num) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	const int *vert_idx = RefHex::get_edge_vertices(edge_num);
	// 0 - the edge orientation on the physical domain agrees with the orientation on reference domain
	// 1 - egde orientation is opposite
	if (vtcs[vert_idx[0]] < vtcs[vert_idx[1]])
		return 0;
	else
		return 1;
}

int Hex::get_face_orientation(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	unsigned int v[4];
	get_face_vertices(face_num, v);

	unsigned int minval = 1000;
	int min = 0;
	for (int i = 0; i < 4; i++) {
		if (v[i] < minval) {
			minval = v[i];
			min = i;
		}
	}

	assert(min >= 0 && min <= 3); // Error in face orientation
	// FIXME: replace magic numbers
	if (min == 0) return (v[1] < v[3]) ? 0 : 4;
	else if (min == 1) return (v[0] < v[2]) ? 1 : 6;
	else if (min == 2) return (v[3] < v[1]) ? 3 : 7;
	else if (min == 3) return (v[2] < v[0]) ? 2 : 5;
	else return -1;
}

Element *Hex::copy() {
	_F_
	return new Hex(*this);
}

Element *Hex::copy_base() {
	_F_
#ifdef WITH_HEX
	Hex *copy = new Hex(vtcs);
	MEM_CHECK(copy);
	copy->id = id;

	return copy;
#else
	return NULL;
#endif
}

void Hex::ref_all_nodes() {
}

void Hex::unref_all_nodes() {
}

// for debugging
void Hex::dump() {
	_F_
	printf("id = %u (%u, %u, %d), vertices(%u, %u, %u, %u, %u, %u, %u, %u), ", id, active, used, reft,
		vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
	printf("sons(%d, %d, %d, %d, %d, %d, %d, %d), ",
		sons[0], sons[1], sons[2], sons[3], sons[4], sons[5], sons[6], sons[7]);
	printf("marker = %d\n", marker);
}

// Tetra //////////////////////////////////////////////////////////////////////

Tetra::Tetra() {
	_F_
#ifdef WITH_TETRA
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

Tetra::Tetra(unsigned int v[]) {
	_F_
#ifdef WITH_TETRA
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

Tetra::Tetra(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4) {
	_F_
#ifdef WITH_TETRA
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

Tetra::Tetra(const Tetra &o) :
	Element(o)
{
	_F_
#ifdef WITH_TETRA
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = o.vtcs[i];
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

Tetra::~Tetra() {
	_F_
}

int Tetra::get_edge_vertices(int edge_num, unsigned int *vtcs) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	vtcs[0] = this->vtcs[RefTetra::get_edge_vertices(edge_num)[0]];
	vtcs[1] = this->vtcs[RefTetra::get_edge_vertices(edge_num)[1]];
	return Edge::NUM_VERTICES;
}

const int *Tetra::get_edge_vertices(int edge_num) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	return RefTetra::get_edge_vertices(edge_num);
}

int Tetra::get_face_vertices(int face_num, unsigned int *vtcs) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	const int *local_numbers = RefTetra::get_face_vertices(face_num);
	for (int i = 0; i < Tri::NUM_VERTICES; i++)
		vtcs[i] = this->vtcs[local_numbers[i]];
	return Tri::NUM_VERTICES;
}

const int *Tetra::get_face_vertices(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	return RefTetra::get_face_vertices(face_num);
}

const int *Tetra::get_face_edges(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	return RefTetra::get_face_edges(face_num);
}

int Tetra::get_edge_orientation(int edge_num) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	static int map[6][2] = { { 0, 1 }, { 1, 2 }, { 0, 2 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };

	// 0 - the edge orientation on the physical domain agrees with the orientation on reference domain
	// 1 - egde orientation is opposite
	if (vtcs[map[edge_num][0]] < vtcs[map[edge_num][1]])
		return 0;
	else
		return 1;
}

int Tetra::get_face_orientation(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	static int map[4][3] = { { 0, 1, 3 }, { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 2 } };

	unsigned int v0 = vtcs[map[face_num][0]];
	unsigned int v1 = vtcs[map[face_num][1]];
	unsigned int v2 = vtcs[map[face_num][2]];

	if (v0 < v1 && v1 < v2) return 0;
	else if (v1 < v2 && v2 < v0) return 1;
	else if (v2 < v0 && v0 < v1) return 2;
	else if (v0 < v2 && v2 < v1) return 3;
	else if (v1 < v0 && v0 < v2) return 4;
	else if (v2 < v1 && v1 < v0) return 5;
	else return -1;
}

Element *Tetra::copy() {
	_F_
	return new Tetra(*this);
}

Element *Tetra::copy_base() {
	_F_
#ifdef WITH_TETRA
	Tetra *copy = new Tetra();
	for (int i = 0; i < NUM_VERTICES; i++)
		copy->vtcs[i] = vtcs[i];

	return copy;
#else
	return NULL;
#endif
}

void Tetra::ref_all_nodes() {
}

void Tetra::unref_all_nodes() {
}

// for debugging
void Tetra::dump() {
	printf("id = %u, vertices(%u, %u, %u, %u), ", id, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
	printf("marker = %d\n", marker);
}

// Prism //////////////////////////////////////////////////////////////////////

Prism::Prism() {
	_F_
#ifdef WITH_PRISM
#else
	EXIT(H3D_ERR_PRISM_NOT_COMPILED);
#endif
}

Prism::Prism(unsigned int v[]) {
	_F_
#ifdef WITH_PRISM
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
#else
	EXIT(H3D_ERR_PRISM_NOT_COMPILED);
#endif
}

Prism::Prism(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, unsigned int v5, unsigned int v6) {
	_F_
#ifdef WITH_PRISM
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
	vtcs[4] = v5;
	vtcs[5] = v6;
#else
	EXIT(H3D_ERR_PRISM_NOT_COMPILED);
#endif
}

Prism::Prism(const Prism &o) :
	Element(o)
{
	_F_
#ifdef WITH_PRISM
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = o.vtcs[i];
#else
	EXIT(H3D_ERR_PRISM_NOT_COMPILED);
#endif
}

Prism::~Prism() {
}

int Prism::get_edge_vertices(int edge_num, unsigned int *vtcs) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	vtcs[0] = this->vtcs[RefPrism::get_edge_vertices(edge_num)[0]];
	vtcs[1] = this->vtcs[RefPrism::get_edge_vertices(edge_num)[1]];
	return Edge::NUM_VERTICES;
}

const int *Prism::get_edge_vertices(int edge_num) const {
	_F_
	assert((edge_num >= 0) && (edge_num < NUM_EDGES));
	return RefPrism::get_edge_vertices(edge_num);
}

ElementMode2D Prism::get_face_mode(int face_num) const {
	_F_
	assert((face_num >= 0) && (face_num < NUM_FACES));
	return RefPrism::get_face_mode(face_num);
}

int Prism::get_num_face_vertices(int face_num) const {
	_F_
	assert((face_num >= 0) && (face_num < NUM_FACES));
	return RefPrism::get_num_face_vertices(face_num);
}

int Prism::get_num_face_edges(int face_num) const {
	_F_
	assert((face_num >= 0) && (face_num < NUM_FACES));
	return RefPrism::get_num_face_edges(face_num);
}

int Prism::get_face_vertices(int face_num, unsigned int *vtcs) const {
	_F_
	assert((face_num >= 0) && (face_num < NUM_FACES));
	int nvert = RefPrism::get_num_face_vertices(face_num);
	const int *local_numbers = RefPrism::get_face_vertices(face_num);
	for (int i = 0; i < nvert; i++)
		vtcs[i] = this->vtcs[local_numbers[i]];
	return nvert;
}

const int *Prism::get_face_vertices(int face_num) const {
	_F_
	assert((face_num >= 0) && (face_num < NUM_FACES));
	return RefPrism::get_face_vertices(face_num);
}

const int *Prism::get_face_edges(int face_num) const {
	_F_
	assert(face_num >= 0 && face_num < NUM_FACES);
	return RefPrism::get_face_edges(face_num);
}

int Prism::get_edge_orientation(int edge_num) const {
	_F_
	EXIT(HERMES_ERR_NOT_IMPLEMENTED);
	// FIXME
	return -1;
}

int Prism::get_face_orientation(int face_num) const {
	_F_
	EXIT(HERMES_ERR_NOT_IMPLEMENTED);
	// FIXME
	return -1;
}

Element *Prism::copy() {
	_F_
	return new Prism(*this);
}

Element *Prism::copy_base() {
	_F_
#ifdef WITH_PRISM
	Prism *copy = new Prism();
	for (int i = 0; i < NUM_VERTICES; i++)
		copy->vtcs[i] = vtcs[i];

	return copy;
#else
	return NULL;
#endif
}

void Prism::ref_all_nodes() {
}

void Prism::unref_all_nodes() {
}

// for debugging
void Prism::dump() {
	printf("id = %u, vertices(%u, %u, %u, %u, %u, %u), ", id,
		vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5]);
	printf("marker = %d\n", marker);
}

// Boundary ///////////////////////////////////////////////////////////////////

Boundary::Boundary(int marker) {
	this->marker = marker;
	this->id = -1;
}

Boundary::Boundary(const Boundary &o) {
	marker = o.marker;
	id = o.id;
}

Boundary::~Boundary() {
}

// for debugging
void Boundary::dump() {
	printf("id = %u, marker = %d\n", id, marker);
}

// BoundaryTri ////////////////////////////////////////////////////////////////

BoundaryTri::BoundaryTri(int marker) :
	Boundary(marker) {
}

BoundaryTri::BoundaryTri(const BoundaryTri &o) :
	Boundary(o) {
}

BoundaryTri::~BoundaryTri() {
}

Boundary *BoundaryTri::copy() {
	_F_
	return new BoundaryTri(*this);
}

// BoundaryQuad ///////////////////////////////////////////////////////////////

BoundaryQuad::BoundaryQuad(int marker) :
	Boundary(marker)
{
}

BoundaryQuad::BoundaryQuad(const BoundaryQuad &o) :
	Boundary(o)
{
}

BoundaryQuad::~BoundaryQuad() {
}

Boundary *BoundaryQuad::copy() {
	_F_
	return new BoundaryQuad(*this);
}

// Mesh ///////////////////////////////////////////////////////////////////////

int g_mesh_seq = 0;

Mesh::Mesh() {
	_F_

	nactive = 0;
	nbase = 0;
	seq = g_mesh_seq++;
}

Mesh::~Mesh() {
	_F_
	free();
}

void Mesh::free() {
	_F_
	for(std::map<unsigned int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++)
    delete it->second;
  vertices.clear();

	for(std::map<unsigned int, Element*>::iterator it = elements.begin(); it != elements.end(); it++)
    delete it->second;
  elements.clear();

	for(std::map<unsigned int, Boundary*>::iterator it = boundaries.begin(); it != boundaries.end(); it++)
    delete it->second;
  boundaries.clear();

  for(std::map<Facet::Key, Facet*>::iterator it = facets.begin(); it != facets.end(); it++)
    delete it->second;
  facets.clear();

  for(std::map<Edge::Key, Edge*>::iterator it = edges.begin(); it != edges.end(); it++)
    delete it->second;
  edges.clear();

	midpoints.clear();
}

void Mesh::copy(const Mesh &mesh) 
{
  _F_
  if (&mesh == this) warning("Copying mesh into itself.");
  free();

  for(std::map<unsigned int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++) delete it->second;
  vertices.clear();

  for(std::map<unsigned int, Element*>::iterator it = elements.begin(); it != elements.end(); it++) delete it->second;
  elements.clear();

  for(std::map<unsigned int, Boundary*>::iterator it = boundaries.begin(); it != boundaries.end(); it++) delete it->second;
  boundaries.clear();

  for(std::map<Edge::Key, Edge*>::iterator it = edges.begin(); it != edges.end(); it++) delete it->second;
  edges.clear();

  for(std::map<Facet::Key, Facet*>::iterator it = facets.begin(); it != facets.end(); it++) delete it->second;
  facets.clear();

  midpoints.clear();

  // copy vertices
  for(std::map<unsigned int, Vertex*>::const_iterator it = mesh.vertices.begin(); it != mesh.vertices.end(); it++)
    if(it->first != INVALID_IDX)
      this->vertices[it->first] = it->second->copy();

  // copy boundaries
  for(std::map<unsigned int, Boundary*>::const_iterator it = mesh.boundaries.begin(); it != mesh.boundaries.end(); it++)
    if(it->first != INVALID_IDX)
      this->boundaries[it->first] = it->second->copy();

  // copy elements, midpoints, facets and edges
  for(std::map<unsigned int, Element*>::const_iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++) {
    if(it->first == INVALID_IDX)
      continue;
    Element *e = it->second;
    this->elements[it->first] = e->copy();

		// copy mid points on edges and edges
		unsigned int *emp = new unsigned int[e->get_num_edges()];
		for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
			unsigned int *edge_vtx = new unsigned int[Edge::NUM_VERTICES];
			e->get_edge_vertices(iedge, edge_vtx);
			emp[iedge] = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
			if (emp[iedge] != INVALID_IDX) {
        MidPointKey key(edge_vtx[0], edge_vtx[1]);
        if(midpoints.find(key) == midpoints.end())
          create_midpoint(edge_vtx[0], edge_vtx[1]);
				set_midpoint(edge_vtx[0], edge_vtx[1], emp[iedge]);
      }

      Edge::Key key(edge_vtx + 0, Edge::NUM_VERTICES);
      if (mesh.edges.at(key) != NULL) {
        edges[key] = new Edge;
				*edges[key] = *mesh.edges.at(key);
      }
      delete [] edge_vtx;
		}

		// copy mid points on faces
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			const int *edge = e->get_face_edges(iface);

			switch (e->get_mode()) {
				case HERMES_MODE_HEX:
					// horz
					if (emp[edge[1]] != INVALID_IDX && emp[edge[3]] != INVALID_IDX) {
						unsigned int edge_vtx[Edge::NUM_VERTICES] = { emp[edge[1]], emp[edge[3]] };

						unsigned int fmp = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
						if (fmp != INVALID_IDX)
							set_midpoint(edge_vtx[0], edge_vtx[1], fmp);
					}
					// vert
					if (emp[edge[0]] != INVALID_IDX && emp[edge[2]] != INVALID_IDX) {
						unsigned int edge_vtx[Edge::NUM_VERTICES] = { emp[edge[0]], emp[edge[2]] };

						unsigned int fmp = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
						if (fmp != INVALID_IDX)
							set_midpoint(edge_vtx[0], edge_vtx[1], fmp);
					}
					break;

				case HERMES_MODE_TET:
				case HERMES_MODE_PRISM:
					break;
			}
		}
    delete [] emp;
  }

  // facets
  for(std::map<Facet::Key, Facet*>::const_iterator it = mesh.facets.begin(); it != mesh.facets.end(); it++) {
    Facet *facet = it->second;

      unsigned int *face_idxs = new unsigned int[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices
      if ((unsigned) facet->left != INVALID_IDX) {
	Element *left_e = mesh.elements.at(facet->left);
	int nvtcs = left_e->get_face_vertices(facet->left_face_num, face_idxs);
        Facet::Key key(face_idxs + 0, nvtcs);
        this->facets[key] = facet->copy();
      }
      else if ((unsigned) facet->right != INVALID_IDX && facet->type == Facet::INNER) {
	Element *right_e = mesh.elements.at(facet->right);
	int nvtcs = right_e->get_face_vertices(facet->right_face_num, face_idxs);
        Facet::Key key(face_idxs + 0, nvtcs);
        this->facets[key] = facet->copy();
      }
      else EXIT("Internal error in Mesh::copy().");		// FIXME
    
      delete [] face_idxs;
  }

  nbase = mesh.nbase;
  nactive = mesh.nactive;
  seq = mesh.seq;
}

void Mesh::copy_base(const Mesh &mesh) {
	_F_
	// TODO: improve this, it is not very nice, especially facet copying

	free();
	// copy elements, facets and edges
	for(std::map<unsigned int, Element*>::const_iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++) {
    if(it->first > mesh.nbase)
      continue;
    Element *e = it->second;

		// vertices
		for (int iv = 0; iv < e->get_num_vertices(); iv++) {
			unsigned int vtx = e->get_vertex(iv);
			if (this->vertices[vtx] == NULL)
				this->vertices[vtx] = mesh.vertices.at(vtx)->copy();
		}

		// edges
		for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
			unsigned int vtx[Edge::NUM_VERTICES];
			e->get_edge_vertices(iedge, vtx);

      Edge::Key key(vtx + 0, Edge::NUM_VERTICES);
      if (mesh.edges.at(key) != NULL) {
        edges[key] = new Edge;
				*edges[key] = *mesh.edges.at(key);
      }

		}

		// facets
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			unsigned int face_idxs[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices

			int nvts = e->get_face_vertices(iface, face_idxs);
			Facet *facet = NULL;
      Facet::Key key(face_idxs + 0, nvts);
      if(mesh.facets.find(key) != mesh.facets.end()) {
        facet = mesh.facets.at(key);
        if(this->facets.find(key) == this->facets.end()) {
					// insert the facet
					Facet *fcopy = facet->copy_base();
					fcopy->left = it->first;
					this->facets[key] = fcopy;

					// boundaries
					if (facet->type == Facet::OUTER) {
						unsigned int bnd_id = facet->right;
						if (this->boundaries[bnd_id] == NULL) 
              this->boundaries[bnd_id] = mesh.boundaries.at(bnd_id)->copy();
					}
				}
				else {
					assert(facet->type == Facet::INNER);
          this->facets.find(key)->second->set_right_info(it->first, iface);
				}
			}
		}

		this->elements[it->first] = e->copy_base();
	}

	this->nbase = this->nactive = mesh.nbase;
	this->seq = g_mesh_seq++;
}

Facet::Key Mesh::get_facet_id(Element *e, int face_num) const {
	_F_
	assert(e != NULL);
	unsigned int facet_idxs[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices
	int nvts = e->get_face_vertices(face_num, facet_idxs);
  return Facet::Key(facet_idxs + 0, nvts);
}

Edge::Key Mesh::get_edge_id(Element *e, int edge_num) const {
	_F_
	assert(e != NULL);
	unsigned int edge_idxs[Edge::NUM_VERTICES];
	int nvtcs = e->get_edge_vertices(edge_num, edge_idxs);
  return Edge::Key(edge_idxs + 0, nvtcs);
}

void Mesh::dump() {
    _F_
    printf("Vertices (count = %lu)\n", (unsigned long int)vertices.size());
    for(std::map<unsigned int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++) {
		Vertex *v = it->second;
    printf("  id = %d, ", it->first);
		v->dump();
	}

	printf("Elements (count = %lu)\n", (unsigned long int)elements.size());
  for(std::map<unsigned int, Element*>::iterator it = elements.begin(); it != elements.end(); it++) {
		Element *e = it->second;
		printf("  ");
		e->dump();
	}

	printf("Boundaries (count = %lu)\n", (unsigned long int)boundaries.size());
  for(std::map<unsigned int, Boundary*>::iterator it = boundaries.begin(); it != boundaries.end(); it++) {
		Boundary *b = it->second;
		printf("  ");
		b->dump();
	}

	printf("Facets (count = %lu)\n", (unsigned long int)facets.size());
  for(std::map<Facet::Key, Facet*>::iterator it = facets.begin(); it != facets.end(); it++) {
    Facet *f = it->second;
    if(it->first.size > 0)
      printf("Vertices: \n");
    for(unsigned int i = 0; i < it->first.size; i++)
      printf("%u: %u\t", i, it->first.vtcs[i]);
		f->dump();
	}
}

unsigned int Mesh::add_vertex(double x, double y, double z) {
	_F_
  vertices[vertices.size() + 1] = new Vertex(x, y, z);
  return vertices.size();
}

Tetra *Mesh::create_tetra(unsigned int vtcs[]) {
	_F_
	Tetra *tetra = new Tetra(vtcs);
	MEM_CHECK(tetra);
	
  unsigned int i;
  for(i = 1; ; i++)
    if(elements[i] == NULL)
      break;
  elements[i] = tetra;

	tetra->id = i;

	tetra->ref_all_nodes();

	return tetra;
}

Tetra *Mesh::add_tetra(unsigned int vtcs[]) {
	_F_
	Tetra *tetra = create_tetra(vtcs);

	// edges
	ref_edges(tetra);

	// facets
	for (int i = 0; i < Tetra::NUM_FACES; i++) {
		unsigned int facet_idxs[Tri::NUM_VERTICES];
		int nvtcs = tetra->get_face_vertices(i, facet_idxs);
    Facet::Key key(facet_idxs + 0, nvtcs);
    if (facets.find(key) != facets.end()) {
			facets[key]->type = Facet::INNER;
			facets[key]->set_right_info(tetra->id, i);
		}
		else {
			Facet* facet = new Facet(RefTetra::get_face_mode(i));
			facet->set_left_info(tetra->id, i);
			facets[key] = facet;
		}
	}

	return tetra;
}

Hex *Mesh::create_hex(unsigned int vtcs[]) {
	_F_
	// build up the element
	Hex *hex = new Hex(vtcs);
	MEM_CHECK(hex);
	
  unsigned int i;
  for(i = 1; ; i++)
    if(elements[i] == NULL)
      break;
  elements[i] = hex;

	hex->id = i;

	hex->ref_all_nodes();

	return hex;
}

Hex *Mesh::add_hex(unsigned int vtcs[]) 
{
  _F_
  Hex *hex = create_hex(vtcs);

  // edges
  ref_edges(hex);

  // facets
  for (int i = 0; i < Hex::NUM_FACES; i++) {
    unsigned int facet_idxs[Quad::NUM_VERTICES];
    int nvtcs = hex->get_face_vertices(i, facet_idxs);
    Facet::Key key(facet_idxs + 0, nvtcs);
    if (facets.find(key) != facets.end()) {
      facets[key]->type = Facet::INNER;
      facets[key]->set_right_info(hex->id, i);
    }
    else {
      Facet *fct = new Facet(HERMES_MODE_QUAD);
      MEM_CHECK(fct);
      fct->set_left_info(hex->id, i);
      facets[key] = fct;
    }
  }

  return hex;
}

Prism *Mesh::create_prism(unsigned int vtcs[])
{
  _F_
  Prism *prism = new Prism(vtcs);
  MEM_CHECK(prism);

  unsigned int i;
  for(i = 1; ; i++) if(elements[i] == NULL) break;
  elements[i] = prism;

  prism->id = i;

  prism->ref_all_nodes();

  return prism;
}

Prism *Mesh::add_prism(unsigned int vtcs[]) 
{
  _F_
  Prism *prism = create_prism(vtcs);

  // edges
  ref_edges(prism);

  // facets
  for (int i = 0; i < Prism::NUM_FACES; i++) {
    unsigned int facet_idxs[Quad::NUM_VERTICES];
    int nvtcs = prism->get_face_vertices(i, facet_idxs);
    Facet::Key key(facet_idxs + 0, nvtcs);
    if (facets.find(key) != facets.end()) {
      facets[key]->type = Facet::INNER;
      facets[key]->set_right_info(prism->id, i);
    }
    else {
      Facet *fct = new Facet(HERMES_MODE_QUAD);
      MEM_CHECK(fct);
      fct->set_left_info(prism->id, i);
      facets[key] = fct;
    }
  }
  return prism;
}

Boundary *Mesh::add_tri_boundary(unsigned int vtcs[], int marker) {
	_F_
  Facet::Key key(vtcs + 0, Tri::NUM_VERTICES);
  if (facets.find(key) != facets.end()) {
		Boundary *bdr = new BoundaryTri(marker);
		MEM_CHECK(bdr);

    unsigned int i;
    for(i = 1; ; i++)
      if(boundaries[i] == NULL)
        break;
    boundaries[i] = bdr;

		bdr->id = i;

		facets[key]->type = Facet::OUTER;
		facets[key]->set_right_info(bdr->id);

		return bdr;
	}
	else
		return NULL;
}

Boundary *Mesh::add_quad_boundary(unsigned int vtcs[], int marker) {
	_F_
	Facet::Key key(vtcs + 0, Quad::NUM_VERTICES);
  if (facets.find(key) != facets.end()) {
		Boundary *bdr = new BoundaryQuad(marker);
		MEM_CHECK(bdr);

    unsigned int i;
    for(i = 1; ; i++)
      if(boundaries[i] == NULL)
        break;
    boundaries[i] = bdr;

		bdr->id = i;

		facets[key]->type = Facet::OUTER;
		facets[key]->set_right_info(bdr->id);

		return bdr;
	}
	else
		return NULL;
}

void Mesh::ugh()
{
  _F_
  // set the number of active/base elements
  nactive = nbase = elements.size();

  // set bnd flag for boundary edges
  for(std::map<Facet::Key, Facet*>::iterator it = facets.begin(); it != facets.end(); it++) {
    Facet *facet = it->second;
    if (facet->type == Facet::OUTER) {
      Element *elem = elements[facet->left];
      const int *face_edge = elem->get_face_edges(facet->left_face_num);
      for (int iedge = 0; iedge < elem->get_num_face_edges(facet->left_face_num); iedge++) {
	unsigned int vtx[Edge::NUM_VERTICES];
	elem->get_edge_vertices(face_edge[iedge], vtx);

        Edge::Key key(vtx + 0, Edge::NUM_VERTICES);
        if(edges.find(key) == edges.end()) edges[key] = new Edge;
	edges[key]->bnd = 1;
      }
    }
  }

  check_elem_oris();
}

void Mesh::create_faces()
{
  // Go through all elements as if they were read from the file,
  // move the search for faces from the H3DReader's load() function to here.




}

bool Mesh::is_compatible_quad_refinement(Facet *facet, int reft) const 
{
	_F_
	if (facet->type == Facet::INNER) {
		// BOTH or NONE on the facet => refinement makes no problem
		if (facet->ref_mask == H3D_REFT_QUAD_BOTH || facet->ref_mask == H3D_REFT_FACE_NONE) return true;

		// applying BOTH or NONE is also no problem
		if (reft == H3D_REFT_QUAD_BOTH || reft == H3D_REFT_FACE_NONE) return true;

		unsigned int eid;
		int face_num;
		if (facet->ractive) {
			eid = facet->left;
			face_num = facet->left_face_num;
		}
		else if (facet->lactive) {
			eid = facet->right;
			face_num = facet->right_face_num;
		}
		else
			EXIT("Facet data corrupted or not a CED facet.");

		if (eid == INVALID_IDX) {
			// dunno what to do
			return false;
		}

		Element *e = elements.at(eid);
		int nv = e->get_num_face_vertices(face_num);
		unsigned int *face_vtx = new unsigned int[nv];
		e->get_face_vertices(face_num, face_vtx);

		// check if the vertices are there, if so => compatible refinement
		unsigned int emp[2] = { INVALID_IDX, INVALID_IDX };
		if (reft == H3D_REFT_QUAD_HORZ) {
			emp[0] = peek_midpoint(face_vtx[0], face_vtx[3]);
			emp[1] = peek_midpoint(face_vtx[1], face_vtx[2]);
		}
		else if (reft == H3D_REFT_QUAD_VERT) {
			emp[0] = peek_midpoint(face_vtx[0], face_vtx[1]);
			emp[1] = peek_midpoint(face_vtx[2], face_vtx[3]);
		}
  
    delete [] face_vtx;
		return (emp[0] != INVALID_IDX && emp[1] != INVALID_IDX);
	}
	else {
		// no problem with outer facets
		return true;
	}
}

bool Mesh::can_refine_element(unsigned int eid, int reft) const 
{
  _F_
  bool can_refine = false;

  Element *elem = elements.at(eid);
  assert(elem != NULL);
  switch (elem->get_mode()) {
    case HERMES_MODE_HEX: can_refine = can_refine_hex((Hex *) elem, reft); break;
    case HERMES_MODE_TET: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
    case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
    default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
  }

  return can_refine;
}

bool Mesh::can_refine_hex(Hex *elem, int refinement) const {
	_F_
	int nf; // number of faces to check
	int iface[Hex::NUM_FACES]; // face numbers to check
	int face_reft[Hex::NUM_FACES]; // refinements to apply on each face

	switch (refinement) {
		case H3D_REFT_HEX_NONE:
			return true;

		// refine to 2 hexes
		case H3D_REFT_HEX_X:
			nf = 4;
			iface[0] = 2; face_reft[0] = H3D_REFT_QUAD_VERT;
			iface[1] = 3; face_reft[1] = H3D_REFT_QUAD_VERT;
			iface[2] = 4; face_reft[2] = H3D_REFT_QUAD_VERT;
			iface[3] = 5; face_reft[3] = H3D_REFT_QUAD_VERT;
			break;

		case H3D_REFT_HEX_Y:
			nf = 4;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_VERT;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_VERT;
			iface[2] = 4; face_reft[2] = H3D_REFT_QUAD_HORZ;
			iface[3] = 5; face_reft[3] = H3D_REFT_QUAD_HORZ;
			break;

		case H3D_REFT_HEX_Z:
			nf = 4;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_HORZ;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_HORZ;
			iface[2] = 2; face_reft[2] = H3D_REFT_QUAD_HORZ;
			iface[3] = 3; face_reft[3] = H3D_REFT_QUAD_HORZ;
			break;

		// refine to 4 hexes
		case H3D_H3D_REFT_HEX_XY:
			nf = 6;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_VERT;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_VERT;
			iface[2] = 2; face_reft[2] = H3D_REFT_QUAD_VERT;
			iface[3] = 3; face_reft[3] = H3D_REFT_QUAD_VERT;
			iface[4] = 4; face_reft[4] = H3D_REFT_QUAD_BOTH;
			iface[5] = 5; face_reft[5] = H3D_REFT_QUAD_BOTH;
			break;

		case H3D_H3D_REFT_HEX_YZ:
			nf = 6;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_BOTH;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_BOTH;
			iface[2] = 2; face_reft[2] = H3D_REFT_QUAD_HORZ;
			iface[3] = 3; face_reft[3] = H3D_REFT_QUAD_HORZ;
			iface[4] = 4; face_reft[4] = H3D_REFT_QUAD_HORZ;
			iface[5] = 5; face_reft[5] = H3D_REFT_QUAD_HORZ;
			break;

		case H3D_H3D_REFT_HEX_XZ:
			nf = 6;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_HORZ;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_HORZ;
			iface[2] = 2; face_reft[2] = H3D_REFT_QUAD_BOTH;
			iface[3] = 3; face_reft[3] = H3D_REFT_QUAD_BOTH;
			iface[4] = 4; face_reft[4] = H3D_REFT_QUAD_VERT;
			iface[5] = 5; face_reft[5] = H3D_REFT_QUAD_VERT;
			break;

		// refine to 8 hexes
		case H3D_H3D_H3D_REFT_HEX_XYZ:
			nf = 6;
			iface[0] = 0; face_reft[0] = H3D_REFT_QUAD_BOTH;
			iface[1] = 1; face_reft[1] = H3D_REFT_QUAD_BOTH;
			iface[2] = 2; face_reft[2] = H3D_REFT_QUAD_BOTH;
			iface[3] = 3; face_reft[3] = H3D_REFT_QUAD_BOTH;
			iface[4] = 4; face_reft[4] = H3D_REFT_QUAD_BOTH;
			iface[5] = 5; face_reft[5] = H3D_REFT_QUAD_BOTH;
			break;
		default:
			EXIT(HERMES_ERR_UNKNOWN_REFINEMENT_TYPE);
			break;
	}

	bool can_refine = true;
	for (int i = 0; i < nf; i++) {
		Facet::Key fid = get_facet_id(elem, iface[i]);
    Facet *facet = facets.at(fid);
		assert(facet != NULL);
		can_refine &= is_compatible_quad_refinement(facet, face_reft[i]);
	}

	return can_refine;

}

bool Mesh::refine_element(unsigned int id, int refinement) {
	_F_
	bool refined = false;
	Element *elem = elements[id];
	assert(elem != NULL);
	if (can_refine_element(id, refinement)) {
		switch (elem->get_mode()) {
			case HERMES_MODE_HEX: refined = refine_hex((Hex *) elem, refinement); break;
			case HERMES_MODE_TET: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
			case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
			default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
		}

		seq = g_mesh_seq++;
	}
	else
		EXIT("Applying incompatible refinement (elem = %d, reft = %d).", id, refinement);

	return refined;
}

bool Mesh::refine_hex(Hex *elem, int refinement) {
	_F_
	assert(elem->active);		// Refinement already applied to element

	bool refined = false;
	switch (refinement) {
		case H3D_REFT_HEX_NONE: /* do nothing */
			break;
		// refine to 2 hexes
		case H3D_REFT_HEX_X:
		case H3D_REFT_HEX_Y:
		case H3D_REFT_HEX_Z:
			refined = refine_hex_2(elem, refinement);
			break;
		// refine to 4 hexes
		case H3D_H3D_REFT_HEX_XY:
		case H3D_H3D_REFT_HEX_XZ:
		case H3D_H3D_REFT_HEX_YZ:
			refined = refine_hex_4(elem, refinement);
			break;
		// refine to 8 hexes
		case H3D_H3D_H3D_REFT_HEX_XYZ:
			refined = refine_hex_8(elem, refinement);
			break;
		default:
			EXIT(HERMES_ERR_UNKNOWN_REFINEMENT_TYPE);
			break;
	}

	elem->reft = refinement;

	return refined;
}

// NOTE: can not use RefHex, since it has different face numbering on face 3
static const int hex_face_vtcs_0[] = { 0, 3, 7, 4 };
static const int hex_face_vtcs_1[] = { 1, 2, 6, 5 };
static const int hex_face_vtcs_2[] = { 0, 1, 5, 4 };
static const int hex_face_vtcs_3[] = { 3, 2, 6, 7 };
static const int hex_face_vtcs_4[] = { 0, 1, 2, 3 };
static const int hex_face_vtcs_5[] = { 4, 5, 6, 7 };

bool Mesh::refine_hex_2(Hex *parent, int refinement) {
	_F_
	bool refined = true;

	unsigned int vtx[Hex::NUM_VERTICES]; // vertices of parent element
	parent->get_vertices(vtx);

	unsigned int mp[4] = { 0 }; // four midpoints shared when refining to two elements
	const int *left, *right; // vertex indices for "left" and "right" face
	switch (refinement) {
		case H3D_REFT_HEX_X:
			left = hex_face_vtcs_0;
			right = hex_face_vtcs_1;
			break;

		case H3D_REFT_HEX_Y:
			left = hex_face_vtcs_2;
			right = hex_face_vtcs_3;
			break;

		case H3D_REFT_HEX_Z:
			left = hex_face_vtcs_4;
			right = hex_face_vtcs_5;
			break;
	}

	// create midpoint
	for (int i = 0; i < Quad::NUM_VERTICES; i++)
		mp[i] = get_midpoint(vtx[left[i]], vtx[right[i]]);

	// create sons (child elements keep the orientation of the parent element)
	unsigned int son[2][Hex::NUM_VERTICES]; // two hex elements
	for (int i = 0; i < 4; i++) {
		son[0][left[i]] = vtx[left[i]];
		son[0][right[i]] = son[1][left[i]] = mp[i];
		son[1][right[i]] = vtx[right[i]];
	}

	parent->active = false; // make parent element inactive
	parent->unref_all_nodes();
	unref_edges(parent);
	for (int i = 0; i < 2; i++) { // add child elements
		Hex *hex = create_hex(son[i]);
		parent->sons[i] = hex->id;
		hex->active = true;
		hex->marker = parent->marker;
		ref_edges(hex);
	}
	nactive++; // one new active element

	// adjust facets
	int face_0[2]; // two faces that are not refined
	int face_2[4]; // four faces that are splitted
	int face_refts[4]; // refinements that are applied on faces in face_2 array
	switch (refinement) {
		case H3D_REFT_HEX_X: {
			int fr[] = { H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT };
			memcpy(face_refts, fr, sizeof(fr));
			int f0[] = { 0, 1 };
			int f2[] = { 2, 3, 4, 5 };
			memcpy(face_0, f0, sizeof(f0));
			memcpy(face_2, f2, sizeof(f2));
			} break;

		case H3D_REFT_HEX_Y: {
			int fr[] = { H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ };
			memcpy(face_refts, fr, sizeof(fr));
			int f0[] = { 2, 3 };
			int f2[] = { 0, 1, 4, 5 };
			memcpy(face_0, f0, sizeof(f0));
			memcpy(face_2, f2, sizeof(f2));
			} break;

		case H3D_REFT_HEX_Z: {
			int fr[] = { H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ };
			memcpy(face_refts, fr, sizeof(fr));
			int f0[] = { 4, 5 };
			int f2[] = { 0, 1, 2, 3 };
			memcpy(face_0, f0, sizeof(f0));
			memcpy(face_2, f2, sizeof(f2));
			} break;
	}

	for (int i = 0; i < 4; i++)
		refined &= refine_quad_facet(parent, face_2[i], face_refts[i], parent->sons[0], parent->sons[1]);

	refined &= refine_quad_facet(parent, face_0[0], H3D_REFT_FACE_NONE, parent->sons[0]);
	refined &= refine_quad_facet(parent, face_0[1], H3D_REFT_FACE_NONE, parent->sons[1]);

	// add a facet between two new elements
	add_quad_facet(Facet::INNER, parent->sons[0], face_0[1], parent->sons[1], face_0[0]);

	return refined;
}

bool Mesh::refine_hex_4(Hex *parent, int refinement) {
	_F_
	bool refined = true;

	unsigned int vtx[Hex::NUM_VERTICES]; // vertices of parent element
	parent->get_vertices(vtx);

	const int *left, *right; // vertex indices for "left" and "right" face
	switch (refinement) {
		case H3D_H3D_REFT_HEX_XY:
			left = hex_face_vtcs_4;
			right = hex_face_vtcs_5;
			break;

		case H3D_H3D_REFT_HEX_YZ:
			left = hex_face_vtcs_0;
			right = hex_face_vtcs_1;
			break;

		case H3D_H3D_REFT_HEX_XZ:
			left = hex_face_vtcs_2;
			right = hex_face_vtcs_3;
			break;
	}

	unsigned int vp[2][4]; // vertex points on two facing faces ;)
	for (int i = 0; i < 4; i++) {
		vp[0][i] = vtx[left[i]];
		vp[1][i] = vtx[right[i]];
	}

	unsigned int emp[2][4]; // four edge midpoints on both facing faces
	for (int face = 0; face < 2; face++) {
		for (int i = 0; i < 4; i++)
			emp[face][i] = get_midpoint(vp[face][i], vp[face][(i + 1) % 4]);
	}

	unsigned int fmp[2]; // one midpoint on two facing faces
	for (int face = 0; face < 2; face++) {
		fmp[face] = get_midpoint(emp[face][0], emp[face][2]);
		set_midpoint(emp[face][1], emp[face][3], fmp[face]);
	}

	// sons (child elements keep the orientation of the parent element)
	unsigned int son[4][Hex::NUM_VERTICES]; // four hex elements
	switch (refinement) {
		case H3D_H3D_REFT_HEX_XY: {
			unsigned int s[4][Hex::NUM_VERTICES] = {
				{ vp[0][0], emp[0][0], fmp[0], emp[0][3], vp[1][0], emp[1][0], fmp[1], emp[1][3] },
				{ emp[0][0], vp[0][1], emp[0][1], fmp[0], emp[1][0], vp[1][1], emp[1][1], fmp[1] },
				{ fmp[0], emp[0][1], vp[0][2], emp[0][2], fmp[1], emp[1][1], vp[1][2], emp[1][2] },
				{ emp[0][3], fmp[0], emp[0][2], vp[0][3], emp[1][3], fmp[1], emp[1][2], vp[1][3] }
			};
			memcpy(son, s, sizeof(s));
			} break;

		case H3D_H3D_REFT_HEX_YZ: {
			unsigned int s[4][Hex::NUM_VERTICES] = {
				{ vp[0][0], vp[1][0], emp[1][0], emp[0][0], emp[0][3], emp[1][3], fmp[1], fmp[0] },
				{ emp[0][0], emp[1][0], vp[1][1], vp[0][1], fmp[0], fmp[1], emp[1][1], emp[0][1] },
				{ fmp[0], fmp[1], emp[1][1], emp[0][1], emp[0][2], emp[1][2], vp[1][2], vp[0][2] },
				{ emp[0][3], emp[1][3], fmp[1], fmp[0], vp[0][3], vp[1][3], emp[1][2], emp[0][2] }
			};
			memcpy(son, s, sizeof(s));
			} break;

		case H3D_H3D_REFT_HEX_XZ: {
			unsigned int s[4][Hex::NUM_VERTICES] = {
				{ vp[0][0], emp[0][0], emp[1][0], vp[1][0], emp[0][3], fmp[0], fmp[1], emp[1][3] },
				{ emp[0][0], vp[0][1], vp[1][1], emp[1][0], fmp[0], emp[0][1], emp[1][1], fmp[1] },
				{ fmp[0], emp[0][1], emp[1][1], fmp[1], emp[0][2], vp[0][2], vp[1][2], emp[1][2] },
				{ emp[0][3], fmp[0], fmp[1], emp[1][3], vp[0][3], emp[0][2], emp[1][2], vp[1][3] }
			};
			memcpy(son, s, sizeof(s));
			} break;
	}

	parent->active = false; // make parent element inactive
	parent->unref_all_nodes();
	unref_edges(parent);
	for (int i = 0; i < 4; i++) { // add child elements
		Hex *hex = create_hex(son[i]);
		parent->sons[i] = hex->id;
		hex->active = true;
		hex->marker = parent->marker;
		ref_edges(hex);
	}
	nactive += 3; // three new active element

	// adjust facets
	int face_4[2]; // two faces that are splitted into four facets
	int face_2[4]; // four faces that are splitted into two facets
	int face_refts[4]; // refinements that are applied on faces in face_2 array
	switch (refinement) {
		case H3D_H3D_REFT_HEX_XY: {
			int fr[] = { H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT };
			memcpy(face_refts, fr, sizeof(fr));
			int f4[] = { 4, 5 };
			int f2[] = { 0, 1, 2, 3 };
			memcpy(face_2, f2, sizeof(f2));
			memcpy(face_4, f4, sizeof(f4));
			} break;

		case H3D_H3D_REFT_HEX_YZ: {
			int fr[] = { H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ };
			memcpy(face_refts, fr, sizeof(fr));
			int f4[] = { 0, 1 };
			int f2[] = { 2, 3, 4, 5 };
			memcpy(face_2, f2, sizeof(f2));
			memcpy(face_4, f4, sizeof(f4));
			} break;

		case H3D_H3D_REFT_HEX_XZ: {
			int fr[] = { H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_HORZ, H3D_REFT_QUAD_VERT, H3D_REFT_QUAD_VERT };
			memcpy(face_refts, fr, sizeof(fr));
			int f4[] = { 2, 3 };
			int f2[] = { 0, 1, 4, 5 };
			memcpy(face_2, f2, sizeof(f2));
			memcpy(face_4, f4, sizeof(f4));
			} break;
	}

	int eidx[4][2] = { { 0, 3 }, { 1, 2 }, { 0, 1 }, { 3, 2 } };
	for (int i = 0; i < 4; i++)
		refined &= refine_quad_facet(parent, face_2[i], face_refts[i], parent->sons[eidx[i][0]], parent->sons[eidx[i][1]]);

	refined &= refine_quad_facet(parent, face_4[0], H3D_REFT_QUAD_BOTH, parent->sons[0], parent->sons[1], parent->sons[2], parent->sons[3]);
	refined &= refine_quad_facet(parent, face_4[1], H3D_REFT_QUAD_BOTH, parent->sons[0], parent->sons[1], parent->sons[2], parent->sons[3]);

	// add facet between two new elements
	add_quad_facet(Facet::INNER, parent->sons[0], face_2[1], parent->sons[1], face_2[0]);
	add_quad_facet(Facet::INNER, parent->sons[3], face_2[1], parent->sons[2], face_2[0]);
	add_quad_facet(Facet::INNER, parent->sons[0], face_2[3], parent->sons[3], face_2[2]);
	add_quad_facet(Facet::INNER, parent->sons[1], face_2[3], parent->sons[2], face_2[2]);

	return refined;
}

bool Mesh::refine_hex_8(Hex *parent, int refinement) {
	_F_
	bool refined = true;

	unsigned int vtx[Hex::NUM_VERTICES]; // vertices of parent element
	parent->get_vertices(vtx);

	unsigned int emp[12]; // 12 midpoints on edges
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++) {
		const int *edge_vtx_idx = RefHex::get_edge_vertices(edge);
		emp[edge] = get_midpoint(vtx[edge_vtx_idx[0]], vtx[edge_vtx_idx[1]]);
	}

	unsigned int fmp[6]; // 6 midpoints on faces
	for (int face = 0; face < Hex::NUM_FACES; face++) {
		const int *face_edge_idx = RefHex::get_face_edges(face);
		fmp[face] = get_midpoint(emp[face_edge_idx[0]], emp[face_edge_idx[2]]);
		set_midpoint(emp[face_edge_idx[1]], emp[face_edge_idx[3]], fmp[face]);
	}

	unsigned int ctr = get_midpoint(fmp[0], fmp[1]); // midpoint in the center
	set_midpoint(fmp[2], fmp[3], ctr);
	set_midpoint(fmp[4], fmp[5], ctr);

	// sons (child elements keep the orientation of the parent element)
	unsigned int son[8][Hex::NUM_VERTICES] = { // eight hex elements
	    { vtx[0], emp[0], fmp[4], emp[3], emp[4], fmp[2], ctr, fmp[0] },
	    { emp[0], vtx[1], emp[1], fmp[4], fmp[2], emp[5], fmp[1], ctr },
	    { fmp[4], emp[1], vtx[2], emp[2], ctr, fmp[1], emp[6], fmp[3] },
	    { emp[3], fmp[4], emp[2], vtx[3], fmp[0], ctr, fmp[3], emp[7] },
	    { emp[4], fmp[2], ctr, fmp[0],    vtx[4], emp[8], fmp[5], emp[11] },
	    { fmp[2], emp[5], fmp[1], ctr,    emp[8], vtx[5], emp[9], fmp[5] },
	    { ctr, fmp[1], emp[6], fmp[3],    fmp[5], emp[9], vtx[6], emp[10] },
	    { fmp[0], ctr, fmp[3], emp[7],    emp[11], fmp[5], emp[10], vtx[7] }
	   };

	parent->active = false; // make parent element inactive
	parent->unref_all_nodes();
	unref_edges(parent);
	for (int i = 0; i < 8; i++) { // add child elements
		Hex *hex = create_hex(son[i]);
		parent->sons[i] = hex->id;
		hex->active = true;
		hex->marker = parent->marker;
		ref_edges(hex);
	}
	nactive += 7; // seven new active element

	// adjust facets
	refined &= refine_quad_facet(parent, 0, H3D_REFT_QUAD_BOTH, parent->sons[0], parent->sons[3], parent->sons[7], parent->sons[4]);
	refined &= refine_quad_facet(parent, 1, H3D_REFT_QUAD_BOTH, parent->sons[1], parent->sons[2], parent->sons[6], parent->sons[5]);
	refined &= refine_quad_facet(parent, 2, H3D_REFT_QUAD_BOTH, parent->sons[0], parent->sons[1], parent->sons[5], parent->sons[4]);
	refined &= refine_quad_facet(parent, 3, H3D_REFT_QUAD_BOTH, parent->sons[3], parent->sons[2], parent->sons[6], parent->sons[7]);
	refined &= refine_quad_facet(parent, 4, H3D_REFT_QUAD_BOTH, parent->sons[0], parent->sons[1], parent->sons[2], parent->sons[3]);
	refined &= refine_quad_facet(parent, 5, H3D_REFT_QUAD_BOTH, parent->sons[4], parent->sons[5], parent->sons[6], parent->sons[7]);

	// add new facets
	add_quad_facet(Facet::INNER, parent->sons[0], 1, parent->sons[1], 0);
	add_quad_facet(Facet::INNER, parent->sons[3], 1, parent->sons[2], 0);
	add_quad_facet(Facet::INNER, parent->sons[4], 1, parent->sons[5], 0);
	add_quad_facet(Facet::INNER, parent->sons[7], 1, parent->sons[6], 0);

	add_quad_facet(Facet::INNER, parent->sons[0], 3, parent->sons[3], 2);
	add_quad_facet(Facet::INNER, parent->sons[1], 3, parent->sons[2], 2);
	add_quad_facet(Facet::INNER, parent->sons[4], 3, parent->sons[7], 2);
	add_quad_facet(Facet::INNER, parent->sons[5], 3, parent->sons[6], 2);

	add_quad_facet(Facet::INNER, parent->sons[0], 5, parent->sons[4], 4);
	add_quad_facet(Facet::INNER, parent->sons[1], 5, parent->sons[5], 4);
	add_quad_facet(Facet::INNER, parent->sons[2], 5, parent->sons[6], 4);
	add_quad_facet(Facet::INNER, parent->sons[3], 5, parent->sons[7], 4);

	return refined;
}

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, unsigned int eid) {
	_F_
	assert(face_refinement == H3D_REFT_FACE_NONE);

	Facet::Key fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets[fid];
	assert(facet->mode == HERMES_MODE_QUAD);

	//	if (is_compatible_quad_refinement(facet, face_refinement)) {
	if ((unsigned) facet->left == parent_elem->id) facet->set_left_info(eid, iface);
	else if ((unsigned) facet->right == parent_elem->id) facet->set_right_info(eid, iface);
	else assert(false); // Refining facet that does not face with appropriate element/boundary

	return true;
}

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, unsigned int eid0, unsigned int eid1) {
	_F_
	assert(face_refinement == H3D_REFT_QUAD_HORZ || face_refinement == H3D_REFT_QUAD_VERT);

	Facet::Key fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets[fid];
	assert(facet->mode == HERMES_MODE_QUAD);
	if ((unsigned) facet->type == Facet::INNER && (unsigned) facet->left == parent_elem->id) {
		// refine to the left
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == face_refinement) {
			facet->lactive = false;
			facet->ref_mask = face_refinement;

			Facet *f0 = add_quad_facet(Facet::INNER, eid0, facet->left_face_num, INVALID_IDX, -1);
			Facet *f1 = add_quad_facet(Facet::INNER, eid1, facet->left_face_num, INVALID_IDX, -1);

			f0->parent = fid;
			f1->parent = fid;

			if (face_refinement == H3D_REFT_QUAD_HORZ) {
				facet->sons[0] = get_facet_id(elements[eid0], facet->left_face_num);
				facet->sons[1] = get_facet_id(elements[eid1], facet->left_face_num);
			}
			else {
				facet->sons[2] = get_facet_id(elements[eid0], facet->left_face_num);
				facet->sons[3] = get_facet_id(elements[eid1], facet->left_face_num);
			}
		}
		else if (face_refinement == H3D_REFT_QUAD_HORZ) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, eid1, facet->left_face_num, INVALID_IDX, -1);
			Facet *lower_facet = add_quad_facet(Facet::INNER, eid0, facet->left_face_num, INVALID_IDX, -1);

			upper_facet->parent = fid;
			upper_facet->ractive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = facet->sons[3];
			upper_facet->sons[3] = facet->sons[2];
      Facet::Key upper_id = get_facet_id(elements[eid1], facet->left_face_num);

			lower_facet->parent = fid;
			lower_facet->ractive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = facet->sons[0];
			lower_facet->sons[3] = facet->sons[1];
			Facet::Key lower_id = get_facet_id(elements[eid0], facet->left_face_num);

			facets[facet->sons[0]]->parent = facets[facet->sons[1]]->parent = lower_id;
			facets[facet->sons[3]]->parent = facets[facet->sons[2]]->parent = upper_id;

			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_HORZ;
			facet->sons[0] = upper_id;
			facet->sons[1] = lower_id;
      facet->sons[2] = facet->sons[3] = Facet::invalid_key;
		}
		else if (face_refinement == H3D_REFT_QUAD_VERT) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, eid1, facet->left_face_num, INVALID_IDX, -1);
			Facet *lower_facet = add_quad_facet(Facet::INNER, eid0, facet->left_face_num, INVALID_IDX, -1);

			upper_facet->parent = fid;
			upper_facet->ractive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = facet->sons[2];
			upper_facet->sons[0] = facet->sons[1];
			Facet::Key upper_id = get_facet_id(elements[eid1], facet->left_face_num);

			lower_facet->parent = fid;
			lower_facet->ractive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = facet->sons[3];
			lower_facet->sons[0] = facet->sons[0];
			Facet::Key lower_id = get_facet_id(elements[eid0], facet->left_face_num);

			facets[facet->sons[1]]->parent = facets[facet->sons[2]]->parent = upper_id;
			facets[facet->sons[0]]->parent = facets[facet->sons[3]]->parent = lower_id;

			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_VERT;
			facet->sons[2] = lower_id;
			facet->sons[3] = upper_id;
      facet->sons[0] = facet->sons[1] = Facet::invalid_key;
		}
		else
			EXIT("Trying to apply incompatible face refinement to element #%d.", parent_elem->id);
	}
	else if ((unsigned) facet->type == Facet::INNER && (unsigned) facet->right == parent_elem->id) {
		// refine to the right (with element)
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == face_refinement) {
			facet->ractive = false;
			facet->ref_mask = face_refinement;

			Facet *f0 = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid0, facet->right_face_num);
			Facet *f1 = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid1, facet->right_face_num);

			f0->parent = fid;
			f1->parent = fid;

			if (face_refinement == H3D_REFT_QUAD_HORZ) {
				facet->sons[0] = get_facet_id(elements[eid0], facet->right_face_num);
				facet->sons[1] = get_facet_id(elements[eid1], facet->right_face_num);
			}
			else {
				facet->sons[2] = get_facet_id(elements[eid0], facet->right_face_num);
				facet->sons[3] = get_facet_id(elements[eid1], facet->right_face_num);
			}
		}
		else if (face_refinement == H3D_REFT_QUAD_HORZ) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid1, facet->right_face_num);
			Facet *lower_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid0, facet->right_face_num);

			upper_facet->parent = fid;
			upper_facet->lactive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = facet->sons[3];
			upper_facet->sons[3] = facet->sons[2];
			Facet::Key upper_id = get_facet_id(elements[eid1], facet->right_face_num);

			lower_facet->parent = fid;
			lower_facet->lactive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = facet->sons[0];
			lower_facet->sons[3] = facet->sons[1];
			Facet::Key lower_id = get_facet_id(elements[eid0], facet->right_face_num);

			facets[facet->sons[0]]->parent = facets[facet->sons[1]]->parent = lower_id;
			facets[facet->sons[3]]->parent = facets[facet->sons[2]]->parent = upper_id;

			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_HORZ;
			facet->sons[1] = upper_id;
			facet->sons[0] = lower_id;
      facet->sons[2] = facet->sons[3] = Facet::invalid_key;
		}
		else if (face_refinement == H3D_REFT_QUAD_VERT) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid1, facet->right_face_num);
			Facet *lower_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid0, facet->right_face_num);

			upper_facet->parent = fid;
			upper_facet->lactive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = facet->sons[2];
			upper_facet->sons[0] = facet->sons[1];
			Facet::Key upper_id = get_facet_id(elements[eid1], facet->right_face_num);

			lower_facet->parent = fid;
			lower_facet->lactive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = facet->sons[3];
			lower_facet->sons[0] = facet->sons[0];
			Facet::Key lower_id = get_facet_id(elements[eid0], facet->right_face_num);

			facets[facet->sons[1]]->parent = facets[facet->sons[2]]->parent = upper_id;
			facets[facet->sons[0]]->parent = facets[facet->sons[3]]->parent = lower_id;

			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_VERT;
			facet->sons[2] = lower_id;
			facet->sons[3] = upper_id;
			facet->sons[0] = facet->sons[1] = Facet::invalid_key;
		}
		else
			EXIT("Trying to apply incompatible face refinement to element #%d.", parent_elem->id);
	}
	else if (facet->type == Facet::OUTER) {
		// refine to the right (with boundary)
		facet->lactive = false;
		facet->ractive = false;
		facet->ref_mask = face_refinement;

		Facet *f0 = add_quad_facet(Facet::OUTER, eid0, facet->left_face_num, facet->right, facet->right_face_num);
		Facet *f1 = add_quad_facet(Facet::OUTER, eid1, facet->left_face_num, facet->right, facet->right_face_num);

		f0->parent = fid;
		f1->parent = fid;

		if (face_refinement == H3D_REFT_QUAD_HORZ) {
			facet->sons[0] = get_facet_id(elements[eid0], facet->left_face_num);
			facet->sons[1] = get_facet_id(elements[eid1], facet->left_face_num);
		}
		else {
			facet->sons[2] = get_facet_id(elements[eid0], facet->left_face_num);
			facet->sons[3] = get_facet_id(elements[eid1], facet->left_face_num);
		}
	}
	else {
		EXIT("Refining facet that does not face with appropriate element/boundary");
	}

	return true;
}

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, unsigned int eid0, unsigned int eid1, unsigned int eid2, unsigned int eid3) {
	_F_
	assert(face_refinement == H3D_REFT_QUAD_BOTH);

	Facet::Key fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets[fid];
	assert(facet->mode == HERMES_MODE_QUAD);

	//	if (is_compatible_quad_refinement(facet, face_refinement)) {
	if ((unsigned) facet->type == Facet::INNER && (unsigned) facet->left == parent_elem->id) {
		// refine to the left
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == H3D_REFT_QUAD_BOTH) {
			// no refinement OR the same type of the refinement on the right side
			// -> refinement on the left side is safe
			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_BOTH;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			for (int i = 0; i < 4; i++) {
				Facet *child_facet = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);
				child_facet->parent = fid;
				facet->sons[i] = get_facet_id(elements[ei[i]], facet->left_face_num);
			}
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_HORZ) { // FIXME: ignoring the orientation
			facet->lactive = false;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);

			Facet *upper_facet = facets[facet->sons[1]];
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = get_facet_id(elements[ei[3]], facet->left_face_num);
			upper_facet->sons[3] = get_facet_id(elements[ei[2]], facet->left_face_num);
			child_facets[2]->parent = child_facets[3]->parent = get_facet_id(elements[upper_facet->right], upper_facet->right_face_num);

			Facet *lower_facet = facets[facet->sons[0]];
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = get_facet_id(elements[ei[0]], facet->left_face_num);
			lower_facet->sons[3] = get_facet_id(elements[ei[1]], facet->left_face_num);
			child_facets[0]->parent = child_facets[1]->parent = get_facet_id(elements[lower_facet->right], lower_facet->right_face_num);
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_VERT) { // FIXME: ignoring the orientation
			facet->lactive = false;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);

			Facet *upper_facet = facets[facet->sons[3]];
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = get_facet_id(elements[ei[2]], facet->left_face_num);
			upper_facet->sons[0] = get_facet_id(elements[ei[1]], facet->left_face_num);
			child_facets[2]->parent = child_facets[1]->parent = get_facet_id(elements[upper_facet->right], upper_facet->right_face_num);

			Facet *lower_facet = facets[facet->sons[2]];
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = get_facet_id(elements[ei[3]], facet->left_face_num);
			lower_facet->sons[0] = get_facet_id(elements[ei[0]], facet->left_face_num);
			child_facets[0]->parent = child_facets[3]->parent = get_facet_id(elements[lower_facet->right], lower_facet->right_face_num);
		}
	}
	else if ((unsigned) facet->type == Facet::INNER && (unsigned) facet->right == parent_elem->id) {
		// refine to the right (with element)
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == H3D_REFT_QUAD_BOTH) {
			// no refinement OR the same type of the refinement on the left side
			// -> refinement on the right side is safe
			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_BOTH;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			for (int i = 0; i < 4; i++) {
				Facet *child_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);
				child_facet->parent = fid;
				facet->sons[i] = get_facet_id(elements[ei[i]], facet->right_face_num);
			}

		}
		else if (facet->ref_mask == H3D_REFT_QUAD_HORZ) { // FIXME: ignoring the orientation
			facet->ractive = false;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);

			Facet *upper_facet = facets[facet->sons[1]];
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = get_facet_id(elements[ei[3]], facet->right_face_num);
			upper_facet->sons[3] = get_facet_id(elements[ei[2]], facet->right_face_num);
			child_facets[2]->parent = child_facets[3]->parent = get_facet_id(elements[upper_facet->left], upper_facet->left_face_num);

			Facet *lower_facet = facets[facet->sons[0]];
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = get_facet_id(elements[ei[0]], facet->right_face_num);
			lower_facet->sons[3] = get_facet_id(elements[ei[1]], facet->right_face_num);
			child_facets[0]->parent = child_facets[1]->parent = get_facet_id(elements[lower_facet->left], lower_facet->left_face_num);
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_VERT) { // FIXME: ignoring the orientation
			facet->ractive = false;

			unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);

			Facet *upper_facet = facets[facet->sons[3]];
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = get_facet_id(elements[ei[2]], facet->right_face_num);
			upper_facet->sons[0] = get_facet_id(elements[ei[1]], facet->right_face_num);
			child_facets[2]->parent = child_facets[1]->parent = get_facet_id(elements[upper_facet->left], upper_facet->left_face_num);

			Facet *lower_facet = facets[facet->sons[2]];
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = get_facet_id(elements[ei[3]], facet->right_face_num);
			lower_facet->sons[0] = get_facet_id(elements[ei[0]], facet->right_face_num);
			child_facets[0]->parent = child_facets[3]->parent = get_facet_id(elements[lower_facet->left], lower_facet->left_face_num);
		}
	}
	else if (facet->type == Facet::OUTER) {
		// refine to the right (with boundary)
		facet->lactive = false;
		facet->ractive = false;
		facet->ref_mask = H3D_REFT_QUAD_BOTH;

		unsigned int ei[4] = { eid0, eid1, eid2, eid3 };
		for (int i = 0; i < 4; i++) {
			Facet *child_facet = add_quad_facet(Facet::OUTER, ei[i], facet->left_face_num, facet->right, facet->right_face_num);
			child_facet->parent = fid;
			facet->sons[i] = get_facet_id(elements[ei[i]], facet->left_face_num);
		}
	}
	else {
		EXIT("Refining facet that does not face with appropriate element/boundary");
	}

	return true;
}

Facet *Mesh::add_quad_facet(Facet::Type type, unsigned int left_elem, int left_iface, unsigned int right_elem, int right_iface) {
	_F_
	unsigned int elem_id;
	int iface;
	if (left_elem != INVALID_IDX) {
		elem_id = left_elem;
		iface = left_iface;
	}
	else if (right_elem != INVALID_IDX) {
		elem_id = right_elem;
		iface = right_iface;
	}
	else
		assert(false);
  
  Facet* facet = NULL;
  Facet::Key fidx = get_facet_id(elements[elem_id], iface);
  if (facets.find(fidx) != facets.end()) {
		// update info on existing facet
		facet = facets[fidx];
		if (elem_id == left_elem) facet->set_left_info(left_elem, left_iface);
		else facet->set_right_info(right_elem, right_iface);
	}
	else {
		// create new facet
		facet = new Facet(HERMES_MODE_QUAD);
		MEM_CHECK(facet);
		facet->type = type;
		facet->set_left_info(left_elem, left_iface);
		facet->set_right_info(right_elem, right_iface);
	}

	unsigned int facet_idxs[Quad::NUM_VERTICES];
	Element *e;
	if (left_elem != INVALID_IDX) {
		e = elements[left_elem];
		e->get_face_vertices(left_iface, facet_idxs);
	}
	else if (right_elem != INVALID_IDX) {
		e = elements[right_elem];
		e->get_face_vertices(right_iface, facet_idxs);
	}
	else EXIT("Right and left elements on facet are not set.");

	// update bnd mark for boundary edges
	if (facet->type == Facet::OUTER) {
		for (int i = 0; i < Quad::NUM_EDGES; i++) {
			unsigned int vtx[Edge::NUM_VERTICES] = { facet_idxs[i % Quad::NUM_EDGES], facet_idxs[(i + 1) % Quad::NUM_EDGES] };

      Edge::Key edge_key(vtx + 0, Edge::NUM_VERTICES);
      if(edges.find(edge_key) == edges.end())
        edges[edge_key] = new Edge;
			edges[edge_key]->bnd = 1;
		}
	}

	this->facets[fidx] = facet;

	return facet;
}

void Mesh::refine_all_elements(int refinement) {
	_F_
  std::map<unsigned int, Element*> local_elements = elements;
	for(std::map<unsigned int, Element*>::iterator it = local_elements.begin(); it != local_elements.end(); it++)
		if (it->second->used && it->second->active)
      refine_element(it->first, refinement);
}

void Mesh::refine_by_criterion(int(*criterion)(Element* e), int depth) {
	// TODO: implement me
}

void Mesh::unrefine_element(unsigned int id) {
	// TODO: implement me
}

void Mesh::unrefine_all_elements() {
	// TODO: implement me
}

unsigned int Mesh::create_midpoint(unsigned int a, unsigned int b) {
	_F_
	// get vertices
	Vertex *v1 = vertices.at(a);
	Vertex *v2 = vertices.at(b);
	// create middle point
	return add_vertex((v1->x + v2->x) / 2.0, (v1->y + v2->y) / 2.0, (v1->z + v2->z) / 2.0);
}

unsigned int Mesh::get_midpoint(unsigned int a, unsigned int b) {
	_F_
	unsigned int idx = peek_midpoint(a, b);
	if (idx == INVALID_IDX) {
    idx = create_midpoint(a, b);
		MidPointKey key(a, b);
		midpoints[key] = idx;
	}
	return idx;
}

unsigned int Mesh::peek_midpoint(unsigned int a, unsigned int b) const {
	_F_
  MidPointKey key(a, b);
	unsigned int idx = INVALID_IDX;
  if(midpoints.find(key) != midpoints.end())
    idx = midpoints.find(key)->second;
	return idx;
}

void Mesh::set_midpoint(unsigned int a, unsigned int b, unsigned int idx) {
	_F_
  MidPointKey key(a, b);
  midpoints[key] = idx;
}

Edge::Key Mesh::get_edge_id(unsigned int a, unsigned int b) const {
	_F_
	unsigned int pt[] = { a, b };
  return Edge::Key(pt + 0, Edge::NUM_VERTICES);
}

/// referencing edges
void Mesh::ref_edges(Element *e) 
{
  _F_
  assert(e != NULL);

  for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
       unsigned int vtx[Edge::NUM_VERTICES];
       e->get_edge_vertices(iedge, vtx);

    Edge::Key key(vtx + 0, Edge::NUM_VERTICES);
    if (edges.find(key) != edges.end()) {
      edges.find(key)->second->ref++;
    }
    else {
      Edge* edge = new Edge;
      edge->ref = 1;
      edges[key] = edge;
    }
  }
}

void Mesh::unref_edges(Element *e) 
{
  _F_
  assert(e != NULL);

  for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
    unsigned int vtx[Edge::NUM_VERTICES];
    e->get_edge_vertices(iedge, vtx);

    Edge::Key key(vtx + 0, Edge::NUM_VERTICES);
    if (edges.find(key) != edges.end()) {
      edges.find(key)->second->ref--;
    }
    else assert(false); // Unreferencing non-existent edge
  }
}

Facet::Key Mesh::get_facing_facet(Facet::Key fid, unsigned int elem_id) {
	_F_
	Facet *facet = facets[fid];

	if (facet != NULL) {
		if (elem_id == (unsigned) facet->left) {
      while (!facet->ractive && facet->parent != Facet::invalid_key) {
				fid = facet->parent;
				facet = facets[fid];
			}
			return fid;
		}
		else if (elem_id == (unsigned) facet->right) {
			while (!facet->lactive && facet->parent != Facet::invalid_key) {
				fid = facet->parent;
				facet = facets[fid];
			}
			return fid;
		}
		else
			return Facet::invalid_key;
	}
	else
		return Facet::invalid_key;
}
  
void Mesh::regularize() {
	_F_
	// FIXME: implements only 1-irregularity rule (quite dirty hack this is)
	// Assumes only XYZ refinements of elements (i.e. no anisotropic refinements, no incompatible refinements)

	// NOTE: this is very stupid implementation. It checks the parent facet of a half-active facet. If this parent facet is
	// active on the opposite side, we found a hanging node of 1. order. If it is inactive, we check super parent (parent of
	// this parent facet) the same way. If it is active, we found hanging node of a  2. order and we refine this super parent.
	// If it is inactive, hanging node of a higher order was found and we report an error.

  for(std::map<unsigned int, Element*>::iterator it = elements.begin(); it != elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *elem = elements[it->first];
		  for (int iface = 0; iface < elem->get_num_faces(); iface++) {
        Facet::Key fid = get_facet_id(elem, iface);
			  Facet *facet = facets[fid];
			  assert(facet != NULL);
			  if (facet->lactive && !facet->ractive) {
				  if (facet->parent != Facet::invalid_key) {
            if(facets.find(facet->parent) != facets.end()) {
              Facet *parent = facets.find(facet->parent)->second;
					    if (parent->ractive) {
						    // OK: 1. order hanging node
					    }
					    else {
						    if (parent->parent != Facet::invalid_key) {
                  Facet *super_parent = facets.find(parent->parent)->second;
							    if (super_parent->ractive) {
								    refine_element(super_parent->right, H3D_H3D_H3D_REFT_HEX_XYZ);
							    }
							    else {
								    EXIT("Cannot handle hanging node of order > 1");
							    }
						    }
					    }
				    }
          }
			  }
			  else if (!facet->lactive && facet->ractive) {
				  if (facet->parent != Facet::invalid_key) {
            Facet *parent = facets.find(facet->parent)->second;
					  if (parent->lactive) {
						  // OK: 1. order hanging node
					  }
					  else {
						  if (parent->parent != Facet::invalid_key) {
                Facet *super_parent = facets.find(parent->parent)->second;
							  if (super_parent->lactive) {
								  refine_element(super_parent->left, H3D_H3D_H3D_REFT_HEX_XYZ);
							  }
							  else {
								  EXIT("Cannot handle hanging node of order > 1");
							  }
						  }
					  }
				  }
			  }
		  }
	  }
}

void Mesh::refine_towards_boundary(int marker, int depth) {
	_F_

	if (depth == 0) return;
  std::map<unsigned int, Element*> local_elements = elements;
	for(std::map<unsigned int, Element*>::iterator it = local_elements.begin(); it != local_elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = elements[it->first];

		  int split = H3D_SPLIT_NONE;
		  for (int iface = 0; iface < e->get_num_faces(); iface++) {
        Facet::Key fid = get_facet_id(e, iface);
			  Facet *facet = facets[fid];

			  if (facet->type == Facet::OUTER) {
				  Boundary *bnd = boundaries[facet->right];
				  if (bnd->marker == marker) {
					  if (iface == 0 || iface == 1) split |= H3D_SPLIT_HEX_X;
					  else if (iface == 2 || iface == 3) split |= H3D_SPLIT_HEX_Y;
					  else if (iface == 4 || iface == 5) split |= H3D_SPLIT_HEX_Z;
				  }
			  }
		  }

		  int reft[] = {
			  H3D_REFT_HEX_NONE, H3D_REFT_HEX_X, H3D_REFT_HEX_Y, H3D_H3D_REFT_HEX_XY,
			  H3D_REFT_HEX_Z, H3D_H3D_REFT_HEX_XZ, H3D_H3D_REFT_HEX_YZ, H3D_H3D_H3D_REFT_HEX_XYZ
		  };
      refine_element(it->first, reft[split]);
	  }

	refine_towards_boundary(marker, depth - 1);
}

void Mesh::check_elem_oris()
{
	_F_
	RefMap refmap(this);
  for(std::map<unsigned int, Element*>::iterator it = elements.begin(); it != elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = it->second;
		  refmap.set_active_element(e);

		  Ord3 ord;
		  if (e->get_mode() == HERMES_MODE_HEX) ord = Ord3(1, 1, 1);
		  else if (e->get_mode() == HERMES_MODE_TET) ord = Ord3(2);
		  else warning(HERMES_ERR_NOT_IMPLEMENTED);
		  Quad3D *quad = get_quadrature(e->get_mode());
		  int np = quad->get_num_points(ord);
		  QuadPt3D *pt = quad->get_points(ord);

		  double3x3 *m = refmap.get_ref_map(np, pt);
		  for (int i = 0; i < np; i++) {
        if (det(m[i]) <= 0) error("Element #%ld has an incorrect orientation.", it->first);
		  }
		  delete [] m;
	  }
}
