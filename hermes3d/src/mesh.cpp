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
#include <common/error.h>
#include <common/bitarray.h>
#include <common/callstack.h>

#include "mesh.h"
#include "meshloader.h"
#include "refdomain.h"
#include "refmap.h"
#include "determinant.h"

const int TOP_LEVEL_REF = -1;

// to print out the banner
static bool print_banner = true;
// forward declarations
extern bool verbose;
extern void banner();


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

Edge::Edge() {
	bnd = 0;
	ref = 0;
}

// for debugging
void Edge::dump() {
	printf("bnd = %d, ref = %d\n", bnd, ref);
}

// Facet //////////////////////////////////////////////////////////////////////

Facet::Facet(EMode2D mode) {
	_F_
	this->mode = mode;
	this->type = INNER;
	this->left = INVALID_IDX;
	this->right = INVALID_IDX;
	this->left_face_num = -1;
	this->right_face_num = -1;
	this->lactive = false;
	this->ractive = false;
	this->ref_mask = H3D_REFT_FACE_NONE;

	this->parent = INVALID_IDX;
	for (int i = 0; i < MAX_SONS; i++)
		this->sons[i] = INVALID_IDX;
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

bool Facet::ced(Word_t idx, int iface)
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

	printf("type = %s (%s), [%d, %d], left (elem = %ld, face = %d), ", s_type[type], s_mode[mode], lactive, ractive, left, left_face_num);
	if (type == INNER) printf(" right (elem = %ld, face = %d)", right, right_face_num);
	else printf(" right (bdr = %ld)", right);
	printf(", ref_mask = %d, sons = [%ld, %ld, %ld, %ld], ", ref_mask, sons[0], sons[1], sons[2], sons[3]);
	if (parent != INVALID_IDX) printf("parent = %ld", parent);
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
	printf("id = %ld\n", id);
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

Hex::Hex(Word_t v[]) {
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

Hex::Hex(Word_t v1, Word_t v2, Word_t v3, Word_t v4, Word_t v5, Word_t v6, Word_t v7, Word_t v8) {
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

int Hex::get_edge_vertices(int edge_num, Word_t *vtcs) const {
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

int Hex::get_face_vertices(int face_num, Word_t *vtcs) const {
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
	Word_t v[4];
	get_face_vertices(face_num, v);

	Word_t minval = 1000;
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
	printf("id = %ld (%d, %d, %d), vertices(%ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld), ", id, active, used, reft,
		vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
	printf("sons(%ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld), ",
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

Tetra::Tetra(Word_t v[]) {
	_F_
#ifdef WITH_TETRA
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

Tetra::Tetra(Word_t v1, Word_t v2, Word_t v3, Word_t v4) {
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

int Tetra::get_edge_vertices(int edge_num, Word_t *vtcs) const {
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

int Tetra::get_face_vertices(int face_num, Word_t *vtcs) const {
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

	Word_t v0 = vtcs[map[face_num][0]];
	Word_t v1 = vtcs[map[face_num][1]];
	Word_t v2 = vtcs[map[face_num][2]];

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
	printf("id = %ld, vertices(%ld, %ld, %ld, %ld), ", id, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
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

Prism::Prism(Word_t v[]) {
	_F_
#ifdef WITH_PRISM
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
#else
	EXIT(H3D_ERR_PRISM_NOT_COMPILED);
#endif
}

Prism::Prism(Word_t v1, Word_t v2, Word_t v3, Word_t v4, Word_t v5, Word_t v6) {
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

int Prism::get_edge_vertices(int edge_num, Word_t *vtcs) const {
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

EMode2D Prism::get_face_mode(int face_num) const {
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

int Prism::get_face_vertices(int face_num, Word_t *vtcs) const {
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
	EXIT(H3D_ERR_NOT_IMPLEMENTED);
	// FIXME
	return -1;
}

int Prism::get_face_orientation(int face_num) const {
	_F_
	EXIT(H3D_ERR_NOT_IMPLEMENTED);
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
	printf("id = %ld, vertices(%ld, %ld, %ld, %ld, %ld, %ld), ", id,
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
	printf("id = %ld, marker = %d\n", id, marker);
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

	if (print_banner) {
		banner();
		print_banner = false;
	}

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
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i))
		delete vertices[i];
	vertices.remove_all();

	for (Word_t i = boundaries.first(); i != INVALID_IDX; i = boundaries.next(i))
		delete boundaries[i];
	boundaries.remove_all();

	for (Word_t i = elements.first(); i != INVALID_IDX; i = elements.next(i))
		delete elements[i];
	elements.remove_all();

	for (Word_t i = facets.first(); i != INVALID_IDX; i = facets.next(i))
		delete facets.get(i);
	facets.remove_all();

	midpoints.remove_all();
	edges.remove_all();
}

void Mesh::copy(const Mesh &mesh) {
	_F_
	if (&mesh == this) warning("Copying mesh into itself.");
	free();

	// copy vertices
	for (Word_t i = mesh.vertices.first(); i != INVALID_IDX; i = mesh.vertices.next(i))
		this->vertices.set(i, mesh.vertices[i]->copy());

	// copy boundaries
	for (Word_t i = mesh.boundaries.first(); i != INVALID_IDX; i = mesh.boundaries.next(i))
		this->boundaries.set(i, mesh.boundaries[i]->copy());

	// copy elements, facets and edges
	for (Word_t i = mesh.elements.first(); i != INVALID_IDX; i = mesh.elements.next(i)) {
		Element *e = mesh.elements[i];
		this->elements.set(i, e->copy());

		// copy mid points on edges and edges
		Word_t emp[e->get_num_edges()];
		for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
			Word_t edge_vtx[Edge::NUM_VERTICES];
			e->get_edge_vertices(iedge, edge_vtx);
			emp[iedge] = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
			if (emp[iedge] != INVALID_IDX)
				set_midpoint(edge_vtx[0], edge_vtx[1], emp[iedge]);

			Edge edge;
			if (mesh.edges.lookup(edge_vtx + 0, Edge::NUM_VERTICES, edge))
				edges.set(edge_vtx, Edge::NUM_VERTICES, edge);
		}

		// copy mid points on faces
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			const int *edge = e->get_face_edges(iface);

			switch (e->get_mode()) {
				case MODE_HEXAHEDRON:
					// horz
					if (emp[edge[1]] != INVALID_IDX && emp[edge[3]] != INVALID_IDX) {
						Word_t edge_vtx[Edge::NUM_VERTICES] = { emp[edge[1]], emp[edge[3]] };

						Word_t fmp = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
						if (fmp != INVALID_IDX)
							set_midpoint(edge_vtx[0], edge_vtx[1], fmp);
					}
					// vert
					if (emp[edge[0]] != INVALID_IDX && emp[edge[2]] != INVALID_IDX) {
						Word_t edge_vtx[Edge::NUM_VERTICES] = { emp[edge[0]], emp[edge[2]] };

						Word_t fmp = mesh.peek_midpoint(edge_vtx[0], edge_vtx[1]);
						if (fmp != INVALID_IDX)
							set_midpoint(edge_vtx[0], edge_vtx[1], fmp);
					}
					break;

				case MODE_TETRAHEDRON:
				case MODE_PRISM:
					break;
			}
		}
	}

	// facets
	for (Word_t fid = mesh.facets.first(); fid != INVALID_IDX; fid = mesh.facets.next(fid)) {
		Facet *facet = mesh.facets[fid];

		Word_t face_idxs[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices
		if (facet->left != INVALID_IDX) {
			Element *left_e = mesh.elements[facet->left];
			int nvtcs = left_e->get_face_vertices(facet->left_face_num, face_idxs);
			this->facets.set(face_idxs + 0, nvtcs, facet->copy());
		}
		else if (facet->right != INVALID_IDX && facet->type == Facet::INNER) {
			Element *right_e = mesh.elements[facet->right];
			int nvtcs = right_e->get_face_vertices(facet->right_face_num, face_idxs);
			this->facets.set(face_idxs + 0, nvtcs, facet->copy());
		}
		else
			EXIT("WTF?");		// FIXME
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
	for (Word_t eid = mesh.elements.first(); eid <= mesh.nbase; eid = mesh.elements.next(eid)) {
		Element *e = mesh.elements[eid];

		// vertices
		for (int iv = 0; iv < e->get_num_vertices(); iv++) {
			Word_t vtx = e->get_vertex(iv);
			if (!this->vertices.exists(vtx))
				this->vertices.set(vtx, mesh.vertices[vtx]->copy());
		}

		// edges
		for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
			Word_t vtx[Edge::NUM_VERTICES];
			e->get_edge_vertices(iedge, vtx);

			Edge edge;
			if (mesh.edges.lookup(vtx + 0, Edge::NUM_VERTICES, edge))
				this->edges.set(vtx, Edge::NUM_VERTICES, edge);
		}

		// facets
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Word_t face_idxs[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices

			int nvts = e->get_face_vertices(iface, face_idxs);
			Facet *facet = NULL;
			if (mesh.facets.lookup(face_idxs + 0, nvts, facet)) {
				if (!this->facets.lookup(face_idxs + 0, nvts, facet)) {
					// insert the facet
					Facet *fcopy = facet->copy_base();
					fcopy->left = eid;
					this->facets.set(face_idxs + 0, nvts, fcopy);

					// boundaries
					if (facet->type == Facet::OUTER) {
						Word_t bnd_id = facet->right;
						if (!this->boundaries.exists(bnd_id)) this->boundaries.set(bnd_id, mesh.boundaries[bnd_id]->copy());
					}
				}
				else {
					assert(facet->type == Facet::INNER);
					facet->set_right_info(eid, iface);
				}
			}
		}

		this->elements.set(eid, e->copy_base());
	}

	this->nbase = this->nactive = mesh.nbase;
	this->seq = g_mesh_seq++;
}

Word_t Mesh::get_facet_id(Element *e, int face_num) const {
	_F_
	assert(e != NULL);
	Word_t facet_idxs[Quad::NUM_VERTICES]; // quad is shape with the largest number of vertices
	int nvts = e->get_face_vertices(face_num, facet_idxs);
	return facets.get_idx(facet_idxs + 0, nvts);
}

Word_t Mesh::get_edge_id(Element *e, int edge_num) const {
	_F_
	assert(e != NULL);
	Word_t edge_idxs[Edge::NUM_VERTICES];
	int nvtcs = e->get_edge_vertices(edge_num, edge_idxs);
	return edges.get_idx(edge_idxs + 0, nvtcs);
}

void Mesh::dump() {
	_F_
	printf("Vertices (count = %ld)\n", vertices.count());
	for (Word_t i = vertices.first(); i != INVALID_IDX; i = vertices.next(i)) {
		Vertex *v = vertices[i];
		printf("  id = %ld, ", i);
		v->dump();
	}

	printf("Elements (count = %ld)\n", elements.count());
	for (Word_t i = elements.first(); i != INVALID_IDX; i = elements.next(i)) {
		Element *e = elements[i];
		printf("  ");
		e->dump();
	}

	printf("Boundaries (count = %ld)\n", boundaries.count());
	for (Word_t i = boundaries.first(); i != INVALID_IDX; i = boundaries.next(i)) {
		Boundary *b = boundaries[i];
		printf("  ");
		b->dump();
	}

	printf("Facets (count = %ld)\n", facets.count());
	for (Word_t i = facets.first(); i != INVALID_IDX; i = facets.next(i)) {
		Facet *f = facets.get(i);
		printf("  id = %ld, ", i);
		f->dump();
	}
}

Word_t Mesh::add_vertex(double x, double y, double z) {
	_F_
	Word_t idx = vertices.count() + 1;
	vertices.set(idx, new Vertex(x, y, z));
	return idx;
}

Tetra *Mesh::create_tetra(Word_t vtcs[]) {
	_F_
	Tetra *tetra = new Tetra(vtcs);
	MEM_CHECK(tetra);
	Word_t id = elements.count() + 1;
	elements.set(id, tetra);
	tetra->id = id;

	tetra->ref_all_nodes();

	return tetra;
}

Tetra *Mesh::add_tetra(Word_t vtcs[]) {
	_F_
	Tetra *tetra = create_tetra(vtcs);

	// edges
	ref_edges(tetra);

	// facets
	for (int i = 0; i < Tetra::NUM_FACES; i++) {
		Word_t facet_idxs[Tri::NUM_VERTICES];
		int nvtcs = tetra->get_face_vertices(i, facet_idxs);
		Facet *facet = NULL;
		if (facets.lookup(facet_idxs + 0, nvtcs, facet)) {
			facet->type = Facet::INNER;
			facet->set_right_info(tetra->id, i);
		}
		else {
			facet = new Facet(RefTetra::get_face_mode(i));
			facet->set_left_info(tetra->id, i);
			facets.set(facet_idxs + 0, nvtcs, facet);
		}
	}

	return tetra;
}

Hex *Mesh::create_hex(Word_t vtcs[]) {
	_F_
	// build up the element
	Hex *hex = new Hex(vtcs);
	MEM_CHECK(hex);
	Word_t id = elements.count() + 1;
	elements.set(id, hex);
	hex->id = id;

	hex->ref_all_nodes();

	return hex;
}

Hex *Mesh::add_hex(Word_t vtcs[]) {
	_F_
	Hex *hex = create_hex(vtcs);

	// edges
	ref_edges(hex);

	// facets
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		Word_t facet_idxs[Quad::NUM_VERTICES];
		int nvtcs = hex->get_face_vertices(i, facet_idxs);
		Facet *facet = NULL;
		if (facets.lookup(facet_idxs + 0, nvtcs, facet)) {
			facet->type = Facet::INNER;
			facet->set_right_info(hex->id, i);
		}
		else {
			Facet *fct = new Facet(MODE_QUAD);
			MEM_CHECK(fct);
			fct->set_left_info(hex->id, i);
			facets.set(facet_idxs + 0, nvtcs, fct);
		}
	}

	return hex;
}

Prism *Mesh::create_prism(Word_t vtcs[]) {
	_F_
	Prism *prism = new Prism(vtcs);
	MEM_CHECK(prism);
	Word_t id = elements.count() + 1;
	elements.set(id, prism);
	prism->id = id;

	prism->ref_all_nodes();

	return prism;
}

Prism *Mesh::add_prism(Word_t vtcs[]) {
	_F_
	Prism *prism = create_prism(vtcs);

	// edges
	ref_edges(prism);

	// facets
	for (int i = 0; i < Prism::NUM_FACES; i++) {
		Word_t facet_idxs[Quad::NUM_VERTICES];
		int nvtcs = prism->get_face_vertices(i, facet_idxs);
		Facet *facet = NULL;
		if (facets.lookup(facet_idxs + 0, nvtcs, facet)) {
			facet->type = Facet::INNER;
			facet->set_right_info(prism->id, i);
		}
		else {
			facet = new Facet(RefPrism::get_face_mode(i));
			MEM_CHECK(facet);
			facet->set_left_info(prism->id, i);
			facets.set(facet_idxs + 0, nvtcs, facet);
		}
	}

	return prism;
}

Boundary *Mesh::add_tri_boundary(Word_t vtcs[], int marker) {
	_F_
	Facet *facet = NULL;
	if (facets.lookup(vtcs + 0, Tri::NUM_VERTICES, facet)) {
		Boundary *bdr = new BoundaryTri(marker);
		MEM_CHECK(bdr);
		Word_t pos = boundaries.count() + 1;
		boundaries.set(pos, bdr);
		bdr->id = pos;

		facet->type = Facet::OUTER;
		facet->set_right_info(bdr->id);

		return bdr;
	}
	else
		return NULL;
}

Boundary *Mesh::add_quad_boundary(Word_t vtcs[], int marker) {
	_F_
	Facet *facet = NULL;
	if (facets.lookup(vtcs + 0, Quad::NUM_VERTICES, facet)) {
		Boundary *bdr = new BoundaryQuad(marker);
		MEM_CHECK(bdr);
		Word_t pos = boundaries.count() + 1;
		boundaries.set(pos, bdr);
		bdr->id = pos;

		facet->type = Facet::OUTER;
		facet->set_right_info(bdr->id);

		return bdr;
	}
	else
		return NULL;
}

void Mesh::ugh()
{
	_F_
	// set the number of active/base elements
	nactive = nbase = elements.count();

	// set bnd flag for boundary edges
	FOR_ALL_FACETS(idx, this){
		Facet *facet = facets[idx];
		if (facet->type == Facet::OUTER) {
			Element *elem = elements[facet->left];
			const int *face_edge = elem->get_face_edges(facet->left_face_num);
			for (int iedge = 0; iedge < elem->get_num_face_edges(facet->left_face_num); iedge++) {
				Word_t vtx[Edge::NUM_VERTICES];
				elem->get_edge_vertices(face_edge[iedge], vtx);

				Edge edge;
				edges.lookup(vtx + 0, Edge::NUM_VERTICES, edge);
				edge.bnd = 1;
				edges.set(vtx + 0, Edge::NUM_VERTICES, edge);
			}
		}
	}

	check_elem_oris();
}

bool Mesh::is_compatible_quad_refinement(Facet *facet, int reft) const {
	_F_
	if (facet->type == Facet::INNER) {
		// BOTH or NONE on the facet => refinement makes no problem
		if (facet->ref_mask == H3D_REFT_QUAD_BOTH || facet->ref_mask == H3D_REFT_FACE_NONE) return true;

		// applying BOTH or NONE is also no problem
		if (reft == H3D_REFT_QUAD_BOTH || reft == H3D_REFT_FACE_NONE) return true;

		Word_t eid;
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

		Element *e = elements[eid];
		int nv = e->get_num_face_vertices(face_num);
		Word_t face_vtx[nv];
		e->get_face_vertices(face_num, face_vtx);

		// check if the vertices are there, if so => compatible refinement
		Word_t emp[2] = { INVALID_IDX, INVALID_IDX };
		if (reft == H3D_REFT_QUAD_HORZ) {
			emp[0] = peek_midpoint(face_vtx[0], face_vtx[3]);
			emp[1] = peek_midpoint(face_vtx[1], face_vtx[2]);
		}
		else if (reft == H3D_REFT_QUAD_VERT) {
			emp[0] = peek_midpoint(face_vtx[0], face_vtx[1]);
			emp[1] = peek_midpoint(face_vtx[2], face_vtx[3]);
		}

		return (emp[0] != INVALID_IDX && emp[1] != INVALID_IDX);
	}
	else {
		// no problem with outer facets
		return true;
	}
}

bool Mesh::can_refine_element(Word_t eid, int reft) const {
	_F_
	bool can_refine = false;

	Element *elem = elements.get(eid);
	assert(elem != NULL);
	switch (elem->get_mode()) {
		case MODE_HEXAHEDRON: can_refine = can_refine_hex((Hex *) elem, reft); break;
		case MODE_TETRAHEDRON: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
		case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
		default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
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

		// WTF?
		default:
			EXIT(H3D_ERR_UNKNOWN_REFINEMENT_TYPE);
			break;
	}

	bool can_refine = true;
	for (int i = 0; i < nf; i++) {
		Word_t fid = get_facet_id(elem, iface[i]);
		Facet *facet = facets.get(fid);
		assert(facet != NULL);
		can_refine &= is_compatible_quad_refinement(facet, face_reft[i]);
	}

	return can_refine;

}

bool Mesh::refine_element(Word_t id, int refinement) {
	_F_
	bool refined = false;
	Element *elem = elements.get(id);
	assert(elem != NULL);
	if (can_refine_element(id, refinement)) {
		switch (elem->get_mode()) {
			case MODE_HEXAHEDRON: refined = refine_hex((Hex *) elem, refinement); break;
			case MODE_TETRAHEDRON: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
			case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
			default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
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
		// WTF?
		default:
			EXIT(H3D_ERR_UNKNOWN_REFINEMENT_TYPE);
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

	Word_t vtx[Hex::NUM_VERTICES]; // vertices of parent element
	parent->get_vertices(vtx);

	Word_t mp[4] = { 0 }; // four midpoints shared when refining to two elements
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
	Word_t son[2][Hex::NUM_VERTICES]; // two hex elements
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

	Word_t vtx[Hex::NUM_VERTICES]; // vertices of parent element
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

	Word_t vp[2][4]; // vertex points on two facing faces ;)
	for (int i = 0; i < 4; i++) {
		vp[0][i] = vtx[left[i]];
		vp[1][i] = vtx[right[i]];
	}

	Word_t emp[2][4]; // four edge midpoints on both facing faces
	for (int face = 0; face < 2; face++) {
		for (int i = 0; i < 4; i++)
			emp[face][i] = get_midpoint(vp[face][i], vp[face][(i + 1) % 4]);
	}

	Word_t fmp[2]; // one midpoint on two facing faces
	for (int face = 0; face < 2; face++) {
		fmp[face] = get_midpoint(emp[face][0], emp[face][2]);
		set_midpoint(emp[face][1], emp[face][3], fmp[face]);
	}

	// sons (child elements keep the orientation of the parent element)
	Word_t son[4][Hex::NUM_VERTICES]; // four hex elements
	switch (refinement) {
		case H3D_H3D_REFT_HEX_XY: {
			Word_t s[4][Hex::NUM_VERTICES] = {
				{ vp[0][0], emp[0][0], fmp[0], emp[0][3], vp[1][0], emp[1][0], fmp[1], emp[1][3] },
				{ emp[0][0], vp[0][1], emp[0][1], fmp[0], emp[1][0], vp[1][1], emp[1][1], fmp[1] },
				{ fmp[0], emp[0][1], vp[0][2], emp[0][2], fmp[1], emp[1][1], vp[1][2], emp[1][2] },
				{ emp[0][3], fmp[0], emp[0][2], vp[0][3], emp[1][3], fmp[1], emp[1][2], vp[1][3] }
			};
			memcpy(son, s, sizeof(s));
			} break;

		case H3D_H3D_REFT_HEX_YZ: {
			Word_t s[4][Hex::NUM_VERTICES] = {
				{ vp[0][0], vp[1][0], emp[1][0], emp[0][0], emp[0][3], emp[1][3], fmp[1], fmp[0] },
				{ emp[0][0], emp[1][0], vp[1][1], vp[0][1], fmp[0], fmp[1], emp[1][1], emp[0][1] },
				{ fmp[0], fmp[1], emp[1][1], emp[0][1], emp[0][2], emp[1][2], vp[1][2], vp[0][2] },
				{ emp[0][3], emp[1][3], fmp[1], fmp[0], vp[0][3], vp[1][3], emp[1][2], emp[0][2] }
			};
			memcpy(son, s, sizeof(s));
			} break;

		case H3D_H3D_REFT_HEX_XZ: {
			Word_t s[4][Hex::NUM_VERTICES] = {
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

	Word_t vtx[Hex::NUM_VERTICES]; // vertices of parent element
	parent->get_vertices(vtx);

	Word_t emp[12]; // 12 midpoints on edges
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++) {
		const int *edge_vtx_idx = RefHex::get_edge_vertices(edge);
		emp[edge] = get_midpoint(vtx[edge_vtx_idx[0]], vtx[edge_vtx_idx[1]]);
	}

	Word_t fmp[6]; // 6 midpoints on faces
	for (int face = 0; face < Hex::NUM_FACES; face++) {
		const int *face_edge_idx = RefHex::get_face_edges(face);
		fmp[face] = get_midpoint(emp[face_edge_idx[0]], emp[face_edge_idx[2]]);
		set_midpoint(emp[face_edge_idx[1]], emp[face_edge_idx[3]], fmp[face]);
	}

	Word_t ctr = get_midpoint(fmp[0], fmp[1]); // midpoint in the center
	set_midpoint(fmp[2], fmp[3], ctr);
	set_midpoint(fmp[4], fmp[5], ctr);

	// sons (child elements keep the orientation of the parent element)
	Word_t son[8][Hex::NUM_VERTICES] = { // eight hex elements
	    { vtx[0], emp[0], fmp[4], emp[3], emp[4], fmp[2], ctr, fmp[0] },
	    { emp[0], vtx[1], emp[1], fmp[4], fmp[2], emp[5], fmp[1], ctr },
	    { fmp[4], emp[1], vtx[2], emp[2], ctr, fmp[1], emp[6], fmp[3] },
	    { emp[3], fmp[4], emp[2], vtx[3], fmp[0], ctr, fmp[3], emp[7] },
	    { emp[4], fmp[2], ctr, fmp[0],    vtx[4], emp[8], fmp[5], emp[11] },
	    { fmp[2], emp[5], fmp[1], ctr,    emp[8], vtx[5], emp[9], fmp[5] },
	    { ctr, fmp[1], emp[6], fmp[3],    fmp[5], emp[9], vtx[6], emp[10] },
	    { fmp[0], ctr, fmp[3], emp[7],    emp[11], fmp[5], emp[10], vtx[7] }
	   };

	// deactivate edges on parent element
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		Word_t edge_vtx[Edge::NUM_VERTICES];
		parent->get_edge_vertices(i, edge_vtx);
	}

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

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, Word_t eid) {
	_F_
	assert(face_refinement == H3D_REFT_FACE_NONE);

	Word_t fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets.get(fid);
	assert(facet->mode == MODE_QUAD);

	//	if (is_compatible_quad_refinement(facet, face_refinement)) {
	if (facet->left == parent_elem->id) facet->set_left_info(eid, iface);
	else if (facet->right == parent_elem->id) facet->set_right_info(eid, iface);
	else assert(false); // Refining facet that does not face with appropriate element/boundary

	return true;
}

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, Word_t eid0, Word_t eid1) {
	_F_
	assert(face_refinement == H3D_REFT_QUAD_HORZ || face_refinement == H3D_REFT_QUAD_VERT);

	Word_t fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets.get(fid);
	assert(facet->mode == MODE_QUAD);
	if (facet->type == Facet::INNER && facet->left == parent_elem->id) {
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
			Word_t upper_id = get_facet_id(elements[eid1], facet->left_face_num);

			lower_facet->parent = fid;
			lower_facet->ractive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = facet->sons[0];
			lower_facet->sons[3] = facet->sons[1];
			Word_t lower_id = get_facet_id(elements[eid0], facet->left_face_num);

			facets.get(facet->sons[0])->parent = facets.get(facet->sons[1])->parent = lower_id;
			facets.get(facet->sons[3])->parent = facets.get(facet->sons[2])->parent = upper_id;

			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_HORZ;
			facet->sons[0] = upper_id;
			facet->sons[1] = lower_id;
			facet->sons[2] = facet->sons[3] = INVALID_IDX;
		}
		else if (face_refinement == H3D_REFT_QUAD_VERT) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, eid1, facet->left_face_num, INVALID_IDX, -1);
			Facet *lower_facet = add_quad_facet(Facet::INNER, eid0, facet->left_face_num, INVALID_IDX, -1);

			upper_facet->parent = fid;
			upper_facet->ractive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = facet->sons[2];
			upper_facet->sons[0] = facet->sons[1];
			Word_t upper_id = get_facet_id(elements[eid1], facet->left_face_num);

			lower_facet->parent = fid;
			lower_facet->ractive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = facet->sons[3];
			lower_facet->sons[0] = facet->sons[0];
			Word_t lower_id = get_facet_id(elements[eid0], facet->left_face_num);

			facets.get(facet->sons[1])->parent = facets.get(facet->sons[2])->parent = upper_id;
			facets.get(facet->sons[0])->parent = facets.get(facet->sons[3])->parent = lower_id;

			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_VERT;
			facet->sons[2] = lower_id;
			facet->sons[3] = upper_id;
			facet->sons[0] = facet->sons[1] = INVALID_IDX;
		}
		else
			EXIT("Trying to apply incompatible face refinement to element #%d.", parent_elem->id);
	}
	else if (facet->type == Facet::INNER && facet->right == parent_elem->id) {
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
			Word_t upper_id = get_facet_id(elements[eid1], facet->right_face_num);

			lower_facet->parent = fid;
			lower_facet->lactive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = facet->sons[0];
			lower_facet->sons[3] = facet->sons[1];
			Word_t lower_id = get_facet_id(elements[eid0], facet->right_face_num);

			facets.get(facet->sons[0])->parent = facets.get(facet->sons[1])->parent = lower_id;
			facets.get(facet->sons[3])->parent = facets.get(facet->sons[2])->parent = upper_id;

			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_HORZ;
			facet->sons[1] = upper_id;
			facet->sons[0] = lower_id;
			facet->sons[2] = facet->sons[3] = INVALID_IDX;
		}
		else if (face_refinement == H3D_REFT_QUAD_VERT) {
			Facet *upper_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid1, facet->right_face_num);
			Facet *lower_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, eid0, facet->right_face_num);

			upper_facet->parent = fid;
			upper_facet->lactive = false;
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = facet->sons[2];
			upper_facet->sons[0] = facet->sons[1];
			Word_t upper_id = get_facet_id(elements[eid1], facet->right_face_num);

			lower_facet->parent = fid;
			lower_facet->lactive = false;
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = facet->sons[3];
			lower_facet->sons[0] = facet->sons[0];
			Word_t lower_id = get_facet_id(elements[eid0], facet->right_face_num);

			facets.get(facet->sons[1])->parent = facets.get(facet->sons[2])->parent = upper_id;
			facets.get(facet->sons[0])->parent = facets.get(facet->sons[3])->parent = lower_id;

			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_VERT;
			facet->sons[2] = lower_id;
			facet->sons[3] = upper_id;
			facet->sons[0] = facet->sons[1] = INVALID_IDX;
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
		// WTF?
		EXIT("Refining facet that does not face with appropriate element/boundary");
	}

	return true;
}

bool Mesh::refine_quad_facet(Hex *parent_elem, int iface, unsigned int face_refinement, Word_t eid0, Word_t eid1, Word_t eid2, Word_t eid3) {
	_F_
	assert(face_refinement == H3D_REFT_QUAD_BOTH);

	Word_t fid = get_facet_id(parent_elem, iface);
	Facet *facet = facets.get(fid);
	assert(facet->mode == MODE_QUAD);

	//	if (is_compatible_quad_refinement(facet, face_refinement)) {
	if (facet->type == Facet::INNER && facet->left == parent_elem->id) {
		// refine to the left
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == H3D_REFT_QUAD_BOTH) {
			// no refinement OR the same type of the refinement on the right side
			// -> refinement on the left side is safe
			facet->lactive = false;
			facet->ref_mask = H3D_REFT_QUAD_BOTH;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			for (int i = 0; i < 4; i++) {
				Facet *child_facet = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);
				child_facet->parent = fid;
				facet->sons[i] = get_facet_id(elements[ei[i]], facet->left_face_num);
			}
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_HORZ) { // FIXME: ignoring the orientation
			facet->lactive = false;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);

			Facet *upper_facet = facets.get(facet->sons[1]);
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = get_facet_id(elements[ei[3]], facet->left_face_num);
			upper_facet->sons[3] = get_facet_id(elements[ei[2]], facet->left_face_num);
			child_facets[2]->parent = child_facets[3]->parent = get_facet_id(elements[upper_facet->right], upper_facet->right_face_num);

			Facet *lower_facet = facets.get(facet->sons[0]);
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = get_facet_id(elements[ei[0]], facet->left_face_num);
			lower_facet->sons[3] = get_facet_id(elements[ei[1]], facet->left_face_num);
			child_facets[0]->parent = child_facets[1]->parent = get_facet_id(elements[lower_facet->right], lower_facet->right_face_num);
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_VERT) { // FIXME: ignoring the orientation
			facet->lactive = false;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, ei[i], facet->left_face_num, INVALID_IDX, -1);

			Facet *upper_facet = facets.get(facet->sons[3]);
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = get_facet_id(elements[ei[2]], facet->left_face_num);
			upper_facet->sons[0] = get_facet_id(elements[ei[1]], facet->left_face_num);
			child_facets[2]->parent = child_facets[1]->parent = get_facet_id(elements[upper_facet->right], upper_facet->right_face_num);

			Facet *lower_facet = facets.get(facet->sons[2]);
			lower_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			lower_facet->sons[1] = get_facet_id(elements[ei[3]], facet->left_face_num);
			lower_facet->sons[0] = get_facet_id(elements[ei[0]], facet->left_face_num);
			child_facets[0]->parent = child_facets[3]->parent = get_facet_id(elements[lower_facet->right], lower_facet->right_face_num);
		}
	}
	else if (facet->type == Facet::INNER && facet->right == parent_elem->id) {
		// refine to the right (with element)
		if (facet->ref_mask == H3D_REFT_FACE_NONE || facet->ref_mask == H3D_REFT_QUAD_BOTH) {
			// no refinement OR the same type of the refinement on the left side
			// -> refinement on the right side is safe
			facet->ractive = false;
			facet->ref_mask = H3D_REFT_QUAD_BOTH;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			for (int i = 0; i < 4; i++) {
				Facet *child_facet = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);
				child_facet->parent = fid;
				facet->sons[i] = get_facet_id(elements[ei[i]], facet->right_face_num);
			}

		}
		else if (facet->ref_mask == H3D_REFT_QUAD_HORZ) { // FIXME: ignoring the orientation
			facet->ractive = false;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);

			Facet *upper_facet = facets.get(facet->sons[1]);
			upper_facet->ref_mask = H3D_REFT_QUAD_VERT;
			upper_facet->sons[2] = get_facet_id(elements[ei[3]], facet->right_face_num);
			upper_facet->sons[3] = get_facet_id(elements[ei[2]], facet->right_face_num);
			child_facets[2]->parent = child_facets[3]->parent = get_facet_id(elements[upper_facet->left], upper_facet->left_face_num);

			Facet *lower_facet = facets.get(facet->sons[0]);
			lower_facet->ref_mask = H3D_REFT_QUAD_VERT;
			lower_facet->sons[2] = get_facet_id(elements[ei[0]], facet->right_face_num);
			lower_facet->sons[3] = get_facet_id(elements[ei[1]], facet->right_face_num);
			child_facets[0]->parent = child_facets[1]->parent = get_facet_id(elements[lower_facet->left], lower_facet->left_face_num);
		}
		else if (facet->ref_mask == H3D_REFT_QUAD_VERT) { // FIXME: ignoring the orientation
			facet->ractive = false;

			Word_t ei[4] = { eid0, eid1, eid2, eid3 };
			Facet *child_facets[4];
			for (int i = 0; i < 4; i++)
				child_facets[i] = add_quad_facet(Facet::INNER, INVALID_IDX, -1, ei[i], facet->right_face_num);

			Facet *upper_facet = facets.get(facet->sons[3]);
			upper_facet->ref_mask = H3D_REFT_QUAD_HORZ;
			upper_facet->sons[1] = get_facet_id(elements[ei[2]], facet->right_face_num);
			upper_facet->sons[0] = get_facet_id(elements[ei[1]], facet->right_face_num);
			child_facets[2]->parent = child_facets[1]->parent = get_facet_id(elements[upper_facet->left], upper_facet->left_face_num);

			Facet *lower_facet = facets.get(facet->sons[2]);
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

		Word_t ei[4] = { eid0, eid1, eid2, eid3 };
		for (int i = 0; i < 4; i++) {
			Facet *child_facet = add_quad_facet(Facet::OUTER, ei[i], facet->left_face_num, facet->right, facet->right_face_num);
			child_facet->parent = fid;
			facet->sons[i] = get_facet_id(elements[ei[i]], facet->left_face_num);
		}
	}
	else {
		// WTF?
		EXIT("Refining facet that does not face with appropriate element/boundary");
	}

	return true;
}

Facet *Mesh::add_quad_facet(Facet::Type type, Word_t left_elem, int left_iface, Word_t right_elem, int right_iface) {
	_F_
	Word_t elem_id;
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

	Word_t fidx = get_facet_id(elements[elem_id], iface);
	Facet *facet = NULL;
	if (fidx != INVALID_IDX) {
		// update info on existing facet
		facet = facets.get(fidx);
		if (elem_id == left_elem) facet->set_left_info(left_elem, left_iface);
		else facet->set_right_info(right_elem, right_iface);
	}
	else {
		// create new facet
		facet = new Facet(MODE_QUAD);
		MEM_CHECK(facet);
		facet->type = type;
		facet->set_left_info(left_elem, left_iface);
		facet->set_right_info(right_elem, right_iface);
	}

	Word_t facet_idxs[Quad::NUM_VERTICES];
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
			Word_t vtx[Edge::NUM_VERTICES] = { facet_idxs[i % Quad::NUM_EDGES], facet_idxs[(i + 1) % Quad::NUM_EDGES] };

			Edge edge;
			edges.lookup(vtx + 0, Edge::NUM_VERTICES, edge);
			edge.bnd = 1;
			edges.set(vtx, Edge::NUM_VERTICES, edge);
		}
	}

	this->facets.set(facet_idxs + 0, Quad::NUM_VERTICES, facet);

	return facet;
}

void Mesh::refine_all_elements(int refinement) {
	_F_
	FOR_ALL_ACTIVE_ELEMENTS(idx, this) {
		refine_element(idx, refinement);
	}
}

void Mesh::refine_by_criterion(int(*criterion)(Element* e), int depth) {
	// TODO: implement me
}

void Mesh::unrefine_element(Word_t id) {
	// TODO: implement me
}

void Mesh::unrefine_all_elements() {
	// TODO: implement me
}

Word_t Mesh::create_midpoint(Word_t a, Word_t b) {
	_F_
	// get vertices
	Vertex *v1 = vertices.get(a);
	Vertex *v2 = vertices.get(b);
	// create middle point
	return add_vertex((v1->x + v2->x) / 2.0, (v1->y + v2->y) / 2.0, (v1->z + v2->z) / 2.0);
}

Word_t Mesh::get_midpoint(Word_t a, Word_t b) {
	_F_
	Word_t idx = peek_midpoint(a, b);
	if (idx == INVALID_IDX) {
		idx = create_midpoint(a, b);
		Word_t pt[] = { a, b };
		midpoints.set(pt, Edge::NUM_VERTICES, idx);
	}
	return idx;
}

Word_t Mesh::peek_midpoint(Word_t a, Word_t b) const {
	_F_
	Word_t pt[] = { a, b };
	Word_t idx = INVALID_IDX;
	midpoints.lookup(pt, Edge::NUM_VERTICES, idx);
	return idx;
}

void Mesh::set_midpoint(Word_t a, Word_t b, Word_t idx) {
	_F_
	Word_t pt[] = { a, b };
	midpoints.set(pt, Edge::NUM_VERTICES, idx);
}

Word_t Mesh::get_edge_id(Word_t a, Word_t b) const {
	_F_
	Word_t pt[] = { a, b };
	return edges.get_idx(pt + 0, Edge::NUM_VERTICES);
}

/// referencing edges
void Mesh::ref_edges(Element *e) {
	_F_
	assert(e != NULL);

	for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
		Word_t vtx[Edge::NUM_VERTICES];
		e->get_edge_vertices(iedge, vtx);

		Edge edge;
		if (edges.lookup(vtx + 0, Edge::NUM_VERTICES, edge)) {
			edge.ref++;
			edges.set(vtx, Edge::NUM_VERTICES, edge);
		}
		else {
			edge.ref = 1;
			edges.set(vtx, Edge::NUM_VERTICES, edge);
		}
	}
}

void Mesh::unref_edges(Element *e) {
	_F_
	assert(e != NULL);

	for (int iedge = 0; iedge < e->get_num_edges(); iedge++) {
		Word_t vtx[Edge::NUM_VERTICES];
		e->get_edge_vertices(iedge, vtx);

		Edge edge;
		if (edges.lookup(vtx + 0, Edge::NUM_VERTICES, edge)) {
			edge.ref--;
			edges.set(vtx, Edge::NUM_VERTICES, edge);
		}
		else assert(false); // Unreferencing non-existent edge
	}
}

Word_t Mesh::get_facing_facet(Word_t fid, Word_t elem_id) {
	_F_
	Facet *facet = facets[fid];

	if (facet != NULL) {
		if (elem_id == facet->left) {
			while (!facet->ractive && facet->parent != INVALID_IDX) {
				fid = facet->parent;
				facet = facets[fid];
			}
			return fid;
		}
		else if (elem_id == facet->right) {
			while (!facet->lactive && facet->parent != INVALID_IDX) {
				fid = facet->parent;
				facet = facets[fid];
			}
			return fid;
		}
		else
			return INVALID_IDX;
	}
	else
		return INVALID_IDX;
}

Word_t Mesh::get_facet_id(int nv, ...) const {
	_F_
	Word_t k[nv];

	va_list ap;
	va_start(ap, nv);
	for (int i = 0; i < nv; i++)
		k[i] = va_arg(ap, Word_t);
	va_end(ap);

	return facets.get_idx(k + 0, nv);
}

void Mesh::regularize() {
	_F_
	// FIXME: implements only 1-irregularity rule (quite dirty hack this is)
	// Assumes only XYZ refinements of elements (i.e. no anisotropic refinements, no incompatible refinements)

	// NOTE: this is very stupid implementation. It checks the parent facet of a half-active facet. If this parent facet is
	// active on the opposite side, we found a hanging node of 1. order. If it is inactive, we check super parent (parent of
	// this parent facet) the same way. If it is active, we found hanging node of a  2. order and we refine this super parent.
	// If it is inactive, hanging node of a higher order was found and we report an error.

	FOR_ALL_ACTIVE_ELEMENTS(eid, this) {
		Element *elem = elements[eid];
		for (int iface = 0; iface < elem->get_num_faces(); iface++) {
			Word_t fid = get_facet_id(elem, iface);
			Facet *facet = facets[fid];
			assert(facet != NULL);
			if (facet->lactive && !facet->ractive) {
				if (facet->parent != INVALID_IDX) {
					Facet *parent = facets.get(facet->parent);
					if (parent->ractive) {
						// OK: 1. order hanging node
					}
					else {
						if (parent->parent != INVALID_IDX) {
							Facet *super_parent = facets.get(parent->parent);
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
			else if (!facet->lactive && facet->ractive) {
				if (facet->parent != INVALID_IDX) {
					Facet *parent = facets.get(facet->parent);
					if (parent->lactive) {
						// OK: 1. order hanging node
					}
					else {
						if (parent->parent != INVALID_IDX) {
							Facet *super_parent = facets.get(parent->parent);
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

	FOR_ALL_ACTIVE_ELEMENTS(eid, this) {
		Element *e = elements[eid];

		int split = H3D_SPLIT_NONE;
		for (int iface = 0; iface < e->get_num_faces(); iface++) {
			Word_t fid = get_facet_id(e, iface);
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
		refine_element(eid, reft[split]);
	}

	refine_towards_boundary(marker, depth - 1);
}

void Mesh::check_elem_oris()
{
	_F_
	RefMap refmap(this);
	FOR_ALL_ACTIVE_ELEMENTS(eid, this) {
		Element *e = elements[eid];
		refmap.set_active_element(e);

		order3_t ord;
		if (e->get_mode() == MODE_HEXAHEDRON) ord = order3_t(1, 1, 1);
		else if (e->get_mode() == MODE_TETRAHEDRON) ord = order3_t(2);
		else warning(H3D_ERR_NOT_IMPLEMENTED);
		Quad3D *quad = get_quadrature(e->get_mode());
		int np = quad->get_num_points(ord);
		QuadPt3D *pt = quad->get_points(ord);

		double3x3 *m = refmap.get_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			if (det(m[i]) <= 0) error("Element #%ld has an incorrect orientation.", eid);
		}
		delete [] m;
	}
}
