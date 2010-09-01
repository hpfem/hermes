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

//
// mesh.cc
//
// Implemetation of mesh
//

#include "config.h"
#include "mesh.h"
#include "mesh3dreader.h"



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

// HEX ////////////////////////////////////////////////////////////////////////

Hex::Hex() {
}

Hex::Hex(uint v[]) {
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
}

Hex::Hex(uint v1, uint v2, uint v3, uint v4, uint v5, uint v6, uint v7, uint v8) {
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
	vtcs[4] = v5;
	vtcs[5] = v6;
	vtcs[6] = v7;
	vtcs[7] = v8;
}

// TETRA  /////////////////////////////////////////////////////////////////////

Tetra::Tetra() {
}

Tetra::Tetra(uint v[]) {
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
}

Tetra::Tetra(uint v1, uint v2, uint v3, uint v4) {
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
}

// PRISM //////////////////////////////////////////////////////////////////////

Prism::Prism() {
}

Prism::Prism(uint v[]) {
	for (int i = 0; i < NUM_VERTICES; i++)
		vtcs[i] = v[i];
}

Prism::Prism(uint v1, uint v2, uint v3, uint v4, uint v5, uint v6) {
	vtcs[0] = v1;
	vtcs[1] = v2;
	vtcs[2] = v3;
	vtcs[3] = v4;
	vtcs[4] = v5;
	vtcs[5] = v6;
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

// BoundaryTri ////////////////////////////////////////////////////////////////

BoundaryTri::BoundaryTri(uint vtcs[], int marker) :
	Boundary(marker)
{
	for (int i = 0; i < NUM_VERTICES; i++)
		this->vtcs[i] = vtcs[i];
}

BoundaryTri::BoundaryTri(const BoundaryTri &o) :
	Boundary(o)
{
}

BoundaryTri::~BoundaryTri() {
}

Boundary *BoundaryTri::copy() {
	return new BoundaryTri(*this);
}

// BoundaryQuad ///////////////////////////////////////////////////////////////

BoundaryQuad::BoundaryQuad(uint vtcs[], int marker) :
	Boundary(marker)
{
	for (int i = 0; i < NUM_VERTICES; i++)
		this->vtcs[i] = vtcs[i];
}

BoundaryQuad::BoundaryQuad(const BoundaryQuad &o) :
	Boundary(o)
{
}

BoundaryQuad::~BoundaryQuad() {
}

Boundary *BoundaryQuad::copy() {
	return new BoundaryQuad(*this);
}

// Mesh ///////////////////////////////////////////////////////////////////////

Mesh::Mesh() {
}

int Mesh::add_vertex(double x, double y, double z) {
	return vertices.add(new Vertex(x, y, z));
}

Tetra *Mesh::create_tetra(uint vtcs[]) {
	Tetra *tetra = new Tetra(vtcs);
	int id = elements.add(tetra);
	tetra->id = id;
	return tetra;
}

Tetra *Mesh::add_tetra(uint vtcs[]) {
	Tetra *tetra = create_tetra(vtcs);
	return tetra;
}

Hex *Mesh::create_hex(uint vtcs[]) {
	// build up the element
	Hex *hex = new Hex(vtcs);
	int id = elements.add(hex);
	hex->id = id;
	return hex;
}

Hex *Mesh::add_hex(uint vtcs[]) {
	Hex *hex = create_hex(vtcs);
	return hex;
}

Prism *Mesh::create_prism(uint vtcs[]) {
	Prism *prism = new Prism(vtcs);
	int id = elements.add(prism);
	prism->id = id;
	return prism;
}

Prism *Mesh::add_prism(uint vtcs[]) {
	Prism *prism = create_prism(vtcs);
	return prism;
}

Boundary *Mesh::add_tri_boundary(uint vtcs[], int marker) {
	Facet *facet = NULL;
//	if (facets.lookup(vtcs + 0, Tri::NUM_VERTICES, facet)) {
		Boundary *bdr = new BoundaryTri(vtcs, marker);
		int pos = boundaries.add(bdr);
		bdr->id = pos;
//	}
}

Boundary *Mesh::add_quad_boundary(uint vtcs[], int marker) {
	Facet *facet = NULL;
//	if (facets.lookup(vtcs + 0, Quad::NUM_VERTICES, facet)) {
		Boundary *bdr = new BoundaryQuad(vtcs, marker);
		int pos = boundaries.add(bdr);
		bdr->id = pos;
//	}
}

bool Mesh::load(const char *file_name, Mesh3DLoader *loader) {
	if (loader->load(file_name, this))
		return true;
	else
		return false;
}
