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
// vtkoutputengine.cc
//
//

#include "../h3dconfig.h"
#include "vtk.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../common.h"
#include "../traverse.h"

#include <stdio.h>
#include <errno.h>
#include <common/utils.h>
#include <common/callstack.h>
#include <common/error.h>

// size of the buffer that is used for copying files
#define BUFLEN							8192
#define FORMAT							"%.17g"

#define VTK_TRIANGLE					5
#define VTK_QUAD						9

#define VTK_TETRA						10
#define VTK_HEXAHEDRON					12
#define VTK_WEDGE						13
#define VTK_PYRAMID						14

#define EPS								10e-15

static int divs[] = { 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6 };

namespace Vtk {

//// OutputQuad ////////////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class OutputQuad : public Quad3D {
public:
	virtual QuadPt3D *get_points(const order3_t &order) {
		_F_
		if (!tables.exists(order.get_idx())) calculate_view_points(order);
		return tables[order.get_idx()];
	}

	virtual int get_num_points(const order3_t &order) {
		_F_
		if (!np.exists(order.get_idx())) calculate_view_points(order);
		return np[order.get_idx()];
	}

protected:
	virtual void calculate_view_points(order3_t order) = 0;
};

//// OutputQuadTetra ///////////////////////////////////////////////////////////////////////////////

/// Quadrature for visualizing the solution on tetrahedron
///
/// @ingroup visualization
class OutputQuadTetra : public OutputQuad {
public:
	OutputQuadTetra();
	virtual ~OutputQuadTetra();

protected:
	virtual void calculate_view_points(order3_t order);
};

OutputQuadTetra::OutputQuadTetra()
{
#ifdef WITH_TETRA
	mode = MODE_TETRAHEDRON;
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

OutputQuadTetra::~OutputQuadTetra()
{
	_F_
#ifdef WITH_HEX
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];
#endif
}

void OutputQuadTetra::calculate_view_points(order3_t order)
{
	_F_
#ifdef WITH_TETRA
	int o = order.get_idx();
	np[o] = Tetra::NUM_VERTICES;
	tables[o] = new QuadPt3D[np[o]];

	const Point3D *ref_vtcs = RefTetra::get_vertices();

	for (int i = 0; i < Tetra::NUM_VERTICES; i++) {
		tables[o][i].x = ref_vtcs[i].x;
		tables[o][i].y = ref_vtcs[i].y;
		tables[o][i].z = ref_vtcs[i].z;
		tables[o][i].w = 1.0;	// not used
	}
#endif
}

/// Quadrature for visualizing the solution on hexahedron
///
/// @ingroup visualization
class OutputQuadHex : public OutputQuad {
public:
	OutputQuadHex();
	virtual ~OutputQuadHex();

protected:
	virtual void calculate_view_points(order3_t order);
};

OutputQuadHex::OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];
#endif
}

void OutputQuadHex::calculate_view_points(order3_t order) {
	_F_
#ifdef WITH_HEX
	int o = order.get_idx();
	np[o] = (divs[order.x] + 1) * (divs[order.y] + 1) * (divs[order.z] + 1);

	tables[o] = new QuadPt3D[np[o]];
	double step_x, step_y, step_z;
	step_x = 2.0 / divs[order.x];
	step_y = 2.0 / divs[order.y];
	step_z = 2.0 / divs[order.z];

	int n = 0;
	for (int k = 0; k < divs[order.z] + 1; k++) {
		for (int l = 0; l < divs[order.y] + 1; l++) {
			for (int m = 0; m < divs[order.x] + 1; m++, n++) {
				assert(n < np[o]);
				tables[o][n].x = (step_x * m) - 1;
				tables[o][n].y = (step_y * l) - 1;
				tables[o][n].z = (step_z * k) - 1;
				tables[o][n].w = 1.0;   // not used
			}
		}
	}
#endif
}

//// OutputQuadPrism ///////////////////////////////////////////////////////////////////////////////

/// TODO: output quad for prisms

////

/// A class to linearize a higher-order functions to be able to visualize them
///
/// NOTE: Not really a linearizer, so far it just stores vertices, elements and values
/// It is in Vtk namespace so that it does not interfere with the rest of the code, but
/// when a true linearizer is implemented this class will move from this namespace.
class Linearizer {
public:
	Linearizer() { }
	virtual ~Linearizer();

	struct Cell {
		enum EType {
			Hex, Tetra, Prism, Quad, Tri
		};

		int n;			// number of vertices for this cell
		int *idx;		// vertex indices
		EType type;		// type of the cell
	};

	Array<Vertex *> &get_points() { return points; }
	Array<Cell *> &get_cells() { return cells; }
	Array<double> &get_cell_data() { return cell_data; }
	Array<double> &get_point_data(int idx = 0) { return pt_data[idx]; }

	int add_point(double x, double y, double z);
	/// Adding cells prepared by VtkOuputEngine (0-based indexing)
	int add_cell(Linearizer::Cell::EType type, int n, int *vtcs);
	void set_cell_data(int i, double v) { cell_data[i] = v; }
	void set_point_data(int i, double v) { pt_data[0][i] = v; }
	void set_point_data(int i, double v0, double v1, double v2) {
		pt_data[0][i] = v0;
		pt_data[1][i] = v1;
		pt_data[2][i] = v2;
	}

public:
	Map<double3, int> vertex_id;					// mapping: CEDKey => ced function index
	Array<Vertex *> points;
	Array<Cell *> cells;
	Array<double> cell_data;
	Array<double> pt_data[3];
};

Linearizer::~Linearizer()
{
	_F_
	for (Word_t i = points.first(); i != INVALID_IDX; i = points.next(i))
		delete points[i];
	for (Word_t i = cells.first(); i != INVALID_IDX; i = cells.next(i)) {
		delete [] cells[i]->idx;
		delete cells[i];
	}
}

int Linearizer::add_point(double x, double y, double z)
{
	_F_
	double3 key = { x, y, z };
	int idx;
	if (vertex_id.lookup(key, idx))
		return idx;
	else {
		idx = points.add(new Vertex(x, y, z));
		vertex_id.set(key, idx);
		return idx;
	}
}

int Linearizer::add_cell(Linearizer::Cell::EType type, int n, int *vtcs)
{
	_F_
	Cell *cell = new Cell;
	cell->type = type;
	cell->n = n;
	cell->idx = new int [n];
	for (int i = 0; i < n; i++)
		cell->idx[i] = vtcs[i];

	return cells.add(cell);
}

//// FileFormatter /////////////////////////////////////////////////////////////////////////////////

/// Produces a files in VTK format
///
/// This class taker a Linearizer class, reads the info stored in there and produces a VTK
/// legacy file.
class FileFormatter {
public:
	FileFormatter(Linearizer *l) {
		lin = l;
	}

	/// Write the file
	/// @param[in] file - output file to write to
	/// @param[in] name - name of the variable we are putting out
	void write(FILE *file, const char *name);

protected:
	Linearizer *lin;
};

void FileFormatter::write(FILE *file, const char *name)
{
	_F_
	fprintf(file, "# vtk DataFile Version 2.0\n");
	fprintf(file, "\n");
	fprintf(file, "ASCII\n");

	// dataset
	Array<Vertex *> &points = lin->get_points();
	Array<Linearizer::Cell *> &cells = lin->get_cells();
	int sz_cells = 0;
	for (Word_t i = cells.first(); i != INVALID_IDX; i = cells.next(i)) {
		Linearizer::Cell *cell = cells[i];
		switch (cell->type) {
			case Linearizer::Cell::Hex: sz_cells += Hex::NUM_VERTICES; break;
			case Linearizer::Cell::Tetra: sz_cells += Tetra::NUM_VERTICES; break;
			case Linearizer::Cell::Prism: sz_cells += Prism::NUM_VERTICES; break;
			case Linearizer::Cell::Quad: sz_cells += Quad::NUM_VERTICES; break;
			case Linearizer::Cell::Tri: sz_cells += Tri::NUM_VERTICES; break;
		}
		sz_cells++;
	}
	Array<double> &cell_data = lin->get_cell_data();
	Array<double> &pt_data0 = lin->get_point_data(0);
	Array<double> &pt_data1 = lin->get_point_data(1);
	Array<double> &pt_data2 = lin->get_point_data(2);

	fprintf(file, "\n");
	fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(file, "POINTS %ld %s\n", points.count(), "float");
	for (Word_t i = points.first(); i != INVALID_IDX; i = points.next(i)) {
		Vertex *v = points[i];
		fprintf(file, "%e %e %e\n", v->x, v->y, v->z);
	}

	fprintf(file, "\n");
	fprintf(file, "CELLS %ld %d\n", cells.count(), sz_cells);
	for (Word_t i = cells.first(); i != INVALID_IDX; i = cells.next(i)) {
		Linearizer::Cell *cell = cells[i];

		fprintf(file, "%d", cell->n);
		for (int j = 0; j < cell->n; j++)
			fprintf(file, " %d", cell->idx[j]);
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	fprintf(file, "CELL_TYPES %ld\n", cells.count());
	for (Word_t i = cells.first(); i != INVALID_IDX; i = cells.next(i)) {
		Linearizer::Cell *cell = cells[i];

		int vtk_type = 0;
		switch (cell->type) {
			case Linearizer::Cell::Hex: vtk_type = VTK_HEXAHEDRON; break;
			case Linearizer::Cell::Tetra: vtk_type = VTK_TETRA; break;
			case Linearizer::Cell::Prism: vtk_type = VTK_WEDGE; break;
			case Linearizer::Cell::Quad: vtk_type = VTK_QUAD; break;
			case Linearizer::Cell::Tri: vtk_type = VTK_TRIANGLE; break;
		}
		fprintf(file, "%d\n", vtk_type);
	}

	fprintf(file, "\n");
	if (pt_data0.count() > 0) {
		if (pt_data2.count() > 0) {
			// point data
			fprintf(file, "POINT_DATA %ld\n", pt_data0.count());
			fprintf(file, "VECTORS %s %s\n", name, "float");
			for (Word_t i = pt_data0.first(); i != INVALID_IDX; i = pt_data0.next(i))
				fprintf(file, "%e %e %e\n", pt_data0[i], pt_data1[i], pt_data2[i]);
		}
		else {
			// point data
			fprintf(file, "POINT_DATA %ld\n", pt_data0.count());
			fprintf(file, "SCALARS %s %s %d\n", name, "float", 1);
			fprintf(file, "LOOKUP_TABLE %s\n", "default");
			for (Word_t i = pt_data0.first(); i != INVALID_IDX; i = pt_data0.next(i))
				fprintf(file, "%e\n", pt_data0[i]);
		}
	}
	else if (cell_data.count()) {
		// cell data
		fprintf(file, "CELL_DATA %ld\n", cell_data.count());
		fprintf(file, "SCALARS %s %s %d\n", name, "float", 1);
		fprintf(file, "LOOKUP_TABLE %s\n", "default");
		for (Word_t i = cell_data.first(); i != INVALID_IDX; i = cell_data.next(i))
			fprintf(file, "%e\n", cell_data[i]);
	}
}

} // namespace

//
#ifdef WITH_TETRA
static Vtk::OutputQuadTetra output_quad_tetra;
#define OUTPUT_QUAD_TETRA		&output_quad_tetra
#else
#define OUTPUT_QUAD_TETRA		NULL
#endif

#ifdef WITH_HEX
static Vtk::OutputQuadHex output_quad_hex;
#define OUTPUT_QUAD_HEX			&output_quad_hex
#else
#define OUTPUT_QUAD_HEX			NULL
#endif

static Vtk::OutputQuad *output_quad[] = { OUTPUT_QUAD_TETRA, OUTPUT_QUAD_HEX, NULL };

VtkOutputEngine::VtkOutputEngine(FILE *file, int outprec)
{
	_F_
	this->out_file = file;
	this->out_prec = outprec;
}

VtkOutputEngine::~VtkOutputEngine()
{
	_F_
}

void VtkOutputEngine::out(MeshFunction *fn, const char *name, int item)
{
	_F_
	assert(fn->get_num_components() == 1 || fn->get_num_components() == 3);

	int comp[COMPONENTS];		// components to output
	int nc;						// number of components to output
	int b = 0;
	if (fn->get_num_components() == COMPONENTS) {
		int a = 0;
		if ((item & FN_COMPONENT_0) && (item & FN_COMPONENT_1) && (item & FN_COMPONENT_2)) {
			mask_to_comp_val(item, a, b);
			for (int i = 0; i < COMPONENTS; i++) comp[i] = i;
			nc = 3;
		}
		else if ((item & FN_COMPONENT_0) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_0, a, b);
			comp[0] = 0;
			nc = 1;
		}
		else if ((item & FN_COMPONENT_1) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_1, a, b);
			comp[0] = 1;
			nc = 1;
		}
		else if ((item & FN_COMPONENT_2) > 0) {
			mask_to_comp_val(item & FN_COMPONENT_2, a, b);
			comp[0] = 2;
			nc = 1;
		}
		else {
			fprintf(this->out_file, "Unable to satisfy your request\n");
			return;					// Do not know what user wants
		}
	}
	else if (fn->get_num_components() == 1) {
		mask_to_comp_val(item & FN_COMPONENT_0, comp[0], b);
		nc = 1;
	}
	else {
		fprintf(this->out_file, "Unable to satisfy your request\n");
		return;					// Do not know what user wants
	}

	Vtk::Linearizer l;
	Mesh *mesh = fn->get_mesh();
	// values
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		fn->set_active_element(element);

		int mode = element->get_mode();
		Vtk::OutputQuad *quad = output_quad[mode];
		order3_t order = fn->get_order();

		int np = quad->get_num_points(order);
		QuadPt3D *pt = quad->get_points(order);

		// get coordinates of all points
		RefMap *refmap = fn->get_refmap();
		double *x = refmap->get_phys_x(np, pt);
		double *y = refmap->get_phys_y(np, pt);
		double *z = refmap->get_phys_z(np, pt);

		int vtx_pt[np];		// indices of the vertices for current element
		for (int i = 0; i < np; i++)
			vtx_pt[i] = l.add_point(x[i], y[i], z[i]);

		int id;
		switch (mode) {
			case MODE_HEXAHEDRON:
				for (int i = 0; i < divs[order.z]; i++) {
					for (int j = 0; j < divs[order.y]; j++) {
						for (int o = 0; o < divs[order.x]; o++) {
							int cell[Hex::NUM_VERTICES];
							int base = ((divs[order.x] + 1) * (divs[order.y] + 1) * i) + ((divs[order.x] + 1) * j) + o;
							cell[0] = vtx_pt[base];
							cell[1] = vtx_pt[base + 1];
							cell[2] = vtx_pt[base + (divs[order.x] + 1) + 1];
							cell[3] = vtx_pt[base + (divs[order.x] + 1)];

							int pl = (divs[order.x] + 1) * (divs[order.y] + 1);
							cell[4] = vtx_pt[base + pl];
							cell[5] = vtx_pt[base + pl + 1];
							cell[6] = vtx_pt[base + pl + (divs[order.x] + 1) + 1];
							cell[7] = vtx_pt[base + pl + (divs[order.x] + 1)];
							id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, cell);
						}
					}
				}
				break;

			case MODE_TETRAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				break;

			case MODE_PRISM:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;

			default:
				EXIT(H3D_ERR_UNKNOWN_MODE);
				break;
		} // switch

		fn->precalculate(np, pt, item);
		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val[COMPONENTS];
		for (int ic = 0; ic < nc; ic++)
			val[ic] = fn->get_values(ic, FN);

		for (int i = 0; i < np; i++) {
#ifndef H3D_COMPLEX
			if (nc == 1) l.set_point_data(vtx_pt[i], val[0][i]);
			else l.set_point_data(vtx_pt[i], val[0][i], val[1][i], val[2][i]);
#else
			if (nc == 1) l.set_point_data(vtx_pt[i], REAL(val[0][i]));
			else l.set_point_data(vtx_pt[i], REAL(val[0][i]), REAL(val[1][i]), REAL(val[2][i]));
#endif
		}

		delete [] x;
		delete [] y;
		delete [] z;
	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out(MeshFunction *fn1, MeshFunction *fn2, MeshFunction *fn3, const char *name,
                          int item)
{
	_F_
	Vtk::Linearizer l;

	UniData **unidata = NULL;
	MeshFunction *fn[] = { fn1, fn2, fn3 };
	Mesh *mesh = NULL;
	Mesh *meshes[] = { fn1->get_mesh(), fn2->get_mesh(), fn3->get_mesh() };
	bool unimesh = false;	// FIXME: do not build union mesh if the meshes are the same
	if (meshes[0]->elements[1]->get_mode() == MODE_TETRAHEDRON)
		unimesh = false;

	// construct an union mesh
	if (unimesh) {
		Traverse trav;
		trav.begin(COMPONENTS, meshes);
		mesh = new Mesh;
		MEM_CHECK(mesh);
		unidata = trav.construct_union_mesh(mesh);
		trav.finish();

	}
	else
		mesh = meshes[0];

	// create a refmap on union mesh
	RefMap refmap;
	refmap.set_mesh(mesh);

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *e = mesh->elements[idx];
		// set active elements
		if (!unimesh) {
			for (int i = 0; i < COMPONENTS; i++)
				fn[i]->set_active_element(e);
		}
		else {
			for (int i = 0; i < COMPONENTS; i++) {
				fn[i]->set_active_element(unidata[i][e->id].e);
				fn[i]->set_transform(unidata[i][e->id].idx);
			}
		}

		int mode = e->get_mode();
		Vtk::OutputQuad *quad = output_quad[mode];
		order3_t order = max(fn1->get_order(), max(fn2->get_order(), fn3->get_order()));

		int np = quad->get_num_points(order);
		QuadPt3D *pt = quad->get_points(order);

		// get coordinates of all points
		refmap.set_active_element(e);
		double *x = refmap.get_phys_x(np, pt);
		double *y = refmap.get_phys_y(np, pt);
		double *z = refmap.get_phys_z(np, pt);

		int vtx_pt[np];		// indices of the vertices for current element
		for (int i = 0; i < np; i++)
			vtx_pt[i] = l.add_point(x[i], y[i], z[i]);

		int id;
		switch (mode) {
			case MODE_HEXAHEDRON:
				for (int i = 0; i < divs[order.z]; i++) {
					for (int j = 0; j < divs[order.y]; j++) {
						for (int o = 0; o < divs[order.x]; o++) {
							int cell[Hex::NUM_VERTICES];
							int base = ((divs[order.x] + 1) * (divs[order.y] + 1) * i) + ((divs[order.x] + 1) * j) + o;
							cell[0] = vtx_pt[base];
							cell[1] = vtx_pt[base + 1];
							cell[2] = vtx_pt[base + (divs[order.x] + 1) + 1];
							cell[3] = vtx_pt[base + (divs[order.x] + 1)];

							int pl = (divs[order.x] + 1) * (divs[order.y] + 1);
							cell[4] = vtx_pt[base + pl];
							cell[5] = vtx_pt[base + pl + 1];
							cell[6] = vtx_pt[base + pl + (divs[order.x] + 1) + 1];
							cell[7] = vtx_pt[base + pl + (divs[order.x] + 1)];
							id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, cell);
						}
					}
				}
				break;

			case MODE_TETRAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				break;

			case MODE_PRISM:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;

			default:
				EXIT(H3D_ERR_UNKNOWN_MODE);
				break;
		} // switch

		for (int i = 0; i < COMPONENTS; i++)
			fn[i]->precalculate(np, pt, item);

		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);
		scalar *val[COMPONENTS];
		for (int ic = 0; ic < COMPONENTS; ic++)
			val[ic] = fn[ic]->get_values(0, b);

		for (int i = 0; i < np; i++) {
#ifndef H3D_COMPLEX
			l.set_point_data(vtx_pt[i], val[0][i], val[1][i], val[2][i]);
#else
			l.set_point_data(vtx_pt[i], REAL(val[0][i]), REAL(val[1][i]), REAL(val[2][i]));
#endif
		}

		delete [] x;
		delete [] y;
		delete [] z;
	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);

	if (unimesh) delete mesh;
}

void VtkOutputEngine::out(Mesh *mesh)
{
	_F_
	Vtk::Linearizer l;
	// add cells
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		int nv = element->get_num_vertices();
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		int vtx_pt[nv];
		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			int idx = l.add_point(v->x, v->y, v->z);
			vtx_pt[i] = idx;
		}

		int id;
		switch (element->get_mode()) {
			case MODE_HEXAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, vtx_pt);
				break;

			case MODE_TETRAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				break;

			default:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;
		}
		l.set_cell_data(id, 0);
	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, "mesh");
}


void VtkOutputEngine::out_bc(Mesh *mesh, const char *name)
{
	_F_
	Vtk::Linearizer l;
	// add cells
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		for (int iface = 0; iface < element->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			int nv = element->get_num_face_vertices(iface);
			Word_t vtcs[nv];
			element->get_face_vertices(iface, vtcs);

			int vtx_pt[nv];		// indices of the vertices for current element
			for (int i = 0; i < nv; i++) {
				Vertex *v = mesh->vertices[vtcs[i]];
				vtx_pt[i] = l.add_point(v->x, v->y, v->z);
			}

			int id;
			switch (facet->mode) {
				case MODE_TRIANGLE:
					id = l.add_cell(Vtk::Linearizer::Cell::Tri, Tri::NUM_VERTICES, vtx_pt);
					break;
				case MODE_QUAD:
					id = l.add_cell(Vtk::Linearizer::Cell::Quad, Quad::NUM_VERTICES, vtx_pt);
					break;
				default:
					EXIT(H3D_ERR_NOT_IMPLEMENTED);
					break;
			}

			Boundary *bnd = mesh->boundaries[facet->right];
			l.set_cell_data(id, bnd->marker);
		}

	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out_orders(Space *space, const char *name)
{
	_F_
	Vtk::Linearizer l;
	Mesh *mesh = space->get_mesh();
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		order3_t ord = space->get_element_order(idx);
		Element *element = mesh->elements[idx];

		int nv = element->get_num_vertices();
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		int vtx_pt[nv];
		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			int idx = l.add_point(v->x, v->y, v->z);
			vtx_pt[i] = idx;
		}

		int id;
		switch (element->get_mode()) {
			case MODE_HEXAHEDRON:
				for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
					Word_t fvtcs[Quad::NUM_VERTICES];
					element->get_face_vertices(iface, fvtcs);

					Vertex *v[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]],
					                 mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
					int fctr = l.add_point((v[0]->x + v[2]->x) / 2.0,
					                       (v[0]->y + v[2]->y) / 2.0,
					                       (v[0]->z + v[2]->z) / 2.0);

					const int *fedges = element->get_face_edges(iface);
					for (int ie = 0; ie < Quad::NUM_EDGES; ie++) {
						int d;
						int iedge = fedges[ie];		// local edge number
						if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) d = ord.x;
						else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) d = ord.y;
						else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) d = ord.z;
						else assert(false);

						int cell_pts[3];		// triangle vertices
						const int *edge_pt_idx = element->get_edge_vertices(iedge);

						cell_pts[0] = vtx_pt[edge_pt_idx[0]];
						cell_pts[1] = vtx_pt[edge_pt_idx[1]];
						cell_pts[2] = fctr;

						id = l.add_cell(Vtk::Linearizer::Cell::Tri, Tri::NUM_VERTICES, cell_pts);
						l.set_cell_data(id, d);
					}
				}
				break;

			case MODE_TETRAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				l.set_cell_data(id, ord.order);
				break;

			default:
				error(H3D_ERR_NOT_IMPLEMENTED);
		}
	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out_elem_markers(Mesh *mesh, const char *name)
{
	_F_
	Vtk::Linearizer l;
	// add cells
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		int nv = element->get_num_vertices();
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		int vtx_pt[nv];		// indices of the vertices for current element
		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			vtx_pt[i] = l.add_point(v->x, v->y, v->z);
		}

		int id;
		switch (element->get_mode()) {
			case MODE_HEXAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, vtx_pt);
				break;

			case MODE_TETRAHEDRON:
				id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				break;

			default:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;
		}
		l.set_cell_data(id, element->marker);
	}

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out(Matrix *mat, bool structure)
{
	_F_
	fprintf(this->out_file, "# vtk DataFile Version 2.0\n");
	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "ASCII\n");

	fprintf(this->out_file, "\n");
	fprintf(this->out_file, "DATASET STRUCTURED_POINTS\n");

	int n = mat->get_size();
	fprintf(this->out_file, "DIMENSIONS %d %d 1\n", n, n);

	fprintf(this->out_file, "ASPECT_RATIO %d %d %d\n", 1, 1, 1);
	fprintf(this->out_file, "ORIGIN %lf %lf %lf", 0.0, 0.0, 0.0);
	fprintf(this->out_file, "POINT_DATA %d\n", n * n);
	fprintf(this->out_file, "SCALARS matrix double 1\n");

	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	SparseMatrix *m = dynamic_cast<SparseMatrix *>(mat);
	if (m == NULL) {
		// dense matrix
		warning(H3D_ERR_NOT_IMPLEMENTED);
	}
	else {
		if (m->row_storage) {
			for (int i = 0; i < n; i++) {
				int n_entries = m->get_num_row_entries(i);
				std::vector<double> vals(n_entries);
				std::vector<int> idxs(n_entries);
				int n_extracted = 0;
				m->extract_row_copy(i, n_entries, n_extracted, &vals[0], &idxs[0]);

				double val;
				for (int j = 0, k = 0; j < n; j++)
					if (idxs[k] == j) {
						if (structure) val = fabs(vals[k]) < EPS ? 1.0 : 0.0;
						else val = vals[k];

						fprintf(this->out_file, "%lf\n", val);
						k++;
					}
					else
						fprintf(this->out_file, "%lf\n", 0.0);
			}
		}
		else if (m->col_storage) {
			for (int i = 0; i < n; i++) {
				int n_entries = m->get_num_col_entries(i);
				std::vector<double> vals(n_entries);
				std::vector<int> idxs(n_entries);
				int n_extracted = 0;
				m->extract_col_copy(i, n_entries, n_extracted, &vals[0], &idxs[0]);

				double val;
				for (int j = 0, k = 0; j < n; j++)
					if (idxs[k] == j) {
						if (structure) val = fabs(vals[k]) < EPS ? 1.0 : 0.0;
						else val = vals[k];

						fprintf(this->out_file, "%lf\n", val);
						k++;
					}
					else
						fprintf(this->out_file, "%lf\n", 0.0);
			}
		}
		else {
			// extract element by element
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					double val;
					if (structure) val = ABS(m->get(i, j)) < EPS ? 1.0 : 0.0;
					else val = REAL(m->get(i, j));

					fprintf(this->out_file, "%lf\n", val);
				}
			}
		}
	}
}
