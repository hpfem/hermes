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
// gmshoutputengine.cc
//
//

#include "../h3dconfig.h"
#include "gmsh.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../common.h"
#include <stdio.h>
#include <errno.h>
#include <common/callstack.h>

// size of the buffer that is used for copying files
#define FORMAT							"%.17g"

///
#define AVGTV(a, b) { 0.5 * (tv[(a)].x + tv[(b)].x), 0.5 * (tv[(a)].y + tv[(b)].y), 0.5 * (tv[(a)].z + tv[(b)].z) }

namespace Gmsh {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

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

	virtual int *get_subdiv_modes(order3_t order) {
		_F_
		if (!subdiv_modes.exists(order.get_idx())) calculate_view_points(order);
		return subdiv_modes[order.get_idx()];
	}

	virtual int get_subdiv_num(order3_t order) {
		_F_
		if (!subdiv_num.exists(order.get_idx())) calculate_view_points(order);
		return subdiv_num[order.get_idx()];
	}

	virtual QuadPt3D *get_face_points(int face, const order2_t &order) {
		_F_
		EXIT(H3D_ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	virtual void set_output_precision(int p) { output_precision = p; }

protected:
	Array<int> subdiv_num;
	Array<int *> subdiv_modes;
	int output_precision;

	virtual void calculate_view_points(order3_t order) = 0;
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx) = 0;
};

//// OutputQuadTetra //////////////////////////////////////////////////////////////////////////////

/// Quadrature for visualizing the solution on tetrahedron
///
/// @ingroup visualization
class OutputQuadTetra : public OutputQuad {
public:
	OutputQuadTetra();
	virtual ~OutputQuadTetra();

protected:
	virtual void calculate_view_points(order3_t order);
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx);
};

OutputQuadTetra::OutputQuadTetra() {
	_F_
#ifdef WITH_TETRA
	mode = MODE_TETRAHEDRON;
	max_order = H3D_MAX_QUAD_ORDER;

	output_precision = 1;
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

OutputQuadTetra::~OutputQuadTetra() {
	_F_
#ifdef WITH_TETRA
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];

	for (Word_t i = subdiv_modes.first(); i != INVALID_IDX; i = subdiv_modes.next(i))
		delete[] subdiv_modes[i];
#endif
}


void OutputQuadTetra::calculate_view_points(order3_t order) {
	_F_
#ifdef WITH_TETRA
	int orderidx = order.get_idx();
	// check if the order is greater than 0, because we are taking log(o)
	if (order.order == 0) order.order++;

	// there should be refinement levels enough to catch the properties of order 'order' functions on an element
	// choose a different formula if this does not behave well
	int levels = int(log(order.order) / log(2)) + 1;

	// each refinement level means that a tetrahedron is divided into 8 subtetrahedra
	// i.e., there are 8^levels resulting tetrahedra
	subdiv_num[orderidx] = (1 << (3 * levels));
	np[orderidx] = subdiv_num[orderidx] * Tetra::NUM_VERTICES;

	// the new subelements are tetrahedra only
	subdiv_modes[orderidx] = new int[subdiv_num[orderidx]];
	MEM_CHECK(subdiv_modes[orderidx]);
	for (int i = 0; i < subdiv_num[orderidx]; i++)
		subdiv_modes[orderidx][i] = MODE_TETRAHEDRON;

	// compute the table of points recursively
	tables[orderidx] = new QuadPt3D[np[orderidx]];
	MEM_CHECK(tables[orderidx]);
	int idx = 0;
	const Point3D *ref_vtcs = RefTetra::get_vertices();
	recursive_division(ref_vtcs, tables[orderidx], levels, idx);
#endif
}

void OutputQuadTetra::recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx) {
	_F_
#ifdef WITH_TETRA
	if (levels == 0) {
		// vertices
		for (int i = 0; i < Tetra::NUM_VERTICES; i++) {
			table[idx].x = tv[i].x;
			table[idx].y = tv[i].y;
			table[idx].z = tv[i].z;
			table[idx].w = 1.0;
			idx++;
		}
	}
	else {
		Point3D div_vtcs[8][Tetra::NUM_VERTICES] = {
			{ { tv[0].x, tv[0].y, tv[0].z }, AVGTV(0,1), AVGTV(0,2), AVGTV(0,3) },
			{ AVGTV(0,1), { tv[1].x, tv[1].y, tv[1].z }, AVGTV(1,2), AVGTV(1,3) },
			{ AVGTV(0,2), AVGTV(1,2), { tv[2].x, tv[2].y, tv[2].z }, AVGTV(2,3) },
			{ AVGTV(0,3), AVGTV(1,3), AVGTV(2,3), { tv[3].x, tv[3].y, tv[3].z } },
			{ AVGTV(0,2), AVGTV(0,1), AVGTV(1,2), AVGTV(1,3) },
			{ AVGTV(0,3), AVGTV(0,1), AVGTV(0,2), AVGTV(1,3) },
			{ AVGTV(0,3), AVGTV(0,2), AVGTV(2,3), AVGTV(1,3) },
			{ AVGTV(2,3), AVGTV(0,2), AVGTV(1,2), AVGTV(1,3) }
		};

		for (int i = 0; i < 8; i++)
			recursive_division(div_vtcs[i], table, levels - 1, idx);
	}
#endif
}


//// OutputQuadHex ////////////////////////////////////////////////////////////////////////////////

int get_principal_order(order3_t order) {
	assert(order.type == MODE_HEXAHEDRON);
	return std::max(order.x, std::max(order.y, order.z));
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
	virtual void recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx);
};

OutputQuadHex::OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;
	max_order = H3D_MAX_QUAD_ORDER;
	output_precision = 0;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	for (Word_t i = tables.first(); i != INVALID_IDX; i = tables.next(i))
		delete[] tables[i];

	for (Word_t i = subdiv_modes.first(); i != INVALID_IDX; i = subdiv_modes.next(i))
		delete[] subdiv_modes[i];
#endif
}

void OutputQuadHex::calculate_view_points(order3_t order) {
	_F_
#ifdef WITH_HEX
//	int o = get_principal_order(order);
//	int levels = int(log(o) / log(2)) + output_precision;
	int o = order.get_idx();
	int levels = 3;

	subdiv_num[o] = (1 << (3 * levels));
	np[o] = subdiv_num[o] * Hex::NUM_VERTICES;

	subdiv_modes[o] = new int[subdiv_num[o]];
	MEM_CHECK(subdiv_modes[o]);
	// the new subelements are hexahedra only
	for (int i = 0; i < subdiv_num[o]; i++)
		subdiv_modes[o][i] = MODE_HEXAHEDRON;

	// compute the table of points recursively
	tables[o] = new QuadPt3D[np[o]];
	MEM_CHECK(tables[o]);
	int idx = 0;
	const Point3D *ref_vtcs = RefHex::get_vertices();
	recursive_division(ref_vtcs, tables[o], levels, idx);
#endif
}

void OutputQuadHex::recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx) {
	_F_
#ifdef WITH_HEX
	if (levels == 0) {
		// vertices
		for (int i = 0; i < Hex::NUM_VERTICES; i++) {
			table[idx].x = tv[i].x;
			table[idx].y = tv[i].y;
			table[idx].z = tv[i].z;
			table[idx].w = 1.0;
			idx++;
		}
	}
	else {
		Point3D div_vtcs[8][Hex::NUM_VERTICES] = {
			{ { tv[0].x, tv[0].y, tv[0].z }, AVGTV(0,1), AVGTV(0,2), AVGTV(0,3), AVGTV(0,4), AVGTV(0,5), AVGTV(0,6), AVGTV(0,7) },
			{ AVGTV(0,1), { tv[1].x, tv[1].y, tv[1].z }, AVGTV(1,2), AVGTV(1,3), AVGTV(1,4), AVGTV(1,5), AVGTV(1,6), AVGTV(1,7) },
			{ AVGTV(0,2), AVGTV(1,2), { tv[2].x, tv[2].y, tv[2].z }, AVGTV(2,3), AVGTV(2,4), AVGTV(2,5), AVGTV(2,6), AVGTV(2,7) },
			{ AVGTV(0,3), AVGTV(1,3), AVGTV(2,3), { tv[3].x, tv[3].y, tv[3].z }, AVGTV(3,4), AVGTV(3,5), AVGTV(3,6), AVGTV(3,7) },
			{ AVGTV(0,4), AVGTV(1,4), AVGTV(2,4), AVGTV(3,4), { tv[4].x, tv[4].y, tv[4].z }, AVGTV(4,5), AVGTV(4,6), AVGTV(4,7) },
			{ AVGTV(0,5), AVGTV(1,5), AVGTV(2,5), AVGTV(3,5), AVGTV(4,5), { tv[5].x, tv[5].y, tv[5].z }, AVGTV(5,6), AVGTV(5,7) },
			{ AVGTV(0,6), AVGTV(1,6), AVGTV(2,6), AVGTV(3,6), AVGTV(4,6), AVGTV(5,6), { tv[6].x, tv[6].y, tv[6].z }, AVGTV(6,7) },
			{ AVGTV(0,7), AVGTV(1,7), AVGTV(2,7), AVGTV(3,7), AVGTV(4,7), AVGTV(5,7), AVGTV(6,7), { tv[7].x, tv[7].y, tv[7].z } }
		};

		for (int i = 0; i < 8; i++)
			recursive_division(div_vtcs[i], table, levels - 1, idx);
	}
#endif
}


//// OutputQuadPrism /////////////////////////////////////////////////////////////////////////////

/// TODO: output quad for prisms

} // namespace

//
#ifdef WITH_TETRA
static Gmsh::OutputQuadTetra output_quad_tetra;
#define OUTPUT_QUAD_TETRA		&output_quad_tetra
#else
#define OUTPUT_QUAD_TETRA		NULL
#endif

#ifdef WITH_HEX
static Gmsh::OutputQuadHex output_quad_hex;
#define OUTPUT_QUAD_HEX			&output_quad_hex
#else
#define OUTPUT_QUAD_HEX			NULL
#endif

static Gmsh::OutputQuad *output_quad[] = { OUTPUT_QUAD_TETRA, OUTPUT_QUAD_HEX, NULL };

///

GmshOutputEngine::GmshOutputEngine(FILE *file) {
	_F_
	this->out_file = file;
}

GmshOutputEngine::~GmshOutputEngine() {
	_F_
}

void GmshOutputEngine::dump_scalars(int mode, int num_pts, Point3D *pts, double *value) {
	_F_
	const char *id;
	switch (mode) {
		case MODE_TETRAHEDRON: id = "SS"; break;
		case MODE_HEXAHEDRON: id = "SH"; break;
		case MODE_PRISM: EXIT("Unsupported mode."); break;
		default: EXIT("Invalid mode."); break;
	}

	// write id
	fprintf(this->out_file, "\t%s(", id);

	// write points
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT ", " FORMAT ", " FORMAT "%s", pts[j].x, pts[j].y, pts[j].z, j == num_pts - 1 ? "" : ", ");
	fprintf(this->out_file, ") { ");
	// write values
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT "%s", value[j], j == num_pts - 1 ? "" : ", ");
	// end the row
	fprintf(this->out_file, " };\n");
}

void GmshOutputEngine::dump_vectors(int mode, int num_pts, Point3D *pts, double *v0, double *v1, double *v2) {
	_F_
	const char *id;
	switch (mode) {
		case MODE_TETRAHEDRON: id = "VS"; break;
		case MODE_HEXAHEDRON: id = "VH"; break;
		case MODE_PRISM: EXIT("Unsupported mode."); break;
		default: EXIT("Invalid mode."); break;
	}

	// write id
	fprintf(this->out_file, "\t%s(", id);

	// write points
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT ", " FORMAT ", " FORMAT "%s", pts[j].x, pts[j].y, pts[j].z, j == num_pts - 1 ? "" : ", ");
	fprintf(this->out_file, ") { ");
	// write values
	for (int j = 0; j < num_pts; j++)
		fprintf(this->out_file, FORMAT ", " FORMAT ", " FORMAT "%s", v0[j], v1[j], v2[j], j == num_pts - 1 ? "" : ", ");
	// end the row
	fprintf(this->out_file, " };\n");
}

void GmshOutputEngine::out(MeshFunction *fn, const char *name, int item/* = FN_VAL*/) {
	_F_
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

	Mesh *mesh = fn->get_mesh();

	// prepare
	fprintf(this->out_file, "View \"%s\" {\n", name);

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();
		// FIXME: get order from the space
		order3_t order;
		switch (mode) {
			case MODE_TETRAHEDRON: order = order3_t(1); break;
			case MODE_HEXAHEDRON: order = order3_t(1, 1, 1); break;
			case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
			default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
		}

		Gmsh::OutputQuad *quad = output_quad[mode];
		int subdiv_num = quad->get_subdiv_num(order);

		fn->set_active_element(element);
		int np2 = quad->get_num_points(order);
		QuadPt3D *pt2 = quad->get_points(order);

		RefMap *refmap = fn->get_refmap();
		double *phys_x = refmap->get_phys_x(np2, pt2);
		double *phys_y = refmap->get_phys_y(np2, pt2);
		double *phys_z = refmap->get_phys_z(np2, pt2);

		fn->precalculate(np2, pt2, item);

		scalar *val[nc];
		for (int ic = 0; ic < nc; ic++)
			val[ic] = fn->get_values(comp[ic], b);

		int pt_idx = 0;
		// iterate through sub-elements and output them
		for (int i = 0; i < subdiv_num; i++) {
			int np;
			switch (mode) {
				case MODE_TETRAHEDRON: np = Tetra::NUM_VERTICES; break;
				case MODE_HEXAHEDRON: np = Hex::NUM_VERTICES; break;
				case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
				default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
			}

			// small buffers to hold values for one sub-element
			Point3D phys_pt[np * nc];
			double v[nc][np];
			for (int j = 0; j < np; j++, pt_idx++) {
				// physical coordinates of sub-element
				phys_pt[j].x = phys_x[pt_idx];
				phys_pt[j].y = phys_y[pt_idx];
				phys_pt[j].z = phys_z[pt_idx];

				for (int ic = 0; ic < nc; ic++) {
#ifndef H3D_COMPLEX
					v[ic][j] = val[ic][pt_idx];
#else
					assert(fabs(IMAG(val[ic][pt_idx])) < 1e-12); // output only for real numbers
					v[ic][j] = REAL(val[ic][pt_idx]);
#endif
				}
			}

			switch (nc) {
				case 1:
					dump_scalars(mode, np, phys_pt, v[0]);
					break;

				case 3:
					dump_vectors(mode, np, phys_pt, v[0], v[1], v[2]);
					break;

				default: assert(false); break;
			}
		}

		delete [] phys_x;
		delete [] phys_y;
		delete [] phys_z;
	}

	// finalize
	fprintf(this->out_file, "};\n");
}

void GmshOutputEngine::out(MeshFunction *fn1, MeshFunction *fn2, MeshFunction *fn3, const char *name, int item) {
	MeshFunction *fn[] = { fn1, fn2, fn3 };
	Mesh *mesh = fn[0]->get_mesh();

	// prepare
	fprintf(this->out_file, "View \"%s\" {\n", name);

	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int mode = element->get_mode();
		// FIXME: get order from the space
		order3_t order;
		switch (mode) {
			case MODE_TETRAHEDRON: order = order3_t(1); break;
			case MODE_HEXAHEDRON: order = order3_t(1, 1, 1); break;
			case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
			default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
		}

		Gmsh::OutputQuad *quad = output_quad[mode];
		int subdiv_num = quad->get_subdiv_num(order);

		for (int i = 0; i < COMPONENTS; i++) {
			fn[i]->set_active_element(element);
		}

		int np2 = quad->get_num_points(order);
		QuadPt3D *pt2 = quad->get_points(order);

		RefMap *refmap = fn[0]->get_refmap();
		double *phys_x = refmap->get_phys_x(np2, pt2);
		double *phys_y = refmap->get_phys_y(np2, pt2);
		double *phys_z = refmap->get_phys_z(np2, pt2);

		for (int i = 0; i < COMPONENTS; i++)
			fn[i]->precalculate(np2, pt2, item);

		int a = 0, b = 0;
		mask_to_comp_val(item, a, b);

		scalar *val[COMPONENTS];
		for (int i = 0; i < COMPONENTS; i++) {
			val[i] = fn[i]->get_values(0, b);
		}

		int pt_idx = 0;
		// iterate through sub-elements and output them
		for (int i = 0; i < subdiv_num; i++) {
			int np;
			switch (mode) {
				case MODE_TETRAHEDRON: np = Tetra::NUM_VERTICES; break;
				case MODE_HEXAHEDRON: np = Hex::NUM_VERTICES; break;
				case MODE_PRISM: EXIT(H3D_ERR_NOT_IMPLEMENTED); break;
				default: EXIT(H3D_ERR_UNKNOWN_MODE); break;
			}

			// small buffers to hold values for one sub-element
			Point3D phys_pt[np * COMPONENTS];
			double v[COMPONENTS][np];
			for (int j = 0; j < np; j++, pt_idx++) {
				// physical coordinates of sub-element
				phys_pt[j].x = phys_x[pt_idx];
				phys_pt[j].y = phys_y[pt_idx];
				phys_pt[j].z = phys_z[pt_idx];

				for (int ic = 0; ic < COMPONENTS; ic++) {
#ifndef H3D_COMPLEX
					v[ic][j] = val[ic][pt_idx];
#else
					v[ic][j] = REAL(val[ic][pt_idx]);
#endif
				}
			}

			dump_vectors(mode, np, phys_pt, v[0], v[1], v[2]);
		}
	}

	// finalize
	fprintf(this->out_file, "};\n");
}

void GmshOutputEngine::out(Mesh *mesh) {
	_F_
	// see Gmsh documentation on details (http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html)

	// header
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, (int) sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// vertices
	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%ld\n", mesh->vertices.count());
	FOR_ALL_VERTICES(idx, mesh) {
		Vertex *v = mesh->vertices[idx];
		fprintf(this->out_file, "%ld %lf %lf %lf\n", idx, v->x, v->y, v->z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	int n_edges = 0;
	int n_faces = 0;
	// elements
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%ld\n", mesh->get_num_active_elements());
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		n_edges += element->get_num_edges();
		n_faces += element->get_num_faces();

		Word_t vtcs[element->get_num_vertices()];
		element->get_vertices(vtcs);

		switch (element->get_mode()) {
			case MODE_TETRAHEDRON:
				fprintf(this->out_file, "%ld 4 0 %ld %ld %ld %ld\n",
					element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
				break;

			case MODE_HEXAHEDRON:
				fprintf(this->out_file, "%ld 5 0 %ld %ld %ld %ld %ld %ld %ld %ld\n",
					element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
				break;

			case MODE_PRISM:
				EXIT(H3D_ERR_NOT_IMPLEMENTED);
				break;

			default:
				EXIT(H3D_ERR_UNKNOWN_MODE);
				break;
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// edges
	// TODO: do not include edges twice or more
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", n_edges);
	FOR_ALL_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[Edge::NUM_VERTICES];
		for (int iedge = 0; iedge < element->get_num_edges(); iedge++) {
			element->get_edge_vertices(iedge, vtcs);
			fprintf(this->out_file, "%ld 1 0 %ld %ld\n", mesh->get_edge_id(vtcs[0], vtcs[1]), vtcs[0], vtcs[1]);
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// faces
	// TODO: do not include faces twice
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", n_faces);
	FOR_ALL_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		for (int iface = 0; iface < element->get_num_faces(); iface++) {
			int nv = element->get_num_face_vertices(iface);
			Word_t vtcs[nv];
			element->get_face_vertices(iface, vtcs);
			switch (element->get_face_mode(iface)) {
				case MODE_TRIANGLE:
					fprintf(this->out_file, "%ld 2 0 %ld %ld %ld\n", mesh->get_facet_id(element, iface), vtcs[0], vtcs[1], vtcs[2]);
					break;

				case MODE_QUAD:
					fprintf(this->out_file, "%ld 3 0 %ld %ld %ld %ld\n", mesh->get_facet_id(element, iface), vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
					break;
			}
		}
	}
	fprintf(this->out_file, "$EndElements\n");
}

void GmshOutputEngine::out_bc(Mesh *mesh, const char *name) {
	_F_
	// see Gmsh documentation on details (http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html)

	int fc = 0; 		// number of outer facets
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		for (int iface = 0; iface < element->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::OUTER) fc++;
		}
	}

	// header
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, (int) sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// vertices
	// TODO: dump only vertices on the boundaries
	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%ld\n", mesh->vertices.count());
	FOR_ALL_VERTICES(idx, mesh) {
		Vertex *v = mesh->vertices[idx];
		fprintf(this->out_file, "%ld %lf %lf %lf\n", idx, v->x, v->y, v->z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	// elements
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", fc);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];

		for (int iface = 0; iface < element->get_num_faces(); iface++) {
			int nv = element->get_num_face_vertices(iface);
			Word_t vtcs[nv];
			element->get_face_vertices(iface, vtcs);
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			switch (facet->mode) {
				case MODE_TRIANGLE:
					fprintf(this->out_file, "%ld 2 0 %ld %ld %ld\n", mesh->get_facet_id(element, iface), vtcs[0], vtcs[1], vtcs[2]);
					break;

				case MODE_QUAD:
					fprintf(this->out_file, "%ld 3 0 %ld %ld %ld %ld\n", mesh->get_facet_id(element, iface), vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
					break;

				default:
					EXIT(H3D_ERR_NOT_IMPLEMENTED);
					break;
			}
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// faces
	// TODO: do not include faces twice
	fprintf(this->out_file, "$ElementNodeData \n");
	fprintf(this->out_file, "1\n\"%s\"\n0\n3\n0\n1\n", name);
	fprintf(this->out_file, "%d\n", fc);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		for (int iface = 0; iface < element->get_num_faces(); iface++) {
			Word_t fid = mesh->get_facet_id(element, iface);
			Facet *facet = mesh->facets[fid];
			if (facet->type == Facet::INNER) continue;

			Boundary *bnd = mesh->boundaries[facet->right];
			int marker = bnd->marker;
			switch (facet->mode) {
				case MODE_TRIANGLE:
					fprintf(this->out_file, "%ld 3 %d %d %d\n", mesh->get_facet_id(element, iface), marker, marker, marker);
					break;

				case MODE_QUAD:
					fprintf(this->out_file, "%ld 4 %d %d %d %d\n", mesh->get_facet_id(element, iface), marker, marker, marker, marker);
					break;

				default:
					EXIT(H3D_ERR_NOT_IMPLEMENTED);
					break;
			}
		}
	}
	fprintf(this->out_file, "$EndElementNodeData\n");
}

void GmshOutputEngine::out_orders(Space *space, const char *name) {
	_F_
	Mesh *mesh = space->get_mesh();

	// prepare
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, (int) sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// HEX specific
	Array<Vertex *> out_vtcs;	// vertices
	Array<int> vtx_pt;			// mapping from mesh vertex id to output vertex id
	MapHSOrd face_pts;			// id of points on faces
	MapHSOrd ctr_pts;			// id of points in the center

	// nodes
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		int nv = Hex::NUM_VERTICES;
		Word_t vtcs[nv];
		element->get_vertices(vtcs);

		for (int i = 0; i < nv; i++) {
			Vertex *v = mesh->vertices[vtcs[i]];
			Word_t idx = out_vtcs.add(new Vertex(*v));
			vtx_pt[vtcs[i]] = idx;
		}

		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			Word_t fvtcs[Quad::NUM_VERTICES];
			element->get_face_vertices(iface, fvtcs);

			Word_t k[] = { fvtcs[0], fvtcs[1], fvtcs[2], fvtcs[3] };
			Word_t idx = INVALID_IDX;
			if (!face_pts.lookup(k, Quad::NUM_VERTICES, idx)) {
				// create new vertex
				Vertex *v[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]], mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
				Vertex *fcenter = new Vertex((v[0]->x + v[2]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[2]->z) / 2.0);
				Word_t idx = out_vtcs.add(fcenter);
				face_pts.set(k, Quad::NUM_VERTICES, idx);
			}
		}

		Word_t c[] = { vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7] };
		Word_t idx = INVALID_IDX;
		if (!ctr_pts.lookup(c, Hex::NUM_VERTICES, idx)) {
			// create new vertex
			Vertex *v[4] = { mesh->vertices[vtcs[0]], mesh->vertices[vtcs[1]], mesh->vertices[vtcs[3]], mesh->vertices[vtcs[4]] };
			Vertex *center = new Vertex((v[0]->x + v[1]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[3]->z) / 2.0);
			Word_t idx = out_vtcs.add(center);
			ctr_pts.set(c, Hex::NUM_VERTICES, idx);
		}
	}

	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%ld\n", out_vtcs.count());
	for (Word_t i = out_vtcs.first(); i != INVALID_IDX; i = out_vtcs.next(i)) {
		Vertex *v = out_vtcs[i];
		fprintf(this->out_file, "%ld %lf %lf %lf\n", i + 1, v->x, v->y, v->z);			// IDs for GMSH are indexed from 1
		delete v;																		// we no longer need the vertex data
	}
	fprintf(this->out_file, "$EndNodes\n");

	// faces that edges coincide with
	int eface[][2] = {
		{ 2, 4 }, { 1, 4 }, { 3, 4 }, { 0, 4 },
		{ 2, 0 }, { 1, 2 }, { 3, 1 }, { 0, 3 },
		{ 2, 5 }, { 1, 5 }, { 3, 5 }, { 0, 5 }
	};
	int id = 1;
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%ld\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		Element *element = mesh->elements[idx];
		Word_t vtcs[element->get_num_vertices()];
		element->get_vertices(vtcs);

		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			Word_t fvtcs[2][Quad::NUM_VERTICES];
			element->get_face_vertices(eface[iedge][0], fvtcs[0]);
			element->get_face_vertices(eface[iedge][1], fvtcs[1]);

			Word_t fidx[2] = { INVALID_IDX, INVALID_IDX };
			face_pts.lookup(fvtcs[0], Quad::NUM_VERTICES, fidx[0]);
			face_pts.lookup(fvtcs[1], Quad::NUM_VERTICES, fidx[1]);

			Word_t evtcs[2];
			element->get_edge_vertices(iedge, evtcs);
			Word_t v[4] = { vtx_pt[evtcs[0]] + 1, fidx[0] + 1, vtx_pt[evtcs[1]] + 1, fidx[1] + 1 };
			fprintf(this->out_file, "%d 3 0 %ld %ld %ld %ld\n", id++, v[0], v[1], v[2], v[3]);
		}
	}
	fprintf(this->out_file, "$EndElements\n");

	// node data
	fprintf(this->out_file, "$ElementNodeData\n");
	fprintf(this->out_file, "%d\n", 1);
	fprintf(this->out_file, "\"%s\"\n", name);
	fprintf(this->out_file, "%d\n", 0);
	fprintf(this->out_file, "%d\n", 3);
	fprintf(this->out_file, "0\n"); // time step (not used, but has to be there)
	fprintf(this->out_file, "1\n"); // 1 value per node
	fprintf(this->out_file, "%ld\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	id = 1;
	FOR_ALL_ACTIVE_ELEMENTS(idx, mesh) {
		assert(mesh->elements[idx]->get_mode() == MODE_HEXAHEDRON);			// HEX-specific
		// get order from the space
		order3_t order = space->get_element_order(idx);

		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			int dord;
			if (iedge == 0 || iedge == 2 || iedge == 8 || iedge == 10) dord = order.x;
			else if (iedge == 1 || iedge == 3 || iedge == 9 || iedge == 11) dord = order.y;
			else if (iedge == 4 || iedge == 5 || iedge == 6 || iedge == 7) dord = order.z;
			else assert(false);
			fprintf(this->out_file, "%d 4 %d %d %d %d\n", id++, dord, dord, dord, dord);
		}

	}
	fprintf(this->out_file, "$EndElementNodeData\n");
}

void GmshOutputEngine::out(Matrix *mat)
{
	_F_
	error(H3D_ERR_NOT_IMPLEMENTED);
}
