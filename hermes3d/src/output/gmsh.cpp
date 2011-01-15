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

#include "gmsh.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../h3d_common.h"
#include <stdio.h>
#include <errno.h>

// size of the buffer that is used for copying files
#define FORMAT							"%.17g"

///
#define AVGTV(a, b) { 0.5 * (tv[(a)].x + tv[(b)].x), 0.5 * (tv[(a)].y + tv[(b)].y), 0.5 * (tv[(a)].z + tv[(b)].z) }

namespace Gmsh {

//// OutputQuad //////////////////////////////////////////////////////////////////////////////

/// Common ancestor for output quadratures. Extends the interface of Quad3D
///
/// @ingroup visualization
class HERMES_API OutputQuad : public Quad3D {
public:
	virtual QuadPt3D *get_points(const Ord3 &order) {
		_F_
    if (tables->find(order.get_idx()) == tables->end()) 
      calculate_view_points(order);
		return (*tables)[order.get_idx()];
	}

	virtual int get_num_points(const Ord3 &order) {
		_F_
    if (np->find(order.get_idx()) == np->end()) 
      calculate_view_points(order);
		return (*np)[order.get_idx()];
	}

	virtual int *get_subdiv_modes(Ord3 order) {
		_F_
    if (subdiv_modes.find(order.get_idx()) == subdiv_modes.end()) 
      calculate_view_points(order);
		return subdiv_modes[order.get_idx()];
	}

	virtual int get_subdiv_num(Ord3 order) {
		_F_
    if (subdiv_num.find(order.get_idx()) == subdiv_num.end()) 
      calculate_view_points(order);
		return subdiv_num[order.get_idx()];
	}

	virtual QuadPt3D *get_face_points(int face, const Ord2 &order) {
		_F_
		EXIT(HERMES_ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	virtual void set_output_precision(int p) { output_precision = p; }

protected:
	std::map<unsigned int, int> subdiv_num;
	std::map<unsigned int, int *> subdiv_modes;
	int output_precision;

	virtual void calculate_view_points(Ord3 order) = 0;
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx) = 0;
};

//// OutputQuadTetra //////////////////////////////////////////////////////////////////////////////

/// Quadrature for visualizing the solution on tetrahedron
///
/// @ingroup visualization
class HERMES_API OutputQuadTetra : public OutputQuad {
public:
	OutputQuadTetra();
	virtual ~OutputQuadTetra();

protected:
	virtual void calculate_view_points(Ord3 order);
	virtual void recursive_division(const Point3D *ref_vtcs, QuadPt3D *table, int levels, int &idx);
};

OutputQuadTetra::OutputQuadTetra() {
	_F_
#ifdef WITH_TETRA
	mode = HERMES_MODE_TET;
	max_order = H3D_MAX_QUAD_ORDER;

	output_precision = 1;
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

OutputQuadTetra::~OutputQuadTetra() {
	_F_
#ifdef WITH_TETRA
	for(std::map<unsigned int, QuadPt3D*>::iterator it = tables->begin(); it != tables->end(); it++)
    delete [] it->second;

  for(std::map<unsigned int, int*>::iterator it = subdiv_modes.begin(); it != subdiv_modes.end(); it++)
    delete [] it->second;
#endif
}


void OutputQuadTetra::calculate_view_points(Ord3 order) {
	_F_
#ifdef WITH_TETRA
	int orderidx = order.get_idx();
	// check if the order is greater than 0, because we are taking log(o)
	if (order.order == 0) order.order++;

	// there should be refinement levels enough to catch the properties of order 'order' functions on an element
	// choose a different formula if this does not behave well
	int levels = int(log(double(order.order))) / log(double(2)) + 1;

	// each refinement level means that a tetrahedron is divided into 8 subtetrahedra
	// i.e., there are 8^levels resulting tetrahedra
	subdiv_num[orderidx] = (1 << (3 * levels));
	(*np)[orderidx] = subdiv_num[orderidx] * Tetra::NUM_VERTICES;

	// the new subelements are tetrahedra only
	subdiv_modes[orderidx] = new int[subdiv_num[orderidx]];
	MEM_CHECK(subdiv_modes[orderidx]);
	for (int i = 0; i < subdiv_num[orderidx]; i++)
		subdiv_modes[orderidx][i] = HERMES_MODE_TET;

	// compute the table of points recursively
	(*tables)[orderidx] = new QuadPt3D[(*np)[orderidx]];
	int idx = 0;
	const Point3D *ref_vtcs = RefTetra::get_vertices();
	recursive_division(ref_vtcs, (*tables)[orderidx], levels, idx);
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

int get_principal_order(Ord3 order) {
	assert(order.type == HERMES_MODE_HEX);
	return std::max(order.x, std::max(order.y, order.z));
}

/// Quadrature for visualizing the solution on hexahedron
///
/// @ingroup visualization
class HERMES_API OutputQuadHex : public OutputQuad {
public:
	OutputQuadHex();
	virtual ~OutputQuadHex();

protected:
	virtual void calculate_view_points(Ord3 order);
	virtual void recursive_division(const Point3D *tv, QuadPt3D *table, int levels, int &idx);
};

OutputQuadHex::OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	mode = HERMES_MODE_HEX;
	max_order = H3D_MAX_QUAD_ORDER;
	output_precision = 0;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	for(std::map<unsigned int, QuadPt3D*>::iterator it = tables->begin(); it != tables->end(); it++)
    delete [] it->second;

  for(std::map<unsigned int, int*>::iterator it = subdiv_modes.begin(); it != subdiv_modes.end(); it++)
    delete [] it->second;
#endif
}

void OutputQuadHex::calculate_view_points(Ord3 order) {
	_F_
#ifdef WITH_HEX
//	int o = get_principal_order(order);
//	int levels = int(log(o) / log(2)) + output_precision;
	int o = order.get_idx();
	int levels = 3;

	subdiv_num[o] = (1 << (3 * levels));
	(*np)[o] = subdiv_num[o] * Hex::NUM_VERTICES;

	subdiv_modes[o] = new int[subdiv_num[o]];
	MEM_CHECK(subdiv_modes[o]);
	// the new subelements are hexahedra only
	for (int i = 0; i < subdiv_num[o]; i++)
		subdiv_modes[o][i] = HERMES_MODE_HEX;

	// compute the table of points recursively
	(*tables)[o] = new QuadPt3D[(*np)[o]];
	int idx = 0;
	const Point3D *ref_vtcs = RefHex::get_vertices();
	recursive_division(ref_vtcs, (*tables)[o], levels, idx);
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
		case HERMES_MODE_TET: id = "SS"; break;
		case HERMES_MODE_HEX: id = "SH"; break;
		case HERMES_MODE_PRISM: EXIT("Unsupported mode."); break;
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
		case HERMES_MODE_TET: id = "VS"; break;
		case HERMES_MODE_HEX: id = "VH"; break;
		case HERMES_MODE_PRISM: EXIT("Unsupported mode."); break;
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

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];
		  int mode = element->get_mode();
		  // FIXME: get order from the space
		  Ord3 order;
		  switch (mode) {
			  case HERMES_MODE_TET: order = Ord3(1); break;
			  case HERMES_MODE_HEX: order = Ord3(1, 1, 1); break;
			  case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
			  default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
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

		  scalar **val = new scalar *[nc];
		  for (int ic = 0; ic < nc; ic++)
			  val[ic] = fn->get_values(comp[ic], b);

		  int pt_idx = 0;
		  // iterate through sub-elements and output them
		  for (int i = 0; i < subdiv_num; i++) {
			  int np;
			  switch (mode) {
				  case HERMES_MODE_TET: np = Tetra::NUM_VERTICES; break;
				  case HERMES_MODE_HEX: np = Hex::NUM_VERTICES; break;
				  case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
				  default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
			  }

			  // small buffers to hold values for one sub-element
        Point3D *phys_pt = new Point3D[np * nc];
        double **v = new double*[nc];
			  for(int i = 0; i < nc; i++)
          v[i] = new double [np];

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
        delete [] phys_pt;
        delete [] v;
		  }
    
      delete [] val;
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

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];
		  int mode = element->get_mode();
		  // FIXME: get order from the space
		  Ord3 order;
		  switch (mode) {
			  case HERMES_MODE_TET: order = Ord3(1); break;
			  case HERMES_MODE_HEX: order = Ord3(1, 1, 1); break;
			  case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
			  default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
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
				  case HERMES_MODE_TET: np = Tetra::NUM_VERTICES; break;
				  case HERMES_MODE_HEX: np = Hex::NUM_VERTICES; break;
				  case HERMES_MODE_PRISM: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
				  default: EXIT(HERMES_ERR_UNKNOWN_MODE); break;
			  }

			  // small buffers to hold values for one sub-element
			  Point3D *phys_pt = new Point3D[np * COMPONENTS];
        double **v = new double*[COMPONENTS];
			  for(int i = 0; i < COMPONENTS; i++)
          v[i] = new double [np];

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
        delete [] phys_pt;
        delete [] v;
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
	fprintf(this->out_file, "%lu\n", (unsigned long int)mesh->vertices.size());
  for(std::map<unsigned int, Vertex*>::iterator it = mesh->vertices.begin(); it != mesh->vertices.end(); it++) {
    Vertex *v = mesh->vertices[it->first];
    fprintf(this->out_file, "%u %lf %lf %lf\n", it->first, v->x, v->y, v->z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	int n_edges = 0;
	int n_faces = 0;
	// elements
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%u\n", mesh->get_num_active_elements());
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];

		  n_edges += element->get_num_edges();
		  n_faces += element->get_num_faces();

		  unsigned int *vtcs = new unsigned int[element->get_num_vertices()];
		  element->get_vertices(vtcs);

		  switch (element->get_mode()) {
			  case HERMES_MODE_TET:
				  fprintf(this->out_file, "%u 4 0 %u %u %u %u\n",
					  element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
				  break;

			  case HERMES_MODE_HEX:
				  fprintf(this->out_file, "%u 5 0 %u %u %u %u %u %u %u %u\n",
					  element->id, vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
				  break;

			  case HERMES_MODE_PRISM:
				  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
				  break;

			  default:
				  EXIT(HERMES_ERR_UNKNOWN_MODE);
				  break;
		  }
      delete [] vtcs;
	  }
	fprintf(this->out_file, "$EndElements\n");

	// edges
	// TODO: do not include edges twice or more
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", n_edges);
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used) {
      Element *element = mesh->elements[it->first];
		  unsigned int vtcs[Edge::NUM_VERTICES];
		  for (int iedge = 0; iedge < element->get_num_edges(); iedge++) {
			  element->get_edge_vertices(iedge, vtcs);
        unsigned int i = 0;
        std::map<Edge::Key, Edge*>::const_iterator it_inner = mesh->edges.begin();
        while(it_inner != mesh->edges.end() && it_inner->first != mesh->get_edge_id(vtcs[0], vtcs[1])) {
          it_inner++;
          i++;
        }
        fprintf(this->out_file, "%u 1 0 %u %u\n", i, vtcs[0], vtcs[1]);
		  }
	  }
	fprintf(this->out_file, "$EndElements\n");

	// faces
	// TODO: do not include faces twice
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", n_faces);
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used) {
      Element *element = mesh->elements[it->first];
		  for (int iface = 0; iface < element->get_num_faces(); iface++) {
			  int nv = element->get_num_face_vertices(iface);
			  unsigned int *vtcs = new unsigned int[nv];
			  element->get_face_vertices(iface, vtcs);
        unsigned int i = 0;
        std::map<Facet::Key, Facet*>::const_iterator it_inner = mesh->facets.begin();
        while(it_inner != mesh->facets.end() && it_inner->first != mesh->get_facet_id(element, iface)) {
          it_inner++;
          i++;
        }
			  switch (element->get_face_mode(iface)) {
				  case HERMES_MODE_TRIANGLE:
            fprintf(this->out_file, "%u 2 0 %u %u %u\n", i, vtcs[0], vtcs[1], vtcs[2]);
					  break;

				  case HERMES_MODE_QUAD:
					  fprintf(this->out_file, "%u 3 0 %u %u %u %u\n", i, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
					  break;
			  }
        delete [] vtcs;
      }
	  }
	fprintf(this->out_file, "$EndElements\n");
}

void GmshOutputEngine::out_bc_gmsh(Mesh *mesh, const char *name) {
	_F_
	// see Gmsh documentation on details (http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html)

	int fc = 0; 		// number of outer facets
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];
		  for (int iface = 0; iface < element->get_num_faces(); iface++) {
        Facet::Key fid = mesh->get_facet_id(element, iface);
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
	fprintf(this->out_file, "%lu\n", (unsigned long int)mesh->vertices.size());
	for(std::map<unsigned int, Vertex*>::iterator it = mesh->vertices.begin(); it != mesh->vertices.end(); it++) {
    Vertex *v = mesh->vertices[it->first];
    fprintf(this->out_file, "%u %lf %lf %lf\n", it->first, v->x, v->y, v->z);
	}
	fprintf(this->out_file, "$EndNodes\n");

	// elements
	fprintf(this->out_file, "$Elements\n");
	fprintf(this->out_file, "%d\n", fc);
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];

		  for (int iface = 0; iface < element->get_num_faces(); iface++) {
			  int nv = element->get_num_face_vertices(iface);
			  unsigned int *vtcs = new unsigned int[nv];
			  element->get_face_vertices(iface, vtcs);
        Facet::Key fid = mesh->get_facet_id(element, iface);
			  Facet *facet = mesh->facets[fid];
			  if (facet->type == Facet::INNER) continue;
        unsigned int i = 0;
        std::map<Facet::Key, Facet*>::const_iterator it_inner = mesh->facets.begin();
        while(it_inner != mesh->facets.end() && it_inner->first != mesh->get_facet_id(element, iface)) {
          it_inner++;
          i++;
        }
			  switch (facet->mode) {
				  case HERMES_MODE_TRIANGLE:
            fprintf(this->out_file, "%u 2 0 %u %u %u\n", i, vtcs[0], vtcs[1], vtcs[2]);
					  break;

				  case HERMES_MODE_QUAD:
					  fprintf(this->out_file, "%u 3 0 %u %u %u %u\n", i, vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
					  break;

				  default:
					  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
					  break;
			  }
        delete [] vtcs;
		  }
	  }
	fprintf(this->out_file, "$EndElements\n");

	// faces
	// TODO: do not include faces twice
	fprintf(this->out_file, "$ElementNodeData \n");
	fprintf(this->out_file, "1\n\"%s\"\n0\n3\n0\n1\n", name);
	fprintf(this->out_file, "%d\n", fc);
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
		  Element *element = mesh->elements[it->first];
		  for (int iface = 0; iface < element->get_num_faces(); iface++) {
        Facet::Key fid = mesh->get_facet_id(element, iface);
			  Facet *facet = mesh->facets[fid];
			  if (facet->type == Facet::INNER) continue;

			  Boundary *bnd = mesh->boundaries[facet->right];
			  int marker = bnd->marker;
        unsigned int i = 0;
        std::map<Facet::Key, Facet*>::const_iterator it_inner = mesh->facets.begin();
        while(it_inner != mesh->facets.end() && it_inner->first != mesh->get_facet_id(element, iface)) {
          it_inner++;
          i++;
        }
			  switch (facet->mode) {
				  case HERMES_MODE_TRIANGLE:
					  fprintf(this->out_file, "%u 3 %d %d %d\n", i, marker, marker, marker);
					  break;

				  case HERMES_MODE_QUAD:
					  fprintf(this->out_file, "%u 4 %d %d %d %d\n", i, marker, marker, marker, marker);
					  break;

				  default:
					  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
					  break;
			  }
		  }
	  }
	fprintf(this->out_file, "$EndElementNodeData\n");
}

void GmshOutputEngine::out_orders_gmsh(Space *space, const char *name) {
	_F_
	Mesh *mesh = space->get_mesh();

	// prepare
	fprintf(this->out_file, "$MeshFormat\n");
	fprintf(this->out_file, "%.1lf %d %d\n", 2.0, 0, (int) sizeof(double));
	fprintf(this->out_file, "$EndMeshFormat\n");

	// HEX specific
	std::map<unsigned int, Vertex *> out_vtcs;	// vertices
	std::map<unsigned int, int> vtx_pt;			// mapping from mesh vertex id to output vertex id
  
  
	std::map<PtsKey, unsigned int> face_pts;			// id of points on faces
  std::map<PtsKey, unsigned int> ctr_pts;			// id of points in the center

	// nodes
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];
		  int nv = Hex::NUM_VERTICES;
		  unsigned int *vtcs = new unsigned int[nv];
		  element->get_vertices(vtcs);

		  for (int i = 0; i < nv; i++) {
			  Vertex *v = mesh->vertices[vtcs[i]];
        unsigned int j;
        for(j = 0; ; j++)
          if(out_vtcs[j] == NULL)
            break;
        out_vtcs[j] = new Vertex(*v);
			  vtx_pt[vtcs[i]] = j;
		  }

		  for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			  unsigned int *fvtcs = new unsigned int[Quad::NUM_VERTICES];
			  element->get_face_vertices(iface, fvtcs);

			  unsigned int k[] = { fvtcs[0], fvtcs[1], fvtcs[2], fvtcs[3] };
        PtsKey key(k, Quad::NUM_VERTICES);
        if (face_pts.find(key) == face_pts.end()) {
				  // create new vertex
				  Vertex *v[4] = { mesh->vertices[fvtcs[0]], mesh->vertices[fvtcs[1]], mesh->vertices[fvtcs[2]], mesh->vertices[fvtcs[3]] };
				  Vertex *fcenter = new Vertex((v[0]->x + v[2]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[2]->z) / 2.0);
          unsigned int j;
          for(j = 0; ; j++)
            if(out_vtcs[j] == NULL)
              break;
          out_vtcs[j] = fcenter;
				  face_pts[key] = j;
			  }
        delete [] fvtcs;
		  }


		  unsigned int c[] = { vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7] };
      PtsKey key(c, Hex::NUM_VERTICES);
      if (ctr_pts.find(key) == ctr_pts.end()) {
			  // create new vertex
			  Vertex *v[4] = { mesh->vertices[vtcs[0]], mesh->vertices[vtcs[1]], mesh->vertices[vtcs[3]], mesh->vertices[vtcs[4]] };
			  Vertex *center = new Vertex((v[0]->x + v[1]->x) / 2.0, (v[0]->y + v[2]->y) / 2.0, (v[0]->z + v[3]->z) / 2.0);
			  unsigned int j;
        for(j = 0; ; j++)
          if(out_vtcs[j] == NULL)
            break;
        out_vtcs[j] = center;
			  ctr_pts[key] = j;
		  }
      delete [] vtcs;
	  }

	fprintf(this->out_file, "$Nodes\n");
	fprintf(this->out_file, "%lu\n", (unsigned long int)out_vtcs.size());
  for(std::map<unsigned int, Vertex*>::iterator it = out_vtcs.begin(); it != out_vtcs.end(); it++) {
    Vertex *v = it->second;
    fprintf(this->out_file, "%d %lf %lf %lf\n", it->first + 1, v->x, v->y, v->z);			// IDs for GMSH are indexed from 1
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
	fprintf(this->out_file, "%u\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];
		  unsigned int *vtcs = new unsigned int[element->get_num_vertices()];
		  element->get_vertices(vtcs);

		  for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
			  unsigned int fvtcs[2][Quad::NUM_VERTICES];
			  element->get_face_vertices(eface[iedge][0], fvtcs[0]);
			  element->get_face_vertices(eface[iedge][1], fvtcs[1]);

        PtsKey key0(fvtcs[0], Quad::NUM_VERTICES);
        PtsKey key1(fvtcs[1], Quad::NUM_VERTICES);

			  unsigned int evtcs[2];
			  element->get_edge_vertices(iedge, evtcs);
			  unsigned int v[4] = { vtx_pt[evtcs[0]] + 1, face_pts[key0] + 1, vtx_pt[evtcs[1]] + 1, face_pts[key1] + 1 };
			  fprintf(this->out_file, "%u 3 0 %u %u %u %u\n", id++, v[0], v[1], v[2], v[3]);
		  }
      delete [] vtcs;
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
	fprintf(this->out_file, "%u\n", mesh->get_num_active_elements() * Hex::NUM_EDGES);
	id = 1;
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      assert(mesh->elements[it->first]->get_mode() == HERMES_MODE_HEX);			// HEX-specific
		  // get order from the space
      Ord3 order = space->get_element_order(it->first);

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
	error(HERMES_ERR_NOT_IMPLEMENTED);
}
