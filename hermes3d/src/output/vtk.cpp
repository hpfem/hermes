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

#include "vtk.h"
#include "../refdomain.h"
#include "../quadstd.h"
#include "../h3d_common.h"
#include "../traverse.h"

#include <stdio.h>
#include <errno.h>
#include "../../../hermes_common/utils.h"
#include "../../../hermes_common/error.h"

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

protected:
	virtual void calculate_view_points(Ord3 order) = 0;
};

//// OutputQuadTetra ///////////////////////////////////////////////////////////////////////////////

/// Quadrature for visualizing the solution on tetrahedron
///
/// @ingroup visualization
class HERMES_API OutputQuadTetra : public OutputQuad {
public:
	OutputQuadTetra();
	virtual ~OutputQuadTetra();

protected:
	virtual void calculate_view_points(Ord3 order);
};

OutputQuadTetra::OutputQuadTetra()
{
#ifdef WITH_TETRA
	mode = HERMES_MODE_TET;
#else
	EXIT(H3D_ERR_TETRA_NOT_COMPILED);
#endif
}

OutputQuadTetra::~OutputQuadTetra()
{
	_F_
#ifdef WITH_HEX
	for(std::map<unsigned int, QuadPt3D*>::iterator it = tables->begin(); it != tables->end(); it++)
    delete [] it->second;
#endif
}

void OutputQuadTetra::calculate_view_points(Ord3 order)
{
	_F_
#ifdef WITH_TETRA
	int o = order.get_idx();
	(*np)[o] = Tetra::NUM_VERTICES;
  if((*tables)[o] != NULL)
    delete [] (*tables)[o];
	(*tables)[o] = new QuadPt3D[(*np)[o]];

	const Point3D *ref_vtcs = RefTetra::get_vertices();

	for (int i = 0; i < Tetra::NUM_VERTICES; i++) {
		(*tables)[o][i].x = ref_vtcs[i].x;
		(*tables)[o][i].y = ref_vtcs[i].y;
		(*tables)[o][i].z = ref_vtcs[i].z;
		(*tables)[o][i].w = 1.0;	// not used
	}
#endif
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
};

OutputQuadHex::OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	mode = HERMES_MODE_HEX;
#else
	EXIT(H3D_ERR_HEX_NOT_COMPILED);
#endif
}

OutputQuadHex::~OutputQuadHex() {
	_F_
#ifdef WITH_HEX
	for(std::map<unsigned int, QuadPt3D*>::iterator it = tables->begin(); it != tables->end(); it++)
    delete [] it->second;
#endif
}

void OutputQuadHex::calculate_view_points(Ord3 order) {
	_F_
#ifdef WITH_HEX
	int o = order.get_idx();
	(*np)[o] = (divs[order.x] + 1) * (divs[order.y] + 1) * (divs[order.z] + 1);

  if((*tables)[o] != NULL)
    delete [] (*tables)[o];
	(*tables)[o] = new QuadPt3D[(*np)[o]];
	double step_x, step_y, step_z;
	step_x = 2.0 / divs[order.x];
	step_y = 2.0 / divs[order.y];
	step_z = 2.0 / divs[order.z];

	int n = 0;
	for (int k = 0; k < divs[order.z] + 1; k++) {
		for (int l = 0; l < divs[order.y] + 1; l++) {
			for (int m = 0; m < divs[order.x] + 1; m++, n++) {
				assert(n < (*np)[o]);
				(*tables)[o][n].x = (step_x * m) - 1;
				(*tables)[o][n].y = (step_y * l) - 1;
				(*tables)[o][n].z = (step_z * k) - 1;
				(*tables)[o][n].w = 1.0;   // not used
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
class HERMES_API Linearizer {
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

	std::map<unsigned int, Vertex *> &get_points() { return points; }
	std::map<unsigned int, Cell *> &get_cells() { return cells; }
	std::map<unsigned int, double> &get_cell_data() { return cell_data; }
	std::map<unsigned int, double> &get_point_data(int idx = 0) { return pt_data[idx]; }

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
  struct VertexIdKey
  {
    double x, y, z;
    VertexIdKey()
    {
      x = y = z = 0.0;
    }
    VertexIdKey(double x, double y, double z)
    {
      this->x = x;
      this->y = y;
      this->z = z;
    };
    bool operator<(const VertexIdKey & other) const
    {
      if(this->x < other.x)
        return true;
      else if(this->x > other.x)
        return false;
      else
        if(this->y < other.y)
          return true;
        else if(this->y > other.y)
          return false;
        else
            if(this->z < other.z)
            return true;
          else
            return false;
    };
  };

	std::map<VertexIdKey, int> vertex_id;
	std::map<unsigned int, Vertex *> points;
	std::map<unsigned int, Cell *> cells;
	std::map<unsigned int, double> cell_data;
	std::map<unsigned int, double> pt_data[3];
};

Linearizer::~Linearizer()
{
	_F_
  for(std::map<unsigned int, Vertex*>::iterator it = points.begin(); it != points.end(); it++)
    delete it->second;
  for(std::map<unsigned int, Cell*>::iterator it = cells.begin(); it != cells.end(); it++) {
    delete [] it->second->idx;
    delete it->second;
	}
}

int Linearizer::add_point(double x, double y, double z)
{
	_F_
	VertexIdKey key( x, y, z);
  if (vertex_id.find(key) != vertex_id.end())
		return vertex_id[key];
	else {
    unsigned int i;
    for(i = 0; ; i++)
      if(points[i] == NULL)
        break;
    points[i] = new Vertex(x, y, z);
		vertex_id[key] = i;
		return i;
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
  unsigned int i;
  for(i = 0; ; i++)
    if(cells[i] == NULL)
      break;
  cells[i] = cell;
	return i;
}

//// FileFormatter /////////////////////////////////////////////////////////////////////////////////

/// Produces a file in VTK format
///
/// This class takes a Linearizer class, reads the info stored in there and produces a VTK
/// legacy file.
class HERMES_API FileFormatter {
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
	std::map<unsigned int, Vertex *> &points = lin->get_points();
	std::map<unsigned int, Linearizer::Cell *> &cells = lin->get_cells();
	int sz_cells = 0;
  for(std::map<unsigned int, Vtk::Linearizer::Cell*>::iterator it = cells.begin(); it != cells.end(); it++) {
    Linearizer::Cell *cell = it->second;
		switch (cell->type) {
			case Linearizer::Cell::Hex: sz_cells += Hex::NUM_VERTICES; break;
			case Linearizer::Cell::Tetra: sz_cells += Tetra::NUM_VERTICES; break;
			case Linearizer::Cell::Prism: sz_cells += Prism::NUM_VERTICES; break;
			case Linearizer::Cell::Quad: sz_cells += Quad::NUM_VERTICES; break;
			case Linearizer::Cell::Tri: sz_cells += Tri::NUM_VERTICES; break;
		}
		sz_cells++;
	}
	std::map<unsigned int, double> &cell_data = lin->get_cell_data();
	std::map<unsigned int, double> &pt_data0 = lin->get_point_data(0);
	std::map<unsigned int, double> &pt_data1 = lin->get_point_data(1);
	std::map<unsigned int, double> &pt_data2 = lin->get_point_data(2);

	fprintf(file, "\n");
	fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(file, "POINTS %lu %s\n", (unsigned long int)points.size(), "float");
  for(std::map<unsigned int, Vertex*>::iterator it = points.begin(); it != points.end(); it++) {
    Vertex *v = it->second;
		fprintf(file, "%e %e %e\n", v->x, v->y, v->z);
	}

	fprintf(file, "\n");
	fprintf(file, "CELLS %lu %d\n", (unsigned long int)cells.size(), sz_cells);
  for(std::map<unsigned int, Vtk::Linearizer::Cell*>::iterator it = cells.begin(); it != cells.end(); it++) {
    Linearizer::Cell *cell = it->second;

		fprintf(file, "%d", cell->n);
		for (int j = 0; j < cell->n; j++)
			fprintf(file, " %d", cell->idx[j]);
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	fprintf(file, "CELL_TYPES %lu\n", (unsigned long int)cells.size());
  for(std::map<unsigned int, Vtk::Linearizer::Cell*>::iterator it = cells.begin(); it != cells.end(); it++) {
    Linearizer::Cell *cell = it->second;

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
	if (pt_data0.size() > 0) {
		if (pt_data2.size() > 0) {
			// point data
			fprintf(file, "POINT_DATA %lu\n", (unsigned long int)pt_data0.size());
			fprintf(file, "VECTORS %s %s\n", name, "float");
      for(std::map<unsigned int, double>::iterator it = pt_data0.begin(); it != pt_data0.end(); it++)
        fprintf(file, "%e %e %e\n", it->second, it->second, it->second);
		}
		else {
			// point data
			fprintf(file, "POINT_DATA %lu\n", (unsigned long int)pt_data0.size());
			fprintf(file, "SCALARS %s %s %d\n", name, "float", 1);
			fprintf(file, "LOOKUP_TABLE %s\n", "default");
      for(std::map<unsigned int, double>::iterator it = pt_data0.begin(); it != pt_data0.end(); it++)
        fprintf(file, "%e\n", it->second);
		}
	}
	else if (cell_data.size()) {
		// cell data
		fprintf(file, "CELL_DATA %lu\n", (unsigned long int)cell_data.size());
		fprintf(file, "SCALARS %s %s %d\n", name, "float", 1);
		fprintf(file, "LOOKUP_TABLE %s\n", "default");
    for(std::map<unsigned int, double>::iterator it = cell_data.begin(); it != cell_data.end(); it++)
      fprintf(file, "%e\n", it->second);
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
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
		  Element *element = mesh->elements[it->first];
		  fn->set_active_element(element);

		  int mode = element->get_mode();
		  Vtk::OutputQuad *quad = output_quad[mode];
		  Ord3 order = fn->get_order();

		  int np = quad->get_num_points(order);
		  QuadPt3D *pt = quad->get_points(order);

		  // get coordinates of all points
		  RefMap *refmap = fn->get_refmap();
		  double *x = refmap->get_phys_x(np, pt);
		  double *y = refmap->get_phys_y(np, pt);
		  double *z = refmap->get_phys_z(np, pt);

		  int *vtx_pt = new int[np];		// indices of the vertices for current element
		  for (int i = 0; i < np; i++)
			  vtx_pt[i] = l.add_point(x[i], y[i], z[i]);

		  int id;
		  switch (mode) {
			  case HERMES_MODE_HEX:
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

			  case HERMES_MODE_TET:
				  id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				  break;

			  case HERMES_MODE_PRISM:
				  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
				  break;

			  default:
				  EXIT(HERMES_ERR_UNKNOWN_MODE);
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

      delete [] vtx_pt;
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
	if (meshes[0]->elements[1]->get_mode() == HERMES_MODE_TET)
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

	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = mesh->elements[it->first];
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
		  Ord3 order = max(fn1->get_order(), max(fn2->get_order(), fn3->get_order()));

		  int np = quad->get_num_points(order);
		  QuadPt3D *pt = quad->get_points(order);

		  // get coordinates of all points
		  refmap.set_active_element(e);
		  double *x = refmap.get_phys_x(np, pt);
		  double *y = refmap.get_phys_y(np, pt);
		  double *z = refmap.get_phys_z(np, pt);

		  int *vtx_pt = new int[np];		// indices of the vertices for current element
		  for (int i = 0; i < np; i++)
			  vtx_pt[i] = l.add_point(x[i], y[i], z[i]);

		  int id;
		  switch (mode) {
			  case HERMES_MODE_HEX:
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

			  case HERMES_MODE_TET:
				  id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				  break;

			  case HERMES_MODE_PRISM:
				  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
				  break;

			  default:
				  EXIT(HERMES_ERR_UNKNOWN_MODE);
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
    
      delete []vtx_pt;
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
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];

		  int nv = element->get_num_vertices();
		  unsigned int *vtcs = new unsigned int[nv];
		  element->get_vertices(vtcs);

		  int *vtx_pt = new int[nv];
		  for (int i = 0; i < nv; i++) {
			  Vertex *v = mesh->vertices[vtcs[i]];
			  int idx = l.add_point(v->x, v->y, v->z);
			  vtx_pt[i] = idx;
		  }

		  int id;
		  switch (element->get_mode()) {
			  case HERMES_MODE_HEX:
				  id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, vtx_pt);
				  break;

			  case HERMES_MODE_TET:
				  id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				  break;

			  default:
				  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
				  break;
		  }
      delete [] vtcs;
      delete [] vtx_pt;
		  l.set_cell_data(id, 0);
	  }

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, "mesh");
}


void VtkOutputEngine::out_bc_vtk(Mesh *mesh, const char *name)
{
	_F_
	Vtk::Linearizer l;
	// add cells
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];

		  for (int iface = 0; iface < element->get_num_faces(); iface++) {
			  Facet::Key fid = mesh->get_facet_id(element, iface);
			  Facet *facet = mesh->facets[fid];
			  if (facet->type == Facet::INNER) continue;

			  int nv = element->get_num_face_vertices(iface);
			  unsigned int *vtcs = new unsigned int[nv];
			  element->get_face_vertices(iface, vtcs);

			  int *vtx_pt = new int[nv];		// indices of the vertices for current element
			  for (int i = 0; i < nv; i++) {
				  Vertex *v = mesh->vertices[vtcs[i]];
				  vtx_pt[i] = l.add_point(v->x, v->y, v->z);
			  }

			  int id;
			  switch (facet->mode) {
				  case HERMES_MODE_TRIANGLE:
					  id = l.add_cell(Vtk::Linearizer::Cell::Tri, Tri::NUM_VERTICES, vtx_pt);
					  break;
				  case HERMES_MODE_QUAD:
					  id = l.add_cell(Vtk::Linearizer::Cell::Quad, Quad::NUM_VERTICES, vtx_pt);
					  break;
				  default:
					  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
					  break;
			  }
      
        delete [] vtcs;
        delete [] vtx_pt;
			  Boundary *bnd = mesh->boundaries[facet->right];
			  l.set_cell_data(id, bnd->marker);
		  }

	  }

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out_orders_vtk(Space *space, const char *name)
{
	_F_
	Vtk::Linearizer l;
	Mesh *mesh = space->get_mesh();
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
		  Ord3 ord = space->get_element_order(it->first);
		  Element *element = mesh->elements[it->first];

		  int nv = element->get_num_vertices();
		  unsigned int *vtcs = new unsigned int[nv];
		  element->get_vertices(vtcs);

		  int *vtx_pt = new int[nv];
		  for (int i = 0; i < nv; i++) {
			  Vertex *v = mesh->vertices[vtcs[i]];
			  int idx = l.add_point(v->x, v->y, v->z);
			  vtx_pt[i] = idx;
		  }

		  int id;
		  switch (element->get_mode()) {
			  case HERMES_MODE_HEX:
				  for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
					  unsigned int fvtcs[Quad::NUM_VERTICES];
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

			  case HERMES_MODE_TET:
				  id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				  l.set_cell_data(id, ord.order);
				  break;
    
			  default:
				  error(HERMES_ERR_NOT_IMPLEMENTED);
		  }
      delete [] vtcs;
      delete [] vtx_pt;
	  }

	Vtk::FileFormatter fmt(&l);
	fmt.write(out_file, name);
}

void VtkOutputEngine::out_elem_markers(Mesh *mesh, const char *name)
{
	_F_
	Vtk::Linearizer l;
	// add cells
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *element = mesh->elements[it->first];

		  int nv = element->get_num_vertices();
		  unsigned int *vtcs = new unsigned int[nv];
		  element->get_vertices(vtcs);

		  int *vtx_pt = new int[nv];		// indices of the vertices for current element
		  for (int i = 0; i < nv; i++) {
			  Vertex *v = mesh->vertices[vtcs[i]];
			  vtx_pt[i] = l.add_point(v->x, v->y, v->z);
		  }

		  int id;
		  switch (element->get_mode()) {
			  case HERMES_MODE_HEX:
				  id = l.add_cell(Vtk::Linearizer::Cell::Hex, Hex::NUM_VERTICES, vtx_pt);
				  break;

			  case HERMES_MODE_TET:
				  id = l.add_cell(Vtk::Linearizer::Cell::Tetra, Tetra::NUM_VERTICES, vtx_pt);
				  break;

			  default:
				  EXIT(HERMES_ERR_NOT_IMPLEMENTED);
				  break;
		  }
		  l.set_cell_data(id, element->marker);
      delete [] vtcs;
      delete [] vtx_pt;
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

	unsigned int n = mat->get_size();
	fprintf(this->out_file, "DIMENSIONS %d %d 1\n", n, n);

	fprintf(this->out_file, "ASPECT_RATIO %d %d %d\n", 1, 1, 1);
	fprintf(this->out_file, "ORIGIN %lf %lf %lf", 0.0, 0.0, 0.0);
	fprintf(this->out_file, "POINT_DATA %d\n", n * n);
	fprintf(this->out_file, "SCALARS matrix double 1\n");

	fprintf(this->out_file, "LOOKUP_TABLE %s\n", "default");

	SparseMatrix *m = dynamic_cast<SparseMatrix *>(mat);
	if (m == NULL) {
		// dense matrix
		warning(HERMES_ERR_NOT_IMPLEMENTED);
	}
	else {
		if (m->row_storage) {
			for (unsigned int i = 0; i < n; i++) {
				int n_entries = m->get_num_row_entries(i);
				std::vector<double> vals(n_entries);
				std::vector<unsigned int> idxs(n_entries);
				unsigned int n_extracted = 0;
				m->extract_row_copy(i, n_entries, n_extracted, &vals[0], &idxs[0]);

				double val;
				for (unsigned int j = 0, k = 0; j < n; j++)
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
			for (unsigned int i = 0; i < n; i++) {
				int n_entries = m->get_num_col_entries(i);
				std::vector<double> vals(n_entries);
				std::vector<unsigned int> idxs(n_entries);
				unsigned int n_extracted = 0;
				m->extract_col_copy(i, n_entries, n_extracted, &vals[0], &idxs[0]);

				double val;
				for (unsigned int j = 0, k = 0; j < n; j++)
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
			for (unsigned int i = 0; i < n; i++) {
				for (unsigned int j = 0; j < n; j++) {
					double val;
					if (structure) val = ABS(m->get(i, j)) < EPS ? 1.0 : 0.0;
					else val = REAL(m->get(i, j));

					fprintf(this->out_file, "%lf\n", val);
				}
			}
		}
	}
}


/// Functions facilitating output in the format displayable by e.g. Paraview.
// Space (polynomial orders) output.
void out_orders_vtk(Space *space, const char *name, int iter)
{
  char fname[1024];
  if(iter == -1)
    sprintf(fname, "%s.vtk", name);
  else
    sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out_orders_vtk(space, name);
    fclose(f);
  }
  else
    warning("Could not open file '%s' for writing.", fname);
};

// Solution output for one solution component.
void out_fn_vtk(MeshFunction *fn, const char *name, int iter)
{
  char fname[1024];
  if(iter == -1)
    sprintf(fname, "%s.vtk", name);
  else
    sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(fn, name);
    fclose(f);
  }
  else warning("Could not open file '%s' for writing.", fname);
};

// Solution output for three solution components.
void out_fn_vtk(MeshFunction *x, MeshFunction *y, MeshFunction *z, const char *name, int iter) 
{
  char fname[1024];
  if(iter == -1)
    sprintf(fname, "%s.vtk", name);
  else
    sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(x, y, z, name);
    fclose(f);
  }
  else warning("Could not open file '%s' for writing.", fname);
};

// Boundary conditions output.
void out_bc_vtk(Mesh *mesh, const char *name, int iter)
{
  char of_name[1024];
  FILE *ofile;
  if(iter == -1)
    sprintf(of_name, "%s.vtk", name);
  else
    sprintf(of_name, "iter-%s-%d.vtk", name, iter);
  ofile = fopen(of_name, "w");
  if (ofile != NULL) {
    VtkOutputEngine output(ofile);
    output.out_bc_vtk(mesh, name);
    fclose(ofile);
  }
  else warning("Can not open '%s' for writing.", of_name);
};

// Mesh output.
void out_mesh_vtk(Mesh *mesh, const char *name, int iter)
{
  char fname[1024];
  if(iter == -1)
    sprintf(fname, "%s.vtk", name);
  else
    sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(mesh);
    fclose(f);
  }
  else warning("Could not open file '%s' for writing.", fname);
};

