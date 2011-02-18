// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include <string.h>
#include "mesh.h"
#include <map>
#include "hash.h"
#include <iostream>
#include "h2d_reader.h"

extern unsigned g_mesh_seq;

H2DReader::H2DReader()
{
}

H2DReader::~H2DReader()
{
}

void H2DReader::load_str(const char* mesh_str, Mesh *mesh)
{
  std::istringstream s(mesh_str);
  this->load_stream(s, mesh, "");
}

//// load_nurbs ////////////////////////////////////////////////////////////////////////////////////

Nurbs* H2DReader::load_nurbs(Mesh *mesh, Python &p, int id, Node** en, int &p1, int &p2)
{
  int i;
  Nurbs* nurbs = new Nurbs;

  p.exec("curve_n = len(curve)");
  int curve_n = p.pull_int("curve_n");
  if (curve_n < 0 || (curve_n != 3 && curve_n != 5))
    error("Invalid curve #%d.", id);
  bool circle = (curve_n == 3);
  nurbs->arc = circle;

  // read the end point indices
  p.exec("p1 = curve[0]");
  p1 = p.pull_int("p1");
  p.exec("p2 = curve[1]");
  p2 = p.pull_int("p2");

  *en = mesh->peek_edge_node(p1, p2);
  if (*en == NULL)
    error("Curve #%d: edge %d-%d does not exist.", id, p1, p2);

  // degree of curved edge
  nurbs->degree = 2;
  if (!circle)
  {
    p.exec("degree = curve[2]");
    nurbs->degree = p.pull_int("degree");
  }

  // get the number of control points
  int inner = 1, outer;
  if (!circle)
  {
    p.exec("inner = curve[3]");
    p.exec("inner_n = len(inner)");
    inner = p.pull_int("inner_n");
  }
  nurbs->np = inner + 2;

  // edge endpoints are also control points, with weight 1.0
  nurbs->pt = new double3[nurbs->np];
  nurbs->pt[0][0] = mesh->nodes[p1].x;
  nurbs->pt[0][1] = mesh->nodes[p1].y;
  nurbs->pt[0][2] = 1.0;
  nurbs->pt[inner+1][0] = mesh->nodes[p2].x;
  nurbs->pt[inner+1][1] = mesh->nodes[p2].y;
  nurbs->pt[inner+1][2] = 1.0;

  if (!circle)
  {
    // read inner control points
    for (i = 0; i < inner; i++)
    {
      p.push_int("i", i);
      p.exec("pt0, pt1, pt2 = inner[i]");
      nurbs->pt[i+1][0] = p.pull_double("pt0");
      nurbs->pt[i+1][1] = p.pull_double("pt1");
      nurbs->pt[i+1][2] = p.pull_double("pt2");
    }
  }
  else
  {
    // read the arc angle
    p.exec("angle = curve[2]");
    nurbs->angle = p.pull_double("angle");
    double a = (180.0 - nurbs->angle) / 180.0 * M_PI;

    // generate one control point
    double x = 1.0 / tan(a * 0.5);
    nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
    nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
    nurbs->pt[1][2] = cos((M_PI - a) * 0.5);
  }

  // get the number of knot vector points
  inner = 0;
  if (!circle)
  {
    p.exec("inner = curve[4]");
    p.exec("inner_n = len(inner)");
    inner = p.pull_int("inner_n");
  }

  nurbs->nk = nurbs->degree + nurbs->np + 1;
  outer = nurbs->nk - inner;
  if ((outer & 1) == 1)
    error("Curve #%d: incorrect number of knot points.", id);

  // knot vector is completed by 0.0 on the left and by 1.0 on the right
  nurbs->kv = new double[nurbs->nk];
  for (i = 0; i < outer/2; i++)
    nurbs->kv[i] = 0.0;
  if (inner) {
    for (i = outer/2; i < inner + outer/2; i++) {
      p.push_int("i", i-outer/2);
      p.exec("val = inner[i]");
      nurbs->kv[i] = p.pull_double("val");
    }
  }
  for (i = outer/2 + inner; i < nurbs->nk; i++)
    nurbs->kv[i] = 1.0;

  nurbs->ref = 0;
  return nurbs;
}


bool H2DReader::load(const char *filename, Mesh *mesh)
{
  std::ifstream s(filename);
  return this->load_stream(s, mesh, filename);
}

std::string read_file(std::istream &is)
{
    std::ostringstream s;
    s << is.rdbuf();
    return s.str();
}

PyMODINIT_FUNC initpython_reader(void); /*proto*/

bool H2DReader::load_stream(std::istream &is, Mesh *mesh,
        const char *filename)
{
  int i, j, k, n;
  Node* en;
  bool debug = false;

  mesh->free();

  std::string mesh_str = read_file(is);

  Python p;
  initpython_reader();
  p.exec("from python_reader import read_hermes_format_str");
  p.push_str("s", mesh_str);
  p.exec("vertices, elements, boundaries, curves, refinements"
          " = read_hermes_format_str(s)");

  //// vertices ////////////////////////////////////////////////////////////////

  p.exec("n = len(vertices)");
  n = p.pull_int("n");
  if (n < 0) error("File %s: 'vertices' must be a list.", filename);
  if (n < 2) error("File %s: invalid number of vertices.", filename);

  // create a hash table large enough
  int size = HashTable::H2D_DEFAULT_HASH_SIZE;
  while (size < 8*n) size *= 2;
  mesh->init(size);

  // create top-level vertex nodes
  for (i = 0; i < n; i++)
  {
    Node* node = mesh->nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = HERMES_TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;
    p.push_int("i", i);
    p.exec("x, y = vertices[i]");
    node->x = p.pull_double("x");
    node->y = p.pull_double("y");
  }
  mesh->ntopvert = n;

  //// elements ////////////////////////////////////////////////////////////////

  p.exec("n = len(elements)");
  n = p.pull_int("n");
  if (n < 0) error("File %s: 'elements' must be a list.", filename);
  if (n < 1) error("File %s: no elements defined.", filename);

  // create elements
  mesh->nactive = 0;
  for (i = 0; i < n; i++)
  {
    // read and check vertex indices
    p.push_int("i", i);
    p.exec("nv = len(elements[i])");
    int nv = p.pull_int("nv");
    int idx[5];
    std::string el_marker;
    if (!nv) { mesh->elements.skip_slot(); continue; }
    if (nv < 4 || nv > 5)
      error("File %s: element #%d: wrong number of vertex indices.", filename, i);
    if (nv == 4) {
      p.exec("n1, n2, n3, b = elements[i]");
      idx[0] = p.pull_int("n1");
      idx[1] = p.pull_int("n2");
      idx[2] = p.pull_int("n3");
      p.exec("b_str = 1 if isinstance(b, str) else 0");
      if (p.pull_int("b_str"))
          el_marker = p.pull_str("b");
      else
          idx[3] = p.pull_int("b");
    } else {
      p.exec("n1, n2, n3, n4, b = elements[i]");
      idx[0] = p.pull_int("n1");
      idx[1] = p.pull_int("n2");
      idx[2] = p.pull_int("n3");
      idx[3] = p.pull_int("n4");
      p.exec("b_str = 1 if isinstance(b, str) else 0");
      if (p.pull_int("b_str"))
          el_marker = p.pull_str("b");
      else
          idx[4] = p.pull_int("b");
    }
    for (j = 0; j < nv-1; j++)
      if (idx[j] < 0 || idx[j] >= mesh->ntopvert)
        error("File %s: error creating element #%d: vertex #%d does not exist.", filename, i, idx[j]);

    Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]];
    int marker;

    // If we are dealing with a string as a marker.
    if (el_marker != "") {
      // This functions check if the user-supplied marker on this element has been
      // already used, and if not, inserts it in the appropriate structure.
      mesh->markers_conversion->insert_element_marker(mesh->markers_conversion->min_element_marker_unused, el_marker);
      marker = mesh->markers_conversion->get_internal_element_marker(el_marker);
    }
    else {
      if(nv == 4) {
        // If we have some string-labeled boundary markers.
        if(mesh->markers_conversion != NULL) {
          // We need to make sure that the internal markers do not collide.
          mesh->markers_conversion->check_element_marker(idx[3]);
          mesh->markers_conversion->insert_element_marker(idx[3], "");
        }
        marker = idx[3];
      }
      else {
        // If we have some string-labeled boundary markers.
        if(mesh->markers_conversion != NULL) {
          // We need to make sure that the internal markers do not collide.
          mesh->markers_conversion->check_element_marker(idx[4]);
          mesh->markers_conversion->insert_element_marker(idx[4], "");
        }
        marker = idx[4];
      }
    }

    if(nv == 4) {
        check_triangle(i, v0, v1, v2);
        mesh->create_triangle(marker, v0, v1, v2, NULL);
      }
      else {
        Node *v3 = &mesh->nodes[idx[3]];
        check_quad(i, v0, v1, v2, v3);
        mesh->create_quad(marker, v0, v1, v2, v3, NULL);
      }

    mesh->nactive++;
  }
  mesh->nbase = n;

  //// boundaries //////////////////////////////////////////////////////////////
  p.exec("have_boundaries = 1 if boundaries else 0");
  if (p.pull_int("have_boundaries"))
  {
    p.exec("n = len(boundaries)");
    n = p.pull_int("n");

    // read boundary data
    for (i = 0; i < n; i++)
    {
      int v1, v2, marker;
      p.push_int("i", i);
      p.exec("v1, v2, marker = boundaries[i]");
      v1 = p.pull_int("v1");
      v2 = p.pull_int("v2");

      en = mesh->peek_edge_node(v1, v2);
      if (en == NULL)
        error("File %s: boundary data #%d: edge %d-%d does not exist", filename, i, v1, v2);

      int marker_to_set;
      p.exec("marker_str = 1 if isinstance(marker, str) else 0");
      int marker_str = p.pull_int("marker_str");

      // If we are dealing with a string as a marker.
      if (marker_str) {
        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        std::string smarker = p.pull_str("marker");
        mesh->markers_conversion->insert_boundary_marker(mesh->markers_conversion->min_boundary_marker_unused, smarker);
        marker_to_set = mesh->markers_conversion->get_internal_boundary_marker(smarker);
      }
      else {
        marker = p.pull_int("marker");
        // If we have some string-labeled boundary markers.
        if(mesh->markers_conversion != NULL) {
          // We need to make sure that the internal markers do not collide.
          mesh->markers_conversion->check_boundary_marker(marker);
          mesh->markers_conversion->insert_boundary_marker(marker, "");
        }
        marker_to_set = marker;
      }

      en->marker = marker_to_set;

      // This is extremely important, as in DG, it is assumed that negative boundary markers are reserved
      // for the inner edges.
      if (marker_to_set > 0)
      {
        mesh->nodes[v1].bnd = 1;
        mesh->nodes[v2].bnd = 1;
        en->bnd = 1;
      }
    }
  }

  // check that all boundary edges have a marker assigned
  for_all_edge_nodes(en, mesh)
    if (en->ref < 2 && en->marker == 0) {
      warn("Boundary edge node does not have a boundary marker");
    }

  //// curves //////////////////////////////////////////////////////////////////
  p.exec("have_curves = 1 if curves else 0");
  if (p.pull_int("have_curves"))
  {
    p.exec("n = len(curves)");
    n = p.pull_int("n");
    if (n < 0) error("File %s: 'curves' must be a list.", filename);

    // load curved edges
    for (i = 0; i < n; i++)
    {
      // load the control points, knot vector, etc.
      Node* en;
      int p1, p2;
      p.push_int("i", i);
      p.exec("curve = curves[i]");
      Nurbs* nurbs = load_nurbs(mesh, p, i, &en, p1, p2);

      // assign the nurbs to the elements sharing the edge node
      for (k = 0; k < 2; k++)
      {
        Element* e = en->elem[k];
        if (e == NULL) continue;

        if (e->cm == NULL)
        {
          e->cm = new CurvMap;
          memset(e->cm, 0, sizeof(CurvMap));
          e->cm->toplevel = 1;
          e->cm->order = 4;
        }

        int idx = -1;
        for (unsigned j = 0; j < e->nvert; j++)
          if (e->en[j] == en) { idx = j; break; }
        assert(idx >= 0);

        if (e->vn[idx]->id == p1)
        {
          e->cm->nurbs[idx] = nurbs;
          nurbs->ref++;
        }
        else
        {
          Nurbs* nurbs_rev = mesh->reverse_nurbs(nurbs);
          e->cm->nurbs[idx] = nurbs_rev;
          nurbs_rev->ref++;
        }
      }
      if (!nurbs->ref) delete nurbs;
    }
  }

  // update refmap coeffs of curvilinear elements
  Element* e;
  for_all_elements(e, mesh)
    if (e->cm != NULL)
      e->cm->update_refmap_coeffs(e);

  //// refinements /////////////////////////////////////////////////////////////
  p.exec("have_refinements = 1 if refinements else 0");
  if (p.pull_int("have_refinements"))
  {
    p.exec("n = len(refinements)");
    n = p.pull_int("n");
    if (n < 0) error("File %s: 'refinements' must be a list.", filename);

    // perform initial refinements
    for (i = 0; i < n; i++)
    {
      int id, ref;
      p.push_int("i", i);
      p.exec("id, ref = refinements[i]");
      id = p.pull_int("id");
      ref = p.pull_int("ref");
      mesh->refine_element(id, ref);
    }
  }
  mesh->ninitial = mesh->elements.get_num_items();

  mesh->seq = g_mesh_seq++;


  return true;
}

//// save ////////////////////////////////////////////////////////////////////////////////////

void H2DReader::save_refinements(Mesh *mesh, FILE* f, Element* e, int id, bool& first)
{
  if (e->active) return;
  fprintf(f, first ? "refinements =\n{\n" : ",\n"); first = false;
  if (e->bsplit())
  {
    fprintf(f, "  { %d, 0 }", id);
    int sid = mesh->seq; mesh->seq += 4;
    for (int i = 0; i < 4; i++)
      save_refinements(mesh, f, e->sons[i], sid+i, first);
  }
  else if (e->hsplit())
  {
    fprintf(f, "  { %d, 1 }", id);
    int sid = mesh->seq; mesh->seq += 2;
    save_refinements(mesh, f, e->sons[0], sid, first);
    save_refinements(mesh, f, e->sons[1], sid+1, first);
  }
  else
  {
    fprintf(f, "  { %d, 2 }", id);
    int sid = mesh->seq; mesh->seq += 2;
    save_refinements(mesh, f, e->sons[2], sid, first);
    save_refinements(mesh, f, e->sons[3], sid+1, first);
  }
}


void H2DReader::save_nurbs(Mesh *mesh, FILE* f, int p1, int p2, Nurbs* nurbs)
{
  if (nurbs->arc)
  {
    fprintf(f, "  { %d, %d, %.16g }", p1, p2, nurbs->angle);
  }
  else
  {
    int inner = nurbs->np - 2;
    int outer = nurbs->nk - inner;
    fprintf(f, "  { %d, %d, %d, { ", p1, p2, nurbs->degree);
    for (int i = 1; i < nurbs->np-1; i++)
      fprintf(f, "{ %.16g, %.16g, %.16g }%s ",
                 nurbs->pt[i][0], nurbs->pt[i][1], nurbs->pt[i][2],
                 i < nurbs->np-2 ? "," : "");

    fprintf(f, "}, { ");
    int max = nurbs->nk - (nurbs->degree+1);
    for (int i = nurbs->degree+1; i < max; i++)
      fprintf(f, "%.16g%s", nurbs->kv[i], i < max-1 ? "," : "");
    fprintf(f, "} }");
  }
}


static bool is_twin_nurbs(Element* e, int i)
{
  // on internal edges, where there are two Nurbs', we only save one of them
  return e->cm->nurbs[i]->twin && e->en[i]->ref == 2;
}

bool H2DReader::save(const char* filename, Mesh *mesh)
{
  int i, mrk;
  Element* e;

  // open output file
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not create mesh file.");
  //fprintf(f, "# hermes2d saved mesh\n\n");

  // save vertices
  fprintf(f, "vertices =\n{\n");
  for (i = 0; i < mesh->ntopvert; i++)
    fprintf(f, "  { %.16g, %.16g }%s\n", mesh->nodes[i].x, mesh->nodes[i].y, (i < mesh->ntopvert-1 ? "," : ""));

  // save elements
  fprintf(f, "}\n\nelements =\n{");
  bool first = true;
  for (i = 0; i < mesh->get_num_base_elements(); i++)
  {
    const char* nl = first ? "\n" : ",\n";  first = false;
    e = mesh->get_element_fast(i);
    if (!e->used)
      fprintf(f, "%s  { }", nl);
    else if (e->is_triangle())
      fprintf(f, "%s  { %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->marker);
    else
      fprintf(f, "%s  { %d, %d, %d, %d, %d }", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->vn[3]->id, e->marker);
  }

  // save boundary markers
  fprintf(f, "\n}\n\nboundaries =\n{");
  first = true;
  for_all_base_elements(e, mesh)
    for (unsigned i = 0; i < e->nvert; i++)
      if ((mrk = mesh->get_base_edge_node(e, i)->marker)) {
        const char* nl = first ? "\n" : ",\n";  first = false;
        fprintf(f, "%s  { %d, %d, %d }", nl, e->vn[i]->id, e->vn[e->next_vert(i)]->id, mrk);
      }
  fprintf(f, "\n}\n\n");

  // save curved edges
  first = true;
  for_all_base_elements(e, mesh)
    if (e->is_curved())
      for (unsigned i = 0; i < e->nvert; i++)
        if (e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i)) {
          fprintf(f, first ? "curves =\n{\n" : ",\n");  first = false;
          save_nurbs(mesh, f, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i]);
        }
  if (!first) fprintf(f, "\n}\n\n");

  // save refinements
  unsigned temp = mesh->seq;
  mesh->seq = mesh->nbase;
  first = true;
  for_all_base_elements(e, mesh)
    save_refinements(mesh, f, e, e->id, first);
  if (!first) fprintf(f, "\n}\n\n");

  mesh->seq = temp;
  fclose(f);

  return true;
}
