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

#include "mesh3d.h"
#include <string.h>
#include "../../../hermes_common/error.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/callstack.h"
#include "../mesh.h"
#include "../refdomain.h"

// maximal row length in bytes (used for reading the mesh3d-file)
#define MAX_ROW_LEN               1024

// exception error codes
#define E_CANT_OPEN_FILE          -1
#define E_READ_ERROR              -2

// number of markers on mesd3d file
#define MARKERS                   1

H3DReader::H3DReader() {
  _F_
}

H3DReader::~H3DReader() {
  _F_
}

static bool read_num(char *row, int &vertex_count) {
  _F_
  int n;
  if (sscanf(row, "%d", &n) == 1) {
    vertex_count = n;
    return true;
  }
  else return false;
}

static int read_n_nums(char *row, int n, double values[]) {
  _F_
  int i = 0;
  char delims[] = " \t\n\r";
  char *token = strtok(row, delims);
  while (token != NULL && i < n) {
  double d;
  sscanf(token, "%lf", &d);
  values[i++] = d;

  token = strtok(NULL, delims);
}

  return i;
}

static int read_n_nums(char *row, int n, unsigned int values[]) {
  _F_
  int i = 0;
  char delims[] = " \t\n\r";
  char *token = strtok(row, delims);
  while (token != NULL && i < n) {
  int n;
  sscanf(token, "%d", &n);
  values[i++] = n;

  token = strtok(NULL, delims);
}

  return i;
}

static bool range_check(unsigned int max_index, unsigned int *vs, int num_vs) {
_F_
for (int i = 0; i < num_vs; i++) {
if (vs[i] <= 0 || vs[i] > max_index) return false;
}
return true;
}

bool H3DReader::load(const char *file_name, Mesh *mesh) {
  _F_
  assert(mesh != NULL);

  FILE *file = fopen(file_name, "r");
  if (file == NULL) error("Could not open the mesh file %s", file_name);

  try {
    line_nr = 0;
    enum EState {
      STATE_VERTICES_NUM,
      STATE_VERTICES,
      STATE_TETRAS_NUM,
      STATE_TETRAS,
      STATE_HEXES_NUM,
      STATE_HEXES,
      STATE_PRISMS_NUM,
      STATE_PRISMS,
      STATE_TRIS_NUM,
      STATE_TRIS,
      STATE_QUADS_NUM,
      STATE_QUADS,
      STATE_OK
    } state = STATE_VERTICES_NUM;

    int vertex_count = 0;
    int tetra_count = 0;
    int hex_count = 0;
    int prism_count = 0;
    int tri_count = 0;
    int quad_count = 0;
    double buffer[10];
    unsigned int vs[10];
    int n;

    int max_vertex_index = 0;

    char row[MAX_ROW_LEN] = { 0 };
    while (fgets(row, MAX_ROW_LEN, file) != NULL && state != STATE_OK) {
      line_nr++;
      if (row[0] == '#') continue; // comment

      // trim right
      int i = strlen(row) - 1;
      for (; i >= 0; i--) {
        if (row[i] != ' ' && row[i] != '\r' && row[i] != '\n' && row[i] != '\t') break;
      }
      row[i + 1] = '\0';
      if (strlen(row) <= 0) continue; // skip empty lines

      switch (state) {
        case STATE_VERTICES_NUM:
          if (read_num(row, vertex_count)) {
            state = STATE_VERTICES;

            if (vertex_count <= 0) throw E_READ_ERROR;
            max_vertex_index = vertex_count; //vertices are counted from 1 in mesh3d format
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_VERTICES:
          if (read_n_nums(row, Vertex::NUM_COORDS, buffer) == Vertex::NUM_COORDS) {
            mesh->add_vertex(buffer[0], buffer[1], buffer[2]);

            vertex_count--;
            if (vertex_count == 0)
            {
              state = STATE_TETRAS_NUM;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_TETRAS_NUM:
          if (read_num(row, tetra_count)) {
            state = STATE_TETRAS;
            if (tetra_count <= 0){
            state = STATE_HEXES_NUM;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_TETRAS:
          if ((n = read_n_nums(row, Tetra::NUM_VERTICES + 1, vs)) >= Tetra::NUM_VERTICES) {
            if (!range_check(max_vertex_index, vs, Tetra::NUM_VERTICES)) {
              fprintf(stderr, "Invalid vertex index found in the section defining tetras (line %d).\n", line_nr);
              throw E_READ_ERROR;
            }

            Tetra *tet = mesh->add_tetra(vs);
            if (n > Tetra::NUM_VERTICES) tet->marker = vs[Tetra::NUM_VERTICES];
            tetra_count--;
            if (tetra_count == 0){
              state = STATE_HEXES_NUM;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_HEXES_NUM:
          if (read_num(row, hex_count)) {
            state = STATE_HEXES;
            if (hex_count <= 0){
            state = STATE_PRISMS_NUM;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_HEXES:
          if ((n = read_n_nums(row, Hex::NUM_VERTICES + 1, vs)) >= Hex::NUM_VERTICES) {
            if (!range_check(max_vertex_index, vs, Hex::NUM_VERTICES)) {
              fprintf(stderr, "Invalid vertex index found in the section defining hexes (line %d).\n", line_nr);
              throw E_READ_ERROR;
            }

            Hex *hex = mesh->add_hex(vs);
            //Hex *hex;
            if (n > Hex::NUM_VERTICES) hex->marker = vs[Hex::NUM_VERTICES];
            hex_count--;
            if (hex_count == 0){
              state = STATE_PRISMS_NUM;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_PRISMS_NUM:
          if (read_num(row, prism_count)) {
            state = STATE_PRISMS;
            if (prism_count <= 0) state = STATE_TRIS_NUM;
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_PRISMS:
          if ((n = read_n_nums(row, Prism::NUM_VERTICES + 1, vs)) >= Prism::NUM_VERTICES) {
            if (!range_check(max_vertex_index, vs, Prism::NUM_VERTICES)) {
              fprintf(stderr, "Invalid vertex index found in the section defining prisms (line %d).\n", line_nr);
              throw E_READ_ERROR;
            }

            Prism *pri = mesh->add_prism(vs);
              if (n > Hex::NUM_VERTICES) pri->marker = vs[Prism::NUM_VERTICES];
              prism_count--;
              if (prism_count == 0) state = STATE_TRIS_NUM;
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_TRIS_NUM:
          if (read_num(row, tri_count)) {
            state = STATE_TRIS;
            if (tri_count <= 0) state = STATE_QUADS_NUM;
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_TRIS:
          if (read_n_nums(row, Tri::NUM_VERTICES + MARKERS, vs) == Tri::NUM_VERTICES + MARKERS) {
            if (!range_check(max_vertex_index, vs, Tri::NUM_VERTICES)) {
              fprintf(stderr, "Invalid vertex index found in the section defining tris (line %d).\n", line_nr);
              throw E_READ_ERROR;
            }

            unsigned int facet_idxs[Tri::NUM_VERTICES] = { vs[0], vs[1], vs[2] };
            mesh->add_tri_boundary(facet_idxs, vs[Tri::NUM_VERTICES]);
            tri_count--;
            if (tri_count == 0) state = STATE_QUADS_NUM;
          }
          else {
            fprintf(stderr, "Not enough information for tris. You probably forgot to define boundary condition (line %d).\n", line_nr);
            throw E_READ_ERROR;
          }
          break;

        case STATE_QUADS_NUM:
          if (read_num(row, quad_count)) {
            state = STATE_QUADS;
            if (quad_count <= 0){
              state = STATE_OK;
            }
          }
          else
            throw E_READ_ERROR;
          break;

        case STATE_QUADS:
          if (read_n_nums(row, Quad::NUM_VERTICES + MARKERS, vs) == Quad::NUM_VERTICES + MARKERS) {
            if (!range_check(max_vertex_index, vs, Quad::NUM_VERTICES)) {
              fprintf(stderr, "Invalid vertex index found in the section defining quads (line %d).\n", line_nr);
              throw E_READ_ERROR;
            }

            unsigned int facet_idxs[Quad::NUM_VERTICES] = { vs[0], vs[1], vs[2], vs[3] };
            mesh->add_quad_boundary(facet_idxs, vs[Quad::NUM_VERTICES]);
            quad_count--;
            if (quad_count == 0){
              state = STATE_OK;
            }
          }
          else {
            fprintf(stderr, "Not enough information for quads. You probably forgot to define boundary condition (line %d).", line_nr);
            throw E_READ_ERROR;
          }
          break;

        case STATE_OK:
          break;
      }
    }

#ifdef HERMES_COMMON_CHECK_BOUNDARY_CONDITIONS
    // check if all "outer" faces have defined boundary condition
    for (std::map<Facet::Key, Facet*>::const_iterator it = mesh->facets.begin(); it != mesh->facets.end(); it++) {
      Facet *facet = it->second;

      if(((unsigned) facet->left == INVALID_IDX) || ((unsigned) facet->right == INVALID_IDX)) {
        fprintf(stderr, "Not all outer faces have defined boundary condition (line %d).", line_nr);
        throw E_READ_ERROR;
      }
    }
#endif

    mesh->ugh();
  }
  catch (int e) {
    fclose(file);
    return false;
  }

  fclose(file);

  return true;
}

bool H3DReader::save(const char *file_name, Mesh *mesh) {
  _F_
  assert(mesh != NULL);

  FILE *file = fopen(file_name, "w");
  if (file == NULL) return false;

  // save vertices
  fprintf(file, "# vertices\n");
  fprintf(file, "%lu\n", (unsigned long int)mesh->vertices.size());
  for(std::map<unsigned int, Vertex*>::const_iterator it = mesh->vertices.begin(); it != mesh->vertices.end(); it++) {
    Vertex *v = it->second;
    fprintf(file, "%lf %lf %lf\n", v->x, v->y, v->z);
  }
  fprintf(file, "\n");

  // elements
  std::map<unsigned int, Element *> tet, hex, pri;
  for(std::map<unsigned int, Element*>::const_iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++) {
    Element *elem = it->second;
    if (elem->active) {
      switch (elem->get_mode()) {
        case HERMES_MODE_TET: tet[it->first] = elem; break;
        case HERMES_MODE_HEX: hex[it->first] = elem; break;
        case HERMES_MODE_PRISM: pri[it->first] = elem; break;
      }
    }
  }

  // save tetras
  fprintf(file, "# tetras\n");
  fprintf(file, "%lu\n", (unsigned long int)tet.size());
  for(std::map<unsigned int, Element*>::const_iterator it = tet.begin(); it != tet.end(); it++) {
    unsigned int vtcs[Tetra::NUM_VERTICES];
    it->second->get_vertices(vtcs);
    fprintf(file, "%u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
  }
  fprintf(file, "\n");

  // save hexes
  fprintf(file, "# hexes\n");
  fprintf(file, "%lu\n", (unsigned long int)hex.size());
  for(std::map<unsigned int, Element*>::const_iterator it = hex.begin(); it != hex.end(); it++) {
    unsigned int vtcs[Hex::NUM_VERTICES];
    it->second->get_vertices(vtcs);
    fprintf(file, "%u %u %u %u %u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
  }
  fprintf(file, "\n");

  // save prisms
  fprintf(file, "# prisms\n");
  fprintf(file, "%lu\n", (unsigned long int)pri.size());
  for(std::map<unsigned int, Element*>::const_iterator it = pri.begin(); it != pri.end(); it++) {
    unsigned int vtcs[Prism::NUM_VERTICES];
    it->second->get_vertices(vtcs);
    fprintf(file, "%u %u %u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5]);
  }
  fprintf(file, "\n");

  // boundaries
  std::map<unsigned int, Facet *> tri_facets, quad_facets;
  for(std::map<Facet::Key, Facet*>::iterator it = mesh->facets.begin(); it != mesh->facets.end(); it++) {
    Facet *facet = it->second;
    if(facet->type == Facet::OUTER && mesh->elements[facet->left]->active) {
      switch (facet->type) {
        case HERMES_MODE_TRIANGLE: 
          unsigned int ii;
          for(ii = 0; ; ii++)
            if(tri_facets[ii] == NULL)
              break;
          tri_facets[ii] = facet;
          break;

        case HERMES_MODE_QUAD: 
          unsigned int ij;
          for(ij = 0; ; ij++)
            if(quad_facets[ij] == NULL)
              break;
          quad_facets[ij] = facet;
          break;
      }
    }
  }

  // tris
  fprintf(file, "# tris\n");
  fprintf(file, "%lu\n", (unsigned long int)tri_facets.size());
  for(std::map<unsigned int, Facet*>::const_iterator it = tri_facets.begin(); it != tri_facets.end(); it++) {
    Facet *facet = it->second;
    Boundary *bnd = mesh->boundaries[facet->right];
    Element *elem = mesh->elements[facet->left];

    unsigned int vtcs[Tri::NUM_VERTICES];
    elem->get_face_vertices(facet->left_face_num, vtcs);

    fprintf(file, "%u %u %u     %d\n", vtcs[0], vtcs[1], vtcs[2], bnd->marker);
  }
  fprintf(file, "\n");

  // quads
  fprintf(file, "# quads\n");
  fprintf(file, "%lu\n", (unsigned long int)quad_facets.size());
  for(std::map<unsigned int, Facet*>::const_iterator it = quad_facets.begin(); it != quad_facets.end(); it++) {
    Facet *facet = it->second;
    Boundary *bnd = mesh->boundaries[facet->right];
    Element *elem = mesh->elements[facet->left];

    unsigned int vtcs[Quad::NUM_VERTICES];
    elem->get_face_vertices(facet->left_face_num, vtcs);

    fprintf(file, "%u %u %u %u     %d\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], bnd->marker);
  }
  fprintf(file, "\n");

  // the end
  fclose(file);

  return true;
}
