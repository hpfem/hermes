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

#include "mesh_reader_h2d_bson.h"

#ifdef WITH_BSON

#include "mesh.h"
#include "api2d.h"

using namespace std;

namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderH2DBSON::MeshReaderH2DBSON()
    {
    }

    MeshReaderH2DBSON::~MeshReaderH2DBSON()
    {
    }

    bool MeshReaderH2DBSON::load(const char *filename, MeshSharedPtr mesh)
    {
      throw Exceptions::MethodNotImplementedException("MeshReaderH2DBSON::load");

      if(!mesh)
        throw Exceptions::NullException(1);

      mesh->free();

      std::map<unsigned int, unsigned int> vertex_is;

      FILE *fpr;
      fpr = fopen(filename, "rb");

      // file size:
      fseek (fpr, 0, SEEK_END);
      int size = ftell(fpr);
      rewind(fpr);

      // allocate memory to contain the whole file:
      char *datar = (char*) malloc (sizeof(char)*size);
      fread(datar, size, 1, fpr);
      fclose(fpr);

      bson br;
      bson_init_finished_data(&br, datar, 0);

      bson_iterator it;
      bson sub;

      // Vertices //
      bson_find(&it, &br, "vertex-count");
      int vertices_count = bson_iterator_int(&it);

      // Initialize mesh.
      int mesh_size = HashTable::H2D_DEFAULT_HASH_SIZE;
      while (mesh_size < 8 * vertices_count)
        mesh_size *= 2;
      mesh->init(mesh_size);

      // Create top-level vertex nodes.
      double* vertex_xes = new double[vertices_count];
      double* vertex_yes = new double[vertices_count];

      // Xes
      bson_iterator it_coeffs;
      bson_find(&it_coeffs, &br, "vertex-x");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      int index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_xes[index_coeff++] = bson_iterator_double(&it);

      // Yes
      bson_find(&it_coeffs, &br, "vertex-y");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_yes[index_coeff++] = bson_iterator_double(&it);

      // Assign.
      for (int vertex_i = 0; vertex_i < vertices_count; vertex_i++)
      {
        Node* node = mesh->nodes.add();
        assert(node->id == vertex_i);
        node->ref = TOP_LEVEL_REF;
        node->type = HERMES_TYPE_VERTEX;
        node->bnd = 0;
        node->p1 = node->p2 = -1;
        node->next_hash = NULL;
        node->x = vertex_xes[vertex_i];
        node->y = vertex_yes[vertex_i];
      }
      mesh->ntopvert = vertices_count;

      delete [] vertex_xes;
      delete [] vertex_yes;


      // Elements //
      bson_find(&it, &br, "element-count");
      int element_count = bson_iterator_int(&it);
      mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

      int* vertex_0 = new int[element_count];
      int* vertex_1 = new int[element_count];
      int* vertex_2 = new int[element_count];
      int* vertex_3 = new int[element_count];
      const char** markers = new const char*[element_count];

      // Vertex_0s
      bson_find(&it_coeffs, &br, "element-1");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_0[index_coeff++] = bson_iterator_int(&it);
      // Vertex_1s
      bson_find(&it_coeffs, &br, "element-2");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_1[index_coeff++] = bson_iterator_int(&it);
      // Vertex_2s
      bson_find(&it_coeffs, &br, "element-3");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_2[index_coeff++] = bson_iterator_int(&it);
      // Vertex_3s
      bson_find(&it_coeffs, &br, "element-4");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_3[index_coeff++] = bson_iterator_int(&it);
      // Markers
      bson_find(&it_coeffs, &br, "element-marker");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        markers[index_coeff++] = bson_iterator_string(&it);

      Element* e;
      for (int element_i = 0; element_i < element_count; element_i++)
      {
        mesh->element_markers_conversion.insert_marker(markers[element_i]);
        if(vertex_3[element_i] != -1)
          e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(markers[element_i]).marker,
          &mesh->nodes[vertex_0[element_i]],
          &mesh->nodes[vertex_1[element_i]],
          &mesh->nodes[vertex_2[element_i]],
          &mesh->nodes[vertex_3[element_i]],
          NULL);
        else
          e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(markers[element_i]).marker,
          &mesh->nodes[vertex_0[element_i]],
          &mesh->nodes[vertex_1[element_i]],
          &mesh->nodes[vertex_2[element_i]],
          NULL);
      }

      delete [] vertex_0;
      delete [] vertex_1;
      delete [] vertex_2;
      delete [] vertex_3;
      delete [] markers;

      // Boundaries //
      bson_find(&it, &br, "edge-count");
      int edges_count = bson_iterator_int(&it);

      int* vertex_0_edge = new int[edges_count];
      int* vertex_1_edge = new int[edges_count];
      const char** edge_markers = new const char*[edges_count];

      // Vertex_0s
      bson_find(&it_coeffs, &br, "edge-1");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_0_edge[index_coeff++] = bson_iterator_int(&it);
      // Vertex_1s
      bson_find(&it_coeffs, &br, "edge-2");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        vertex_1_edge[index_coeff++] = bson_iterator_int(&it);
      // Markers
      bson_find(&it_coeffs, &br, "edge-marker");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        edge_markers[index_coeff++] = bson_iterator_string(&it);

      Node* en;
      for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
      {
        en = mesh->peek_edge_node(vertex_0_edge[edge_i], vertex_1_edge[edge_i]);
        if(en == NULL)
          throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist.", edge_i, vertex_0_edge[edge_i], vertex_1_edge[edge_i]);

        std::string edge_marker = edge_markers[edge_i];

        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker(edge_marker);

        en->marker = mesh->boundary_markers_conversion.get_internal_marker(edge_marker).marker;

        // This is extremely important, as in DG, it is assumed that negative boundary markers are reserved
        // for the inner edges.
        if(en->marker > 0)
        {
          mesh->nodes[vertex_0_edge[edge_i]].bnd = 1;
          mesh->nodes[vertex_1_edge[edge_i]].bnd = 1;
          en->bnd = 1;
        }
      }

      // check that all boundary edges have a marker assigned
      for_all_edge_nodes(en, mesh)
        if(en->ref < 2 && en->marker == 0)
          this->warn("Boundary edge node does not have a boundary marker.");

      delete [] vertex_0_edge;
      delete [] vertex_1_edge;
      delete [] edge_markers;

      // Curves //

      // Save arcs.
      // - count.
      bson_find(&it, &br, "arc-count");
      int arc_count = bson_iterator_int(&it);

      int* p1s = new int[arc_count];
      int* p2s = new int[arc_count];
      double* angles = new double[arc_count];

      // p1s
      bson_find(&it_coeffs, &br, "arcs-1");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        p1s[index_coeff++] = bson_iterator_int(&it);
      // p2s
      bson_find(&it_coeffs, &br, "arcs-2");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        p2s[index_coeff++] = bson_iterator_int(&it);
      // angles
      bson_find(&it_coeffs, &br, "arcs-angle");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        angles[index_coeff++] = bson_iterator_double(&it);

      // Arcs //
      for (unsigned int curves_i = 0; curves_i < arc_count; curves_i++)
      {
        // load the control points, knot vector, etc.
        Node* en;
        int p1, p2;
        double angle;

        Nurbs* nurbs;

        // read the end point indices
        p1 = p1s[curves_i];
        p2 = p2s[curves_i];
        angle = angles[curves_i];

        nurbs = load_arc(mesh, curves_i, &en, p1, p2, angle);

        // assign the arc to the elements sharing the edge node
        for (unsigned int node_i = 0; node_i < 2; node_i++)
        {
          Element* e = en->elem[node_i];
          if(e == NULL) continue;

          if(e->cm == NULL)
          {
            e->cm = new CurvMap;
            memset(e->cm, 0, sizeof(CurvMap));
            e->cm->toplevel = 1;
            e->cm->order = 4;
          }

          int idx = -1;
          for (unsigned j = 0; j < e->get_nvert(); j++)
            if(e->en[j] == en) { idx = j; break; }
            assert(idx >= 0);

            if(e->vn[idx]->id == p1)
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
        if(!nurbs->ref) delete nurbs;
      }

      // update refmap coeffs of curvilinear elements
      for_all_elements(e, mesh)
        if(e->cm != NULL)
          e->cm->update_refmap_coeffs(e);

      delete [] p1s;
      delete [] p2s;
      delete [] angles;

      // Refinements.
      bson_find(&it, &br, "refinement-count");
      int refinement_count = bson_iterator_int(&it);
      int* refined_elements = new int[refinement_count];
      int* refinement_types = new int[refinement_count];

      // refinement-id
      bson_find(&it_coeffs, &br, "refinement-id");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        refined_elements[index_coeff++] = bson_iterator_int(&it);
      // refinement-id
      bson_find(&it_coeffs, &br, "refinement-type");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        refinement_types[index_coeff++] = bson_iterator_int(&it);

      // perform initial refinements
      for (unsigned int refinement_i = 0; refinement_i < refinement_count; refinement_i++)
      {
        int element_id = refined_elements[refinement_i];
        int refinement_type = refinement_types[refinement_i];
        if(refinement_type == -1)
          mesh->unrefine_element_id(element_id);
        else
          mesh->refine_element_id(element_id, refinement_type);
      }

      delete [] refined_elements;
      delete [] refinement_types;

      mesh->initial_single_check();

      return true;
    }

    bool MeshReaderH2DBSON::save(const char *filename, MeshSharedPtr mesh)
    {
      throw Exceptions::MethodNotImplementedException("MeshReaderH2DBSON::save");

      // Utility pointer.
      Element* e;

      // Init bson
      bson bw;
      bson_init(&bw);

      // Save vertices
      // - count.
      bson_append_int(&bw, "vertex-count", mesh->ntopvert);
      // - xes
      bson_append_start_array(&bw, "vertex-x");
      for (int i = 0; i < mesh->ntopvert; i++)
        bson_append_double(&bw, "c", mesh->nodes[i].x);
      bson_append_finish_array(&bw);
      // - yes
      bson_append_start_array(&bw, "vertex-y");
      for (int i = 0; i < mesh->ntopvert; i++)
        bson_append_double(&bw, "c", mesh->nodes[i].y);
      bson_append_finish_array(&bw);

      // Save elements
      // - count.
      bson_append_int(&bw, "element-count", mesh->get_num_base_elements());
      // - first vertex ids  
      bson_append_start_array(&bw, "element-1");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        bson_append_int(&bw, "c", e->vn[0]->id);
      }
      bson_append_finish_array(&bw);
      // - second vertex ids  
      bson_append_start_array(&bw, "element-2");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        bson_append_int(&bw, "c", e->vn[1]->id);
      }
      bson_append_finish_array(&bw);
      // - third vertex ids  
      bson_append_start_array(&bw, "element-3");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        bson_append_int(&bw, "c", e->vn[2]->id);
      }
      bson_append_finish_array(&bw);
      // - first vertex ids  
      bson_append_start_array(&bw, "element-4");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        if(e->is_triangle())
          bson_append_int(&bw, "c", -1);
        else
          bson_append_int(&bw, "c", e->vn[3]->id);
      }
      bson_append_finish_array(&bw);
      // - markers
      bson_append_start_array(&bw, "element-marker");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        bson_append_string(&bw, "c", mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str());
      }
      bson_append_finish_array(&bw);

      // Save boundary markers
      // - count.
      int edge_count = 0;
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_nvert(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            edge_count++;
      bson_append_int(&bw, "edge-count", edge_count);
      // - ids 1
      bson_append_start_array(&bw, "edge-1");
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_nvert(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            bson_append_int(&bw, "c", e->vn[i]->id);
      bson_append_finish_array(&bw);
      // - ids 2
      bson_append_start_array(&bw, "edge-2");
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_nvert(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            bson_append_int(&bw, "c", e->vn[e->next_vert(i)]->id);
      bson_append_finish_array(&bw);
      // - markers
      bson_append_start_array(&bw, "edge-marker");
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        bson_append_string(&bw, "c", mesh->boundary_markers_conversion.get_user_marker(mesh->get_base_edge_node(e, i)->marker).marker.c_str());
      }
      bson_append_finish_array(&bw);

      // Save arcs.
      // - count.
      int arc_count = 0;
      for_all_base_elements(e, mesh)
      {
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
              if(e->cm->nurbs[i]->arc)
                arc_count++;
      }
      bson_append_int(&bw, "arc-count", arc_count);
      // - p1.
      bson_append_start_array(&bw, "arcs-p1");
      for_all_base_elements(e, mesh)
      {
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
            {
              if(!e->cm->nurbs[i]->arc)
                continue;
              bson_append_int(&bw, "c", e->vn[i]->id);
            }
      }
      bson_append_finish_array(&bw);
      // - p2.
      bson_append_start_array(&bw, "arcs-p2");
      for_all_base_elements(e, mesh)
      {
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
            {
              if(!e->cm->nurbs[i]->arc)
                continue;
              bson_append_int(&bw, "c", e->vn[e->next_vert(i)]->id);
            }
      }
      bson_append_finish_array(&bw);
      // - angles.
      bson_append_start_array(&bw, "arcs-angle");
      for_all_base_elements(e, mesh)
      {
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
            {
              if(!e->cm->nurbs[i]->arc)
                continue;
              bson_append_int(&bw, "c", e->cm->nurbs[i]->angle);
            }
      }
      bson_append_finish_array(&bw);

      // Save general NURBS.
      // - count.
      for_all_base_elements(e, mesh)
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
              if(!e->cm->nurbs[i]->arc)
                throw Exceptions::Exception("BSON mesh loader can not operate with general NURBS so far.");

      // Save refinements
      // - count.
      for(unsigned int refinement_i = 0; refinement_i < mesh->refinements.size(); refinement_i++)
        bson_append_int(&bw, "refinement-count", mesh->refinements.size());
      // - ids
      bson_append_start_array(&bw, "refinement-id");
      for(unsigned int refinement_i = 0; refinement_i < mesh->refinements.size(); refinement_i++)
        bson_append_int(&bw, "c", mesh->refinements[refinement_i].first);
      bson_append_finish_array(&bw);
      // - type
      bson_append_start_array(&bw, "refinement-type");
      for(unsigned int refinement_i = 0; refinement_i < mesh->refinements.size(); refinement_i++)
        bson_append_int(&bw, "c", mesh->refinements[refinement_i].second);
      bson_append_finish_array(&bw);

      // Done.
      bson_finish(&bw);

      // Write to disk.
      FILE *fpw;
      fpw = fopen(filename, "wb");
      const char *dataw = (const char *) bson_data(&bw);
      fwrite(dataw, bson_size(&bw), 1, fpw);
      fclose(fpw);

      bson_destroy(&bw);

      return true;
    }

    bool MeshReaderH2DBSON::load(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      throw Exceptions::MethodNotImplementedException("MeshReaderH2DBSON::load");

      return true;
    }

    bool MeshReaderH2DBSON::save(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      throw Exceptions::MethodNotImplementedException("MeshReaderH2DBSON::save");

      return true;
    }

    Nurbs* MeshReaderH2DBSON::load_arc(MeshSharedPtr mesh, int id, Node** en, int p1, int p2, double angle, bool skip_check)
    {
      Nurbs* nurbs = new Nurbs;
      nurbs->arc = true;

      *en = mesh->peek_edge_node(p1, p2);

      if(*en == NULL)
      {
        if(!skip_check)
          throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: edge %d-%d does not exist.", id, p1, p2);
        else
          return NULL;
      }

      // degree of an arc == 2.
      nurbs->degree = 2;
      // there are three control points.
      nurbs->np = 3;
      // there are 6 knots: {0, 0, 0, 1, 1, 1}
      nurbs->nk = 6;
      nurbs->kv = new double[nurbs->nk];

      for (int i = 0; i < 3; i++)
        nurbs->kv[i] = 0.0;

      for (int i = 3; i < nurbs->nk; i++)
        nurbs->kv[i] = 1.0;

      // edge endpoints control points.
      nurbs->pt = new double3[3];
      nurbs->pt[0][0] = mesh->nodes[p1].x;
      nurbs->pt[0][1] = mesh->nodes[p1].y;
      nurbs->pt[0][2] = 1.0;
      nurbs->pt[2][0] = mesh->nodes[p2].x;
      nurbs->pt[2][1] = mesh->nodes[p2].y;
      nurbs->pt[2][2] = 1.0;

      // read the arc angle
      nurbs->angle = angle;
      double a = (180.0 - nurbs->angle) / 180.0 * M_PI;

      // generate one inner control point
      double x = 1.0 / std::tan(a * 0.5);
      nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
      nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
      nurbs->pt[1][2] = Hermes::cos((M_PI - a) * 0.5);

      nurbs->ref = 0;

      return nurbs;
    }
  }
}

#endif