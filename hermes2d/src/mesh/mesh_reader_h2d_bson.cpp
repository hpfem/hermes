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

    void MeshReaderH2DBSON::load(const char *filename, MeshSharedPtr mesh)
    {
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
      bson_find(&it_coeffs, &br, "arcs-p1");
      bson_iterator_subobject_init(&it_coeffs, &sub, 0);
      bson_iterator_init(&it, &sub);
      index_coeff = 0;
      while (bson_iterator_next(&it))
        p1s[index_coeff++] = bson_iterator_int(&it);
      // p2s
      bson_find(&it_coeffs, &br, "arcs-p2");
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
      for_all_used_elements(e, mesh)
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

      bson_destroy(&br);
      ::free(datar);

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

      if(HermesCommonApi.get_integral_param_value(checkMeshesOnLoad))
        mesh->initial_single_check();
    }

    void MeshReaderH2DBSON::save(const char *filename, MeshSharedPtr mesh)
    {
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
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_nvert(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            bson_append_string(&bw, "c", mesh->boundary_markers_conversion.get_user_marker(mesh->get_base_edge_node(e, i)->marker).marker.c_str());
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
    }

    void MeshReaderH2DBSON::save(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      // For mapping of physical coordinates onto top vertices.
      std::map<std::pair<double, double>, unsigned int> points_to_vertices;
      // For mapping of vertex pairs onto boundary edges.
      std::map<std::pair<unsigned int, unsigned int>, unsigned int> vertices_to_boundaries;
      // For mapping of vertex pairs onto curves.
      std::map<std::pair<unsigned int, unsigned int>, bool> vertices_to_curves;

      Hermes::vector<element_BSON> elements;
      Hermes::vector<edge_BSON> edges;
      Hermes::vector<vertex_BSON> vertices;
      Hermes::vector<arc_BSON> arcs;
      Hermes::vector<subdomain_BSON> subdomains;

      bool* baseElementsSaved = new bool[meshes[0]->get_num_base_elements()];
      memset(baseElementsSaved, 0, sizeof(bool) * meshes[0]->get_num_base_elements());

#pragma region Save to structures

      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      {
        bool hasAllElements = (meshes[meshes_i]->get_num_used_base_elements() == meshes[meshes_i]->get_num_base_elements());

        subdomain_BSON subdomain;

        // Mapping of top vertices of subdomains to the global mesh.
        std::map<unsigned int, unsigned int> vertices_to_vertices;

        // Utility pointer.
        Element* e;

        // save vertices
        for (int i = 0; i < meshes[meshes_i]->ntopvert; i++)
        {
          // Look for the coordinates of this vertex.
          // If found, then insert the pair <this vertex number, the found vertex number> into vertices_to_vertices dictionary.
          // If not, insert.
          std::map<std::pair<double, double>, unsigned int>::iterator it = points_to_vertices.find(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y));
          if(it != points_to_vertices.end())
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, it->second));
          else
          {
            int new_i = points_to_vertices.size();
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, new_i));
            points_to_vertices.insert(std::pair<std::pair<double, double>, unsigned int>(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y), points_to_vertices.size()));

            vertices.push_back(vertex_BSON(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y, new_i));
          }
          if(!hasAllElements)
            subdomain.vertices.push_back(vertices_to_vertices.find(i)->second);
        }

        // save elements
        for (int i = 0; i < meshes[meshes_i]->get_num_base_elements(); i++)
        {
          e = &(meshes[meshes_i]->elements[i]);
          if(!e->used)
            continue;
          if(!baseElementsSaved[e->id])
          {
            if(e->nvert == 3)
              elements.push_back(element_BSON(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, -1, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker, e->id));
            else
              elements.push_back(element_BSON(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, vertices_to_vertices.find(e->vn[3]->id)->second, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker, e->id));
            baseElementsSaved[e->id] = true;
          }

          if(!hasAllElements)
            subdomain.elements.push_back(e->id);
        }

        // save boundary edge markers
        bool has_inner_edges = false;
        for_all_base_elements(e, meshes[meshes_i])
        {
          for (unsigned i = 0; i < e->get_nvert(); i++)
          {
            if(meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
            {
              if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
              {
                unsigned int edge_i = edges.size();
                vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                edges.push_back(edge_BSON(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker, edge_i));
              }
              if(!hasAllElements)
                subdomain.boundary_edges.push_back(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)))->second);
            }
            else
              has_inner_edges = true;
          }
        }

        if(has_inner_edges)
        {
          for_all_base_elements(e, meshes[meshes_i])
            for (unsigned i = 0; i < e->get_nvert(); i++)
            {
              if(!meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
              {
                if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
                {
                  unsigned int edge_i = edges.size();
                  vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                  edges.push_back(edge_BSON(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker, edge_i));
                }
                if(!hasAllElements)
                  subdomain.inner_edges.push_back(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)))->second);
              }
            }
        }

        // save curved edges
        for_all_base_elements(e, meshes[meshes_i])
          if(e->is_curved())
            for (unsigned i = 0; i < e->get_nvert(); i++)
              if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
                if(vertices_to_curves.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_curves.end())
                {
                  if(e->cm->nurbs[i]->arc)
                    arcs.push_back(arc_BSON(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, e->cm->nurbs[i]->angle));
                  else
                    throw Exceptions::Exception("BSON mesh loader can not operate with general NURBS so far.");

                  vertices_to_curves.insert(std::pair<std::pair<unsigned int, unsigned int>, bool>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), true));
                }

                // save refinements
                for(unsigned int refinement_i = 0; refinement_i < meshes[meshes_i]->refinements.size(); refinement_i++)
                  subdomain.refinements.push_back(refinement_BSON(meshes[meshes_i]->refinements[refinement_i].first, meshes[meshes_i]->refinements[refinement_i].second));

                subdomains.push_back(subdomain);
      }

      delete [] baseElementsSaved;

      std::sort (elements.begin(), elements.end(), elementCompare);
#pragma endregion

      // Utility pointer.
      Element* e;

      // Init bson
      bson bw;
      bson_init(&bw);

      // Save domain
      // - vertices.
      bson_append_start_array(&bw, "vertices");
      for (int i = 0; i < vertices.size(); i++)
        vertices[i].save_to_BSON(bw);
      bson_append_finish_array(&bw);
      // - elements.
      bson_append_start_array(&bw, "elements");
      for (int i = 0; i < elements.size(); i++)
        elements[i].save_to_BSON(bw);
      bson_append_finish_array(&bw);
      // - edges
      bson_append_start_array(&bw, "edges");
      for (int i = 0; i < edges.size(); i++)
        edges[i].save_to_BSON(bw);
      bson_append_finish_array(&bw);
      // - curves
      bson_append_start_array(&bw, "arcs");
      for (int i = 0; i < arcs.size(); i++)
        arcs[i].save_to_BSON(bw);
      bson_append_finish_array(&bw);

      // Save subdomains
      bson_append_start_array(&bw, "subdomains");
      for (int i = 0; i < subdomains.size(); i++)
        subdomains[i].save_to_BSON(bw);
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
    }

    void MeshReaderH2DBSON::load(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
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

      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
        meshes.at(meshes_i)->free();

      MeshSharedPtr global_mesh(new Mesh);

      std::map<int, int> vertex_is;
      std::map<int, int> element_is;
      std::map<int, int> edge_is;

      Hermes::vector<element_BSON> elements;
      Hermes::vector<edge_BSON> edges;
      Hermes::vector<vertex_BSON> vertices;
      Hermes::vector<arc_BSON> arcs;
      Hermes::vector<subdomain_BSON> subdomains;

      // load
      this->load_domain(br, global_mesh, vertex_is, element_is, edge_is, elements, edges, vertices, arcs, subdomains);

      int max_vertex_i = -1;
      for(std::map<int, int>::iterator it = vertex_is.begin(); it != vertex_is.end(); it++)
        if(it->first > max_vertex_i)
          max_vertex_i = it->first;
      int max_element_i = -1;
      for(std::map<int, int>::iterator it = element_is.begin(); it != element_is.end(); it++)
        if(it->first > max_element_i)
          max_element_i = it->first;
      int max_edge_i = -1;
      for(std::map<int, int>::iterator it = edge_is.begin(); it != edge_is.end(); it++)
        if(it->first > max_edge_i)
          max_edge_i = it->first;

      // Subdomains //
      unsigned int subdomains_count = subdomains.size();
      if(subdomains_count != meshes.size())
        throw Hermes::Exceptions::MeshLoadFailureException("Number of subdomains( = %u) does not equal the number of provided meshes in the vector( = %u).", subdomains_count, meshes.size());

      for(unsigned int subdomains_i = 0; subdomains_i < subdomains_count; subdomains_i++)
      {
        for (int element_i = 0; element_i < elements.size(); element_i++)
        {
          element_BSON element = elements.at(element_i);

          meshes[subdomains_i]->element_markers_conversion.insert_marker(element.marker);
        }
        for(unsigned int edge_i = 0; edge_i < edges.size(); edge_i++)
        {
          edge_BSON edge = edges.at(edge_i);

          meshes[subdomains_i]->boundary_markers_conversion.insert_marker(edge.marker);
        }
      }

      for(unsigned int subdomains_i = 0; subdomains_i < subdomains_count; subdomains_i++)
      {
        unsigned int vertex_number_count = subdomains.at(subdomains_i).vertices.empty() ? 0 : subdomains.at(subdomains_i).vertices.size();
        unsigned int element_number_count = subdomains.at(subdomains_i).elements.empty() ? 0 : subdomains.at(subdomains_i).elements.size();
        unsigned int boundary_edge_number_count = subdomains.at(subdomains_i).boundary_edges.empty() ? 0 : subdomains.at(subdomains_i).boundary_edges.size();
        unsigned int inner_edge_number_count = subdomains.at(subdomains_i).inner_edges.empty() ? 0 : subdomains.at(subdomains_i).inner_edges.size();

        // copy the whole mesh if the subdomain is the whole mesh.
        if(element_number_count == 0 || element_number_count == elements.size())
        {
          meshes[subdomains_i]->copy(global_mesh);
          // refinements.
          if(!subdomains.at(subdomains_i).refinements.empty() && subdomains.at(subdomains_i).refinements.size() > 0)
          {
            // perform initial refinements
            for (unsigned int i = 0; i < subdomains.at(subdomains_i).refinements.size(); i++)
            {
              int element_id = subdomains.at(subdomains_i).refinements.at(i).id;
              int refinement_type = subdomains.at(subdomains_i).refinements.at(i).type;
              if(refinement_type == -1)
                meshes[subdomains_i]->unrefine_element_id(element_id);
              else
                meshes[subdomains_i]->refine_element_id(element_id, refinement_type);
            }
          }
        }
        else
        {
          // Vertex numbers //
          // create a mapping order-in-the-whole-domain <-> order-in-this-subdomain.
          std::map<unsigned int, unsigned int> vertex_vertex_numbers;

          // Initialize mesh.
          int size = HashTable::H2D_DEFAULT_HASH_SIZE;
          while (size < 8 * vertex_number_count)
            size *= 2;
          meshes[subdomains_i]->init(size);

          // Create top-level vertex nodes.
          if(vertex_number_count == 0)
            vertex_number_count = vertices.size();
          for (unsigned int vertex_numbers_i = 0; vertex_numbers_i < vertex_number_count; vertex_numbers_i++)
          {
            unsigned int vertex_number;
            if(vertex_number_count == vertices.size())
              vertex_number = vertex_is[vertex_numbers_i];
            else
            {
              vertex_number =  subdomains.at(subdomains_i).vertices.at(vertex_numbers_i);
              if(vertex_number > max_vertex_i)
                throw Exceptions::MeshLoadFailureException("Wrong vertex number:%u in subdomain %u.", vertex_number, subdomains_i);
            }

            vertex_vertex_numbers.insert(std::pair<unsigned int, unsigned int>(vertex_number, vertex_numbers_i));
            Node* node = meshes[subdomains_i]->nodes.add();
            assert(node->id == vertex_numbers_i);
            node->ref = TOP_LEVEL_REF;
            node->type = HERMES_TYPE_VERTEX;
            node->bnd = 0;
            node->p1 = node->p2 = -1;
            node->next_hash = NULL;

            // assignment.
            node->x = vertices[vertex_number].x;
            node->y = vertices[vertex_number].y;
          }
          meshes[subdomains_i]->ntopvert = vertex_number_count;

          // Element numbers //
          unsigned int element_count = elements.size();
          meshes[subdomains_i]->nbase = element_count;
          meshes[subdomains_i]->nactive = meshes[subdomains_i]->ninitial = element_number_count;

          Element* e;
          int* elements_existing = new int[element_count];
          for(int i = 0; i < element_count; i++)
            elements_existing[i] = -1;
          for (int element_number_i = 0; element_number_i < element_number_count; element_number_i++)
          {
            int elementI = subdomains.at(subdomains_i).elements.at(element_number_i);
            if(elementI > max_element_i)
              throw Exceptions::MeshLoadFailureException("Wrong element number:%i in subdomain %u.", elementI, subdomains_i);

            elements_existing[element_is[subdomains.at(subdomains_i).elements.at(element_number_i)]] = elementI;
          }
          for (int element_i = 0; element_i < element_count; element_i++)
          {
            bool found = false;
            if(element_number_count == 0)
              found = true;
            else
              found = elements_existing[element_i] != -1;

            if(!found)
            {
              meshes[subdomains_i]->elements.skip_slot();
              continue;
            }

            element_BSON element(-1,-1,-1,-1, "", -1);
            for(int searched_element_i = 0; searched_element_i < element_count; searched_element_i++)
            {
              element = elements.at(searched_element_i);
              if(element.i == elements_existing[element_i])
                break;
            }
            if(element.i == -1)
              throw Exceptions::MeshLoadFailureException("Element number wrong in the mesh file.");

            if(element.v4 != -1)
              e = meshes[subdomains_i]->create_quad(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element.marker).marker,
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v1)->second],
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v2)->second],
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v3)->second],
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v4)->second],
              NULL, element_i);
            else
              e = meshes[subdomains_i]->create_triangle(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element.marker).marker,
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v1)->second],
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v2)->second],
              &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(element.v3)->second],
              NULL, element_i);
          }

          // Boundary Edge numbers //
          if(boundary_edge_number_count == 0)
            boundary_edge_number_count = edges.size();

          for (int boundary_edge_number_i = 0; boundary_edge_number_i < boundary_edge_number_count; boundary_edge_number_i++)
          {
            edge_BSON edge(-1,-1, "", -1);
            int domain_edge_count = edges.size();
            for(unsigned int to_find_i = 0; to_find_i < domain_edge_count; to_find_i++)
            {
              if(boundary_edge_number_count != domain_edge_count)
              {
                if(edges.at(to_find_i).i == subdomains.at(subdomains_i).boundary_edges.at(boundary_edge_number_i))
                {
                  edge = edges.at(to_find_i);
                  break;
                }
              }
              else
              {
                if(edges.at(to_find_i).i == edge_is[boundary_edge_number_i])
                {
                  edge = edges.at(to_find_i);
                  break;
                }
              }
            }

            if(edge.i == -1)
              throw Exceptions::MeshLoadFailureException("Wrong boundary-edge number:%i in subdomain %u.", subdomains.at(subdomains_i).boundary_edges.at(boundary_edge_number_i), subdomains_i);

            Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge.v1)->second, vertex_vertex_numbers.find(edge.v2)->second);
            if(en == NULL)
              throw Hermes::Exceptions::MeshLoadFailureException("Boundary data error (edge %i does not exist).", boundary_edge_number_i);

            en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge.marker).marker;

            meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge.v1)->second].bnd = 1;
            meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge.v2)->second].bnd = 1;
            en->bnd = 1;
          }

          // Inner Edge numbers //
          for (int inner_edge_number_i = 0; inner_edge_number_i < inner_edge_number_count; inner_edge_number_i++)
          {
            edge_BSON edge(-1, -1, "", -1);

            for(unsigned int to_find_i = 0; to_find_i < edges.size(); to_find_i++)
            {
              if(edges.at(to_find_i).i == subdomains.at(subdomains_i).inner_edges.at(inner_edge_number_i))
              {
                edge = edges.at(to_find_i);
                break;
              }
            }

            if(edge.i == -1)
              throw Exceptions::MeshLoadFailureException("Wrong inner-edge number:%i in subdomain %u.", subdomains.at(subdomains_i).boundary_edges.at(inner_edge_number_i), subdomains_i);

            Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge.v1)->second, vertex_vertex_numbers.find(edge.v2)->second);
            if(en == NULL)
              throw Hermes::Exceptions::MeshLoadFailureException("Inner data error (edge %i does not exist).", inner_edge_number_i);

            en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge.marker).marker;
            en->bnd = 0;
          }

          // Curves //
          // Arcs & NURBSs //
          unsigned int arc_count = arcs.empty() ? 0 : arcs.size();

          for (unsigned int curves_i = 0; curves_i < arc_count; curves_i++)
          {
            // load the control points, knot vector, etc.
            Node* en;
            int p1, p2;

            // first do arcs, then NURBSs.
            Nurbs* nurbs;
            if(vertex_vertex_numbers.find(arcs.at(curves_i).p1) == vertex_vertex_numbers.end() ||
              vertex_vertex_numbers.find(arcs.at(curves_i).p2) == vertex_vertex_numbers.end())
              continue;
            else
            {
              // read the end point indices
              p1 = vertex_vertex_numbers.find(arcs.at(curves_i).p1)->second;
              p2 = vertex_vertex_numbers.find(arcs.at(curves_i).p2)->second;

              nurbs = load_arc(meshes[subdomains_i], curves_i, &en, p1, p2, arcs.at(curves_i).angle, true);
              if(nurbs == NULL)
                continue;
            }

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
                  Nurbs* nurbs_rev = meshes[subdomains_i]->reverse_nurbs(nurbs);
                  e->cm->nurbs[idx] = nurbs_rev;
                  nurbs_rev->ref++;
                }
            }
            if(!nurbs->ref) delete nurbs;
          }

          // update refmap coeffs of curvilinear elements
          for_all_used_elements(e, meshes[subdomains_i])
            if(e->cm != NULL)
              e->cm->update_refmap_coeffs(e);

          // refinements.
          if(!subdomains.at(subdomains_i).refinements.empty() && subdomains.at(subdomains_i).refinements.size() > 0)
          {
            // perform initial refinements
            for (unsigned int i = 0; i < subdomains.at(subdomains_i).refinements.size(); i++)
            {
              int element_id = subdomains.at(subdomains_i).refinements.at(i).id;
              int refinement_type = subdomains.at(subdomains_i).refinements.at(i).type;
              if(refinement_type == -1)
                meshes[subdomains_i]->unrefine_element_id(element_id);
              else
                meshes[subdomains_i]->refine_element_id(element_id, refinement_type);
            }
          }

          delete [] elements_existing;
        }
        meshes[subdomains_i]->seq = g_mesh_seq++;
        if(HermesCommonApi.get_integral_param_value(checkMeshesOnLoad))
          meshes[subdomains_i]->initial_single_check();
      }
      bson_destroy(&br);
      ::free(datar);
    }

    void MeshReaderH2DBSON::load_domain(bson& br, MeshSharedPtr mesh, std::map<int, int>& vertex_is, std::map<int, int>& element_is, std::map<int, int>& edge_is,
      Hermes::vector<element_BSON>& elements, Hermes::vector<edge_BSON>& edges, Hermes::vector<vertex_BSON>& vertices, Hermes::vector<arc_BSON>& arcs, Hermes::vector<subdomain_BSON>& subdomains)
    {
      bson_iterator it, it_sub;
      bson b_sub, b_sub_sub;
      bson_find(&it, &br, "vertices");
      bson_iterator_subobject_init(&it, &b_sub, 0);
      bson_iterator_init(&it_sub, &b_sub);
      while (bson_iterator_next(&it_sub))
      {
        vertex_BSON vertex;
        bson_iterator_subobject_init(&it_sub, &b_sub_sub, 0);
        vertex.load_from_BSON(b_sub_sub);
        vertices.push_back(vertex);
      }

      bson_find(&it, &br, "elements");
      bson_iterator_subobject_init(&it, &b_sub, 0);
      bson_iterator_init(&it_sub, &b_sub);
      while (bson_iterator_next(&it_sub))
      {
        element_BSON element;
        bson_iterator_subobject_init(&it_sub, &b_sub_sub, 0);
        element.load_from_BSON(b_sub_sub);
        elements.push_back(element);
      }

      bson_find(&it, &br, "edges");
      bson_iterator_subobject_init(&it, &b_sub, 0);
      bson_iterator_init(&it_sub, &b_sub);
      while (bson_iterator_next(&it_sub))
      {
        edge_BSON edge;
        bson_iterator_subobject_init(&it_sub, &b_sub_sub, 0);
        edge.load_from_BSON(b_sub_sub);
        edges.push_back(edge);
      }

      bson_find(&it, &br, "arcs");
      bson_iterator_subobject_init(&it, &b_sub, 0);
      bson_iterator_init(&it_sub, &b_sub);
      while (bson_iterator_next(&it_sub))
      {
        arc_BSON arc;
        bson_iterator_subobject_init(&it_sub, &b_sub_sub, 0);
        arc.load_from_BSON(b_sub_sub);
        arcs.push_back(arc);
      }

      bson_find(&it, &br, "subdomains");
      bson_iterator_subobject_init(&it, &b_sub, 0);
      bson_iterator_init(&it_sub, &b_sub);
      while (bson_iterator_next(&it_sub))
      {
        subdomain_BSON subdomain;
        bson_iterator_subobject_init(&it_sub, &b_sub_sub, 0);
        subdomain.load_from_BSON(b_sub_sub);
        subdomains.push_back(subdomain);
      }

      // Vertices //
      int vertices_count = vertices.size();

      // Initialize mesh.
      int size = HashTable::H2D_DEFAULT_HASH_SIZE;
      while (size < 8 * vertices_count)
        size *= 2;
      mesh->init(size);

      // Create top-level vertex nodes.
      for (int vertex_i = 0; vertex_i < vertices_count; vertex_i++)
      {
        Node* node = mesh->nodes.add();
        assert(node->id == vertex_i);
        node->ref = TOP_LEVEL_REF;
        node->type = HERMES_TYPE_VERTEX;
        node->bnd = 0;
        node->p1 = node->p2 = -1;
        node->next_hash = NULL;

        if(vertices[vertex_i].i > H2D_MAX_NODE_ID - 1)
          throw Exceptions::MeshLoadFailureException("The index 'i' of vertex in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

        // insert
        vertex_is.insert(std::pair<int,int>(vertices[vertex_i].i, vertex_i));

        // assignment.
        node->x = vertices[vertex_i].x;
        node->y = vertices[vertex_i].y;
      }
      mesh->ntopvert = vertices_count;

      // Elements //
      unsigned int element_count = elements.size();
      mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

      Element* e;
      for (int element_i = 0; element_i < element_count; element_i++)
      {
        element_BSON element = elements[element_i];

        // insert.
        if(element.i > H2D_MAX_NODE_ID - 1)
          throw Exceptions::MeshLoadFailureException("The index 'i' of element in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

        element_is.insert(std::pair<int,int>(element.i, element_i));

        mesh->element_markers_conversion.insert_marker(element.marker);

        if(element.v4 != -1)
          e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(element.marker).marker,
          &mesh->nodes[element.v1],
          &mesh->nodes[element.v2],
          &mesh->nodes[element.v3],
          &mesh->nodes[element.v4],
          NULL);
        else
          e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(element.marker).marker,
          &mesh->nodes[element.v1],
          &mesh->nodes[element.v2],
          &mesh->nodes[element.v3],
          NULL);
      }

      // Boundaries //
      unsigned int edges_count = edges.size();

      Node* en;
      for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
      {
        int v1 = edges[edge_i].v1;
        int v2 = edges[edge_i].v2;

        // insert
        if(edges[edge_i].i > H2D_MAX_NODE_ID - 1)
          throw Exceptions::MeshLoadFailureException("The index 'i' of edge in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

        edge_is.insert(std::pair<int, int>(edge_i, edges[edge_i].i));

        en = mesh->peek_edge_node(v1, v2);
        if(en == NULL)
          throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist.", edge_i, v1, v2);

        std::string edge_marker = edges[edge_i].marker;

        // Trim whitespaces.
        unsigned int begin = edge_marker.find_first_not_of(" \t\n");
        unsigned int end = edge_marker.find_last_not_of(" \t\n");
        edge_marker.erase(end + 1, edge_marker.length());
        edge_marker.erase(0, begin);

        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker(edge_marker);
        int marker = mesh->boundary_markers_conversion.get_internal_marker(edge_marker).marker;

        en->marker = marker;

        // This is extremely important, as in DG, it is assumed that negative boundary markers are reserved
        // for the inner edges.
        if(marker > 0)
        {
          mesh->nodes[v1].bnd = 1;
          mesh->nodes[v2].bnd = 1;
          en->bnd = 1;
        }
      }

      // check that all boundary edges have a marker assigned
      for_all_edge_nodes(en, mesh)
        if(en->ref < 2 && en->marker == 0)
          this->warn("Boundary edge node does not have a boundary marker.");

      // Curves //
      // Arcs & NURBSs //
      unsigned int arc_count = arcs.size();

      for (unsigned int curves_i = 0; curves_i < arc_count; curves_i++)
      {
        // load the control points, knot vector, etc.
        Node* en;
        int p1, p2;

        // first do arcs, then NURBSs.
        Nurbs* nurbs;
        // read the end point indices
        p1 = arcs[curves_i].p1;
        p2 = arcs[curves_i].p2;

        nurbs = load_arc(mesh, curves_i, &en, p1, p2, arcs[curves_i].angle);

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
      for_all_used_elements(e, mesh)
        if(e->cm != NULL)
          e->cm->update_refmap_coeffs(e);


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