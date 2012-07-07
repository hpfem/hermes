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

#include "mesh.h"
#include "mesh_reader_h2d_xml.h"
#include <iostream>

using namespace std;

namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderH2DXML::MeshReaderH2DXML()
    {
    }

    MeshReaderH2DXML::~MeshReaderH2DXML()
    {
    }

    bool MeshReaderH2DXML::load(const char *filename, Mesh *mesh)
    {
      mesh->free();

      std::map<unsigned int, unsigned int> vertex_is;

      try
      {
        std::auto_ptr<XMLMesh::mesh> parsed_xml_mesh(XMLMesh::mesh_(filename));

        if(!load(parsed_xml_mesh, mesh, vertex_is))
          return false;

        // refinements.
        if(parsed_xml_mesh->refinements().present() && parsed_xml_mesh->refinements()->refinement().size() > 0)
        {
          // perform initial refinements
          for (unsigned int i = 0; i < parsed_xml_mesh->refinements()->refinement().size(); i++)
          {
            int element_id = parsed_xml_mesh->refinements()->refinement().at(i).element_id();
            int refinement_type = parsed_xml_mesh->refinements()->refinement().at(i).refinement_type();
            if(refinement_type == -1)
              mesh->unrefine_element_id(element_id);
            else
              mesh->refine_element_id(element_id, refinement_type);
          }
        }
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    bool MeshReaderH2DXML::save(const char *filename, Mesh *mesh)
    {
      // Utility pointer.
      Element* e;

      // save vertices
      XMLMesh::vertices_type vertices;
      for (int i = 0; i < mesh->ntopvert; i++)
      {
        std::ostringstream x_stream;
        x_stream << mesh->nodes[i].x;

        std::ostringstream y_stream;
        y_stream << mesh->nodes[i].y;

        vertices.vertex().push_back(std::auto_ptr<XMLMesh::vertex>(new XMLMesh::vertex(x_stream.str(), y_stream.str(), i)));
      }

      // save elements
      XMLMesh::elements_type elements;
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        if(e->used)
          if(e->is_triangle())
            elements.element().push_back(XMLMesh::triangle_type(e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str()));
          else
            elements.element().push_back(XMLMesh::quad_type(e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->vn[3]->id));
      }
      // save boundary markers
      XMLMesh::edges_type edges;
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_num_surf(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            edges.edge().push_back(XMLMesh::edge(e->vn[i]->id, e->vn[e->next_vert(i)]->id, mesh->boundary_markers_conversion.get_user_marker(mesh->get_base_edge_node(e, i)->marker).marker.c_str()));

      // save curved edges
      XMLMesh::curves_type curves;
      for_all_base_elements(e, mesh)
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_num_surf(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
              if(e->cm->nurbs[i]->arc)
                save_arc(mesh, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i], curves);
              else
                save_nurbs(mesh, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i], curves);

      // save refinements
      XMLMesh::refinements_type refinements;
      for(unsigned int refinement_i = 0; refinement_i < mesh->refinements.size(); refinement_i++)
        refinements.refinement().push_back(XMLMesh::refinement(mesh->refinements[refinement_i].first, mesh->refinements[refinement_i].second));

      XMLMesh::mesh xmlmesh(vertices, elements, edges);
      xmlmesh.curves().set(curves);
      xmlmesh.refinements().set(refinements);

      std::string mesh_schema_location(H2D_XML_SCHEMAS_DIRECTORY);
      mesh_schema_location.append("/mesh_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_mesh("XMLMesh", mesh_schema_location);

      ::xml_schema::namespace_infomap namespace_info_map;
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("mesh", namespace_info_mesh));

      std::ofstream out(filename);
      XMLMesh::mesh_(out, xmlmesh, namespace_info_map);
      out.close();

      return true;
    }

    bool MeshReaderH2DXML::load(const char *filename, Hermes::vector<Mesh *> meshes)
    {
      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
        meshes.at(meshes_i)->free();

      Mesh global_mesh;

      try
      {
        std::auto_ptr<XMLSubdomains::domain> parsed_xml_domain (XMLSubdomains::domain_(filename));

        int* vertex_is = new int[H2D_MAX_NODE_ID];
        for(int i = 0; i < H2D_MAX_NODE_ID; i++)
          vertex_is[i] = -1;

        int* element_is = new int[H2D_MAX_NODE_ID];
        for(int i = 0; i < H2D_MAX_NODE_ID; i++)
          element_is[i] = -1;

        int* edge_is = new int[H2D_MAX_NODE_ID];
        for(int i = 0; i < H2D_MAX_NODE_ID; i++)
          edge_is[i] = -1;

        if(!load(parsed_xml_domain, &global_mesh, vertex_is, element_is, edge_is))
          return false;

        // Subdomains //
        unsigned int subdomains_count = parsed_xml_domain->subdomains().subdomain().size();
        if(subdomains_count != meshes.size())
          throw Hermes::Exceptions::MeshLoadFailureException("Number of subdomains( = %u) does not equal the number of provided meshes in the vector( = %u).", subdomains_count, meshes.size());

        for(unsigned int subdomains_i = 0; subdomains_i < subdomains_count; subdomains_i++)
        {
          unsigned int vertex_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices()->i().size() : 0;
          unsigned int element_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements()->i().size() : 0;
          unsigned int boundary_edge_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().size() : 0;
          unsigned int inner_edge_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges()->i().size() : 0;

          // copy nodes and elements
          if(vertex_number_count == 0 && element_number_count == 0 && boundary_edge_number_count == 0)
          {
            meshes[subdomains_i]->copy(&global_mesh);
            continue;
          }
          else
          {
            // Variables //
            unsigned int variables_count = parsed_xml_domain->variables().present() ? parsed_xml_domain->variables()->variable().size() : 0;

            std::map<std::string, double> variables;
            for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
              variables.insert(std::make_pair<std::string, double>(parsed_xml_domain->variables()->variable().at(variables_i).name(), parsed_xml_domain->variables()->variable().at(variables_i).value()));

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
              vertex_number_count = parsed_xml_domain->vertices().vertex().size();
            for (unsigned int vertex_numbers_i = 0; vertex_numbers_i < vertex_number_count; vertex_numbers_i++)
            {
              unsigned int vertex_number;
              if(vertex_number_count == parsed_xml_domain->vertices().vertex().size())
                vertex_number = vertex_is[vertex_numbers_i];
              else
                vertex_number =  vertex_is[parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices()->i().at(vertex_numbers_i)];

              vertex_vertex_numbers.insert(std::pair<unsigned int, unsigned int>(vertex_number, vertex_numbers_i));
              Node* node = meshes[subdomains_i]->nodes.add();
              assert(node->id == vertex_numbers_i);
              node->ref = TOP_LEVEL_REF;
              node->type = HERMES_TYPE_VERTEX;
              node->bnd = 0;
              node->p1 = node->p2 = -1;
              node->next_hash = NULL;

              // variables matching.
              std::string x = parsed_xml_domain->vertices().vertex().at(vertex_number).x();
              std::string y = parsed_xml_domain->vertices().vertex().at(vertex_number).y();
              double x_value;
              double y_value;

              // variables lookup.
              bool x_found = false;
              bool y_found = false;
              if(variables.find(x) != variables.end())
              {
                x_value = variables.find(x)->second;
                x_found = true;
              }
              if(variables.find(y) != variables.end())
              {
                y_value = variables.find(y)->second;
                y_found = true;
              }

              // test of value if no variable found.
              if(!x_found)
                if(std::strtod(x.c_str(), NULL) != 0.0)
                  x_value = std::strtod(x.c_str(), NULL);
                else
                {
                  // This is a hard part, to find out if it is really zero.
                  int dot_position = strchr(x.c_str(), '.') == NULL ? -1 : strchr(x.c_str(), '.') - x.c_str();
                  for(int i = 0; i < dot_position; i++)
                    if(strncmp(x.c_str() + i, "0", 1) != 0)
                      throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_number + 1);
                  for(int i = dot_position + 1; i < x.length(); i++)
                    if(strncmp(x.c_str() + i, "0", 1) != 0)
                      throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_number + 1);
                  x_value = std::strtod(x.c_str(), NULL);
                }

                if(!y_found)
                  if(std::strtod(y.c_str(), NULL) != 0.0)
                    y_value = std::strtod(y.c_str(), NULL);
                  else
                  {
                    // This is a hard part, to find out if it is really zero.
                    int dot_position = strchr(y.c_str(), '.') == NULL ? -1 : strchr(y.c_str(), '.') - y.c_str();
                    for(int i = 0; i < dot_position; i++)
                      if(strncmp(y.c_str() + i, "0", 1) != 0)
                        throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_number + 1);
                    for(int i = dot_position + 1; i < y.length(); i++)
                      if(strncmp(y.c_str() + i, "0", 1) != 0)
                        throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_number + 1);
                    y_value = std::strtod(y.c_str(), NULL);
                  }

                  // assignment.
                  node->x = x_value;
                  node->y = y_value;
            }
            meshes[subdomains_i]->ntopvert = vertex_number_count;

            // Element numbers //
            unsigned int element_count = parsed_xml_domain->elements().element().size();
            meshes[subdomains_i]->nbase = element_count;
            meshes[subdomains_i]->nactive = meshes[subdomains_i]->ninitial = element_number_count;

            Element* e;
            bool* elements_existing = new bool[element_count];
            for(int i = 0; i < element_count; i++)
              elements_existing[i] = false;
            for (int element_number_i = 0; element_number_i < element_number_count; element_number_i++)
              elements_existing[element_is[parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements()->i().at(element_number_i)]] = true;

            for (int element_i = 0; element_i < element_count; element_i++)
            {
              bool found = false;
              if(element_number_count == 0)
                found = true;
              else
                found = elements_existing[element_i];

              if(!found)
              {
                meshes[subdomains_i]->elements.skip_slot();
                continue;
              }

              XMLSubdomains::domain::elements_type::element_type* element = &parsed_xml_domain->elements().element().at(element_i);

              // Trim whitespaces.
              unsigned int begin = element->marker().find_first_not_of(" \t\n");
              unsigned int end = element->marker().find_last_not_of(" \t\n");
              element->marker().erase(end + 1, element->marker().length());
              element->marker().erase(0, begin);

              meshes[subdomains_i]->element_markers_conversion.insert_marker(meshes[subdomains_i]->element_markers_conversion.min_marker_unused, element->marker());

              if(dynamic_cast<XMLSubdomains::quad_type*>(element) != NULL)
                e = meshes[subdomains_i]->create_quad(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element->marker()).marker,
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::quad_type*>(element)->v1())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::quad_type*>(element)->v2())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::quad_type*>(element)->v3())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::quad_type*>(element)->v4())->second],
                NULL);
              else
                e = meshes[subdomains_i]->create_triangle(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element->marker()).marker,
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::triangle_type*>(element)->v1())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::triangle_type*>(element)->v2())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(dynamic_cast<XMLSubdomains::triangle_type*>(element)->v3())->second],
                NULL);
            }

            // Boundary Edge numbers //
            if(boundary_edge_number_count == 0)
              boundary_edge_number_count = parsed_xml_domain->edges().edge().size();

            for (int boundary_edge_number_i = 0; boundary_edge_number_i < boundary_edge_number_count; boundary_edge_number_i++)
            {
              XMLSubdomains::domain::edges_type::edge_type* edge;
              for(unsigned int to_find_i = 0; to_find_i < parsed_xml_domain->edges().edge().size(); to_find_i++)
              {
                if(boundary_edge_number_count != parsed_xml_domain->edges().edge().size())
                {
                  if(parsed_xml_domain->edges().edge().at(to_find_i).i() == parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().at(boundary_edge_number_i))
                  {
                    edge = &parsed_xml_domain->edges().edge().at(to_find_i);
                    break;
                  }
                }
                else
                {
                  if(parsed_xml_domain->edges().edge().at(to_find_i).i() == edge_is[boundary_edge_number_i])
                  {
                    edge = &parsed_xml_domain->edges().edge().at(to_find_i);
                    break;
                  }
                }
              }

              Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge->v1())->second, vertex_vertex_numbers.find(edge->v2())->second);
              if(en == NULL)
                throw Hermes::Exceptions::MeshLoadFailureException("Boundary data error (edge %i does not exist)", boundary_edge_number_i);

              // Trim whitespaces.
              unsigned int begin = edge->marker().find_first_not_of(" \t\n");
              unsigned int end = edge->marker().find_last_not_of(" \t\n");
              edge->marker().erase(end + 1, edge->marker().length());
              edge->marker().erase(0, begin);

              meshes[subdomains_i]->boundary_markers_conversion.insert_marker(meshes[subdomains_i]->boundary_markers_conversion.min_marker_unused, edge->marker());

              en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge->marker()).marker;

              meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge->v1())->second].bnd = 1;
              meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge->v2())->second].bnd = 1;
              en->bnd = 1;
            }

            // Inner Edge numbers //
            for (int inner_edge_number_i = 0; inner_edge_number_i < inner_edge_number_count; inner_edge_number_i++)
            {
              XMLSubdomains::domain::edges_type::edge_type* edge;

              for(unsigned int to_find_i = 0; to_find_i < parsed_xml_domain->edges().edge().size(); to_find_i++)
              {
                if(parsed_xml_domain->edges().edge().at(to_find_i).i() == parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges()->i().at(inner_edge_number_i))
                {
                  edge = &parsed_xml_domain->edges().edge().at(to_find_i);
                  break;
                }
              }

              Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge->v1())->second, vertex_vertex_numbers.find(edge->v2())->second);
              if(en == NULL)
                throw Hermes::Exceptions::MeshLoadFailureException("Inner data error (edge %i does not exist)", inner_edge_number_i);

              // Trim whitespaces.
              unsigned int begin = edge->marker().find_first_not_of(" \t\n");
              unsigned int end = edge->marker().find_last_not_of(" \t\n");
              edge->marker().erase(end + 1, edge->marker().length());
              edge->marker().erase(0, begin);

              meshes[subdomains_i]->boundary_markers_conversion.insert_marker(meshes[subdomains_i]->boundary_markers_conversion.min_marker_unused, edge->marker());

               en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge->marker()).marker;

               en->bnd = 0;
            }

            // Curves //
            // Arcs & NURBSs //
            unsigned int arc_count = parsed_xml_domain->curves().present() ? parsed_xml_domain->curves()->arc().size() : 0;
            unsigned int nurbs_count = parsed_xml_domain->curves().present() ? parsed_xml_domain->curves()->NURBS().size() : 0;

            for (unsigned int curves_i = 0; curves_i < arc_count + nurbs_count; curves_i++)
            {
              // load the control points, knot vector, etc.
              Node* en;
              int p1, p2;

              // first do arcs, then NURBSs.
              Nurbs* nurbs;
              if(curves_i < arc_count)
              {
                if(vertex_vertex_numbers.find(parsed_xml_domain->curves()->arc().at(curves_i).v1()) == vertex_vertex_numbers.end() ||
                  vertex_vertex_numbers.find(parsed_xml_domain->curves()->arc().at(curves_i).v2()) == vertex_vertex_numbers.end())
                  continue;
                else
                {
                  // read the end point indices
                  p1 = vertex_vertex_numbers.find(parsed_xml_domain->curves()->arc().at(curves_i).v1())->second;
                  p2 = vertex_vertex_numbers.find(parsed_xml_domain->curves()->arc().at(curves_i).v2())->second;

                  nurbs = load_arc(meshes[subdomains_i], parsed_xml_domain, curves_i, &en, p1, p2, true);
                  if(nurbs == NULL)
                    continue;
                }
              }
              else
              {
                if(vertex_vertex_numbers.find(parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v1()) == vertex_vertex_numbers.end() ||
                  vertex_vertex_numbers.find(parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v2()) == vertex_vertex_numbers.end())
                  continue;
                else
                {
                  // read the end point indices
                  p1 = vertex_vertex_numbers.find(parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v1())->second;
                  p2 = vertex_vertex_numbers.find(parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v2())->second;

                  nurbs = load_nurbs(meshes[subdomains_i], parsed_xml_domain, curves_i - arc_count, &en, p1, p2, true);
                  if(nurbs == NULL)
                    continue;
                }
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
                for (unsigned j = 0; j < e->get_num_surf(); j++)
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
            for_all_elements(e, meshes[subdomains_i])
              if(e->cm != NULL)
                e->cm->update_refmap_coeffs(e);

            // refinements.
            if(parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements().present() && parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->refinement().size() > 0)
            {
              // perform initial refinements
              for (unsigned int i = 0; i < parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->refinement().size(); i++)
              {
                int element_id = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->refinement().at(i).element_id();
                int refinement_type = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->refinement().at(i).refinement_type();
                if(refinement_type == -1)
                  meshes[subdomains_i]->unrefine_element_id(element_id);
                else
                  meshes[subdomains_i]->refine_element_id(element_id, refinement_type);
              }
            }

            meshes[subdomains_i]->seq = g_mesh_seq++;
            delete [] elements_existing;
          }
        }

        delete [] vertex_is;

        delete [] element_is;

        delete [] edge_is;

        return true;
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    bool MeshReaderH2DXML::save(const char *filename, Hermes::vector<Mesh *> meshes)
    {
      // For mapping of physical coordinates onto top vertices.
      std::map<std::pair<double, double>, unsigned int> points_to_vertices;
      // For mapping of vertex pairs onto boundary edges.
      std::map<std::pair<unsigned int, unsigned int>, unsigned int> vertices_to_boundaries;
      // For mapping of vertex pairs onto curves.
      std::map<std::pair<unsigned int, unsigned int>, bool> vertices_to_curves;

      // Global vertices list.
      XMLMesh::vertices_type vertices;
      // Global elements list.
      XMLSubdomains::elements_type elements;
      // Global boudnary edges list.
      XMLSubdomains::edges_type edges;
      // Global curves list.
      XMLMesh::curves_type curves;

      // Subdomains.
      XMLSubdomains::subdomains subdomains;

      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      {
        // Create a subdomain.
        XMLSubdomains::subdomain subdomain("A subdomain");

        // Refinements.
        XMLMesh::refinements_type refinements;

        // Mapping of top vertices of subdomains to the global mesh.
        std::map<unsigned int, unsigned int> vertices_to_vertices;

        // Utility pointer.
        Element* e;

        // save vertices
        subdomain.vertices().set(XMLSubdomains::subdomain::vertices_type());
        for (int i = 0; i < meshes[meshes_i]->ntopvert; i++)
        {
          // Look for the coordinates of this vertex.
          // If found, then insert the pair <this vertex number, the found vertex number> into vertices_to_vertices dictionary.
          // If not, insert.
          if(points_to_vertices.find(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y)) != points_to_vertices.end())
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, points_to_vertices.find(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y))->second));
          else
          {
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, points_to_vertices.size()));
            points_to_vertices.insert(std::pair<std::pair<double, double>, unsigned int>(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y), points_to_vertices.size()));
            std::ostringstream x_stream;
            x_stream << meshes[meshes_i]->nodes[vertices_to_vertices.find(i)->second].x;

            std::ostringstream y_stream;
            y_stream << meshes[meshes_i]->nodes[vertices_to_vertices.find(i)->second].y;

            vertices.vertex().push_back(std::auto_ptr<XMLMesh::vertex>(new XMLMesh::vertex(x_stream.str(), y_stream.str(), i)));
          }
          subdomain.vertices()->i().push_back(vertices_to_vertices.find(i)->second);
        }

        // save elements
        subdomain.elements().set(XMLSubdomains::subdomain::elements_type());
        for (int i = 0; i < meshes[meshes_i]->get_num_base_elements(); i++)
        {
          e = meshes[meshes_i]->get_element_fast(i);
          if(e->used)
          {
            if(e->is_triangle())
            {
              bool present = false;
              for(unsigned int elements_i = 0; elements_i < elements.element().size(); elements_i++)
                if(dynamic_cast<XMLSubdomains::triangle_type*>(&elements.element().at(elements_i)) != NULL)
                  if(dynamic_cast<XMLSubdomains::triangle_type*>(&elements.element().at(elements_i))->v1() == vertices_to_vertices.find(e->vn[0]->id)->second &&
                    dynamic_cast<XMLSubdomains::triangle_type*>(&elements.element().at(elements_i))->v2() == vertices_to_vertices.find(e->vn[1]->id)->second &&
                    dynamic_cast<XMLSubdomains::triangle_type*>(&elements.element().at(elements_i))->v3() == vertices_to_vertices.find(e->vn[2]->id)->second)
                    {
                      present = true;
                      break;
                    }

              if(!present)
                elements.element().push_back(XMLSubdomains::triangle_type(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->id));
            }
            else
            {
              bool present = false;
              for(unsigned int elements_i = 0; elements_i < elements.element().size(); elements_i++)
                if(dynamic_cast<XMLSubdomains::quad_type*>(&elements.element().at(elements_i)) != NULL)
                  if(dynamic_cast<XMLSubdomains::quad_type*>(&elements.element().at(elements_i))->v1() == vertices_to_vertices.find(e->vn[0]->id)->second &&
                    dynamic_cast<XMLSubdomains::quad_type*>(&elements.element().at(elements_i))->v2() == vertices_to_vertices.find(e->vn[1]->id)->second &&
                    dynamic_cast<XMLSubdomains::quad_type*>(&elements.element().at(elements_i))->v3() == vertices_to_vertices.find(e->vn[2]->id)->second &&
                    dynamic_cast<XMLSubdomains::quad_type*>(&elements.element().at(elements_i))->v4() == vertices_to_vertices.find(e->vn[3]->id)->second)
                    {
                      present = true;
                      break;
                    }

                if(!present)
                  elements.element().push_back(XMLSubdomains::quad_type(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->id, vertices_to_vertices.find(e->vn[3]->id)->second));
              }
            subdomain.elements()->i().push_back(e->id);
          }
        }

        // save boundary edge markers
        subdomain.boundary_edges().set(XMLSubdomains::subdomain::boundary_edges_type());
        bool has_inner_edges = false;
        for_all_base_elements(e, meshes[meshes_i])
        {
          for (unsigned i = 0; i < e->get_num_surf(); i++)
          {
            if(meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
            {
              if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
              {
                unsigned int edge_i = edges.edge().size();
                vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                edges.edge().push_back(XMLSubdomains::edge(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker.c_str(), edge_i));
              }
              subdomain.boundary_edges()->i().push_back(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)))->second);
            }
            else
              has_inner_edges = true;
          }
        }

        if(has_inner_edges)
        {
          subdomain.inner_edges().set(XMLSubdomains::subdomain::inner_edges_type());
          for_all_base_elements(e, meshes[meshes_i])
          for (unsigned i = 0; i < e->get_num_surf(); i++)
          {
            if(!meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
            {
              if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
              {
                unsigned int edge_i = edges.edge().size();
                vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                edges.edge().push_back(XMLSubdomains::edge(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker.c_str(), edge_i));
              }
              subdomain.inner_edges()->i().push_back(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)))->second);
            }
          }
        }

        // save curved edges
        for_all_base_elements(e, meshes[meshes_i])
          if(e->is_curved())
            for (unsigned i = 0; i < e->get_num_surf(); i++)
              if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
                if(vertices_to_curves.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_curves.end())
                {
                  if(e->cm->nurbs[i]->arc)
                    save_arc(meshes[meshes_i], vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, e->cm->nurbs[i], curves);
                  else
                    save_nurbs(meshes[meshes_i], vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, e->cm->nurbs[i], curves);
                  vertices_to_curves.insert(std::pair<std::pair<unsigned int, unsigned int>, bool>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), true));
                }

        // save refinements
        for(unsigned int refinement_i = 0; refinement_i < meshes[meshes_i]->refinements.size(); refinement_i++)
          refinements.refinement().push_back(XMLMesh::refinement(meshes[meshes_i]->refinements[refinement_i].first, meshes[meshes_i]->refinements[refinement_i].second));

        subdomain.refinements().set(refinements);
        subdomains.subdomain().push_back(subdomain);
      }

      XMLSubdomains::domain xmldomain(vertices, elements, edges, subdomains);
      xmldomain.curves().set(curves);

      std::string mesh_schema_location(H2D_XML_SCHEMAS_DIRECTORY);
      mesh_schema_location.append("/mesh_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_mesh("XMLMesh", mesh_schema_location);

      std::string domain_schema_location(H2D_XML_SCHEMAS_DIRECTORY);
      domain_schema_location.append("/subdomains_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_domain("XMLSubdomains", domain_schema_location);

      ::xml_schema::namespace_infomap namespace_info_map;
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("mesh", namespace_info_mesh));
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("domain", namespace_info_domain));

      std::ofstream out(filename);
      XMLSubdomains::domain_(out, xmldomain, namespace_info_map);
      out.close();

      return true;
    }

    bool MeshReaderH2DXML::load(std::auto_ptr<XMLMesh::mesh> & parsed_xml_mesh, Mesh *mesh, std::map<unsigned int, unsigned int>& vertex_is)
    {
      try
      {
        // Variables //
        unsigned int variables_count = parsed_xml_mesh->variables().present() ? parsed_xml_mesh->variables()->variable().size() : 0;
        std::map<std::string, double> variables;
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
          variables.insert(std::make_pair<std::string, double>(parsed_xml_mesh->variables()->variable().at(variables_i).name(), parsed_xml_mesh->variables()->variable().at(variables_i).value()));

        // Vertices //
        int vertices_count = parsed_xml_mesh->vertices().vertex().size();

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

          // variables matching.
          std::string x = parsed_xml_mesh->vertices().vertex().at(vertex_i).x();
          std::string y = parsed_xml_mesh->vertices().vertex().at(vertex_i).y();

          // insert into the map.
          vertex_is.insert(std::pair<unsigned int, unsigned int>(parsed_xml_mesh->vertices().vertex().at(vertex_i).i(), vertex_i));

          double x_value;
          double y_value;

          // variables lookup.
          bool x_found = false;
          bool y_found = false;
          if(variables.find(x) != variables.end())
          {
            x_value = variables.find(x)->second;
            x_found = true;
          }
          if(variables.find(y) != variables.end())
          {
            y_value = variables.find(y)->second;
            y_found = true;
          }

          // test of value if no variable found.
          if(!x_found)
            if(std::strtod(x.c_str(), NULL) != 0.0)
              x_value = std::strtod(x.c_str(), NULL);
            else
            {
              // This is a hard part, to find out if it is really zero.
              int dot_position = strchr(x.c_str(), '.') == NULL ? -1 : strchr(x.c_str(), '.') - x.c_str();
              for(int i = 0; i < dot_position; i++)
                if(strncmp(x.c_str() + i, "0", 1) != 0)
                  throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_i + 1);
              for(int i = dot_position + 1; i < x.length(); i++)
                if(strncmp(x.c_str() + i, "0", 1) != 0)
                  throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_i + 1);
              x_value = std::strtod(x.c_str(), NULL);
            }

            if(!y_found)
              if(std::strtod(y.c_str(), NULL) != 0.0)
                y_value = std::strtod(y.c_str(), NULL);
              else
              {
                // This is a hard part, to find out if it is really zero.
                int dot_position = strchr(y.c_str(), '.') == NULL ? -1 : strchr(y.c_str(), '.') - y.c_str();
                for(int i = 0; i < dot_position; i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_i + 1);
                for(int i = dot_position + 1; i < y.length(); i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_i + 1);
                y_value = std::strtod(y.c_str(), NULL);
              }

              // assignment.
              node->x = x_value;
              node->y = y_value;
        }
        mesh->ntopvert = vertices_count;

        // Elements //
        unsigned int element_count = parsed_xml_mesh->elements().element().size();
        mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

        Element* e;
        for (int element_i = 0; element_i < element_count; element_i++)
        {
          XMLMesh::mesh::elements_type::element_type* element = &parsed_xml_mesh->elements().element().at(element_i);

          // Trim whitespaces.
          unsigned int begin = element->marker().find_first_not_of(" \t\n");
          unsigned int end = element->marker().find_last_not_of(" \t\n");
          element->marker().erase(end + 1, element->marker().length());
          element->marker().erase(0, begin);

          mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, element->marker());

          if(dynamic_cast<XMLMesh::quad_type*>(element) != NULL)
            e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(element->marker()).marker,
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::quad_type*>(element)->v1())->second],
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::quad_type*>(element)->v2())->second],
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::quad_type*>(element)->v3())->second],
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::quad_type*>(element)->v4())->second],
            NULL);
          else
            e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(element->marker()).marker,
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::triangle_type*>(element)->v1())->second],
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::triangle_type*>(element)->v2())->second],
            &mesh->nodes[vertex_is.find(dynamic_cast<XMLMesh::triangle_type*>(element)->v3())->second],
            NULL);
        }

        // Boundaries //
        unsigned int edges_count = parsed_xml_mesh->edges().edge().size();

        Node* en;
        for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
        {
          int v1 = vertex_is.find(parsed_xml_mesh->edges().edge().at(edge_i).v1())->second;
          int v2 = vertex_is.find(parsed_xml_mesh->edges().edge().at(edge_i).v2())->second;

          en = mesh->peek_edge_node(v1, v2);
          if(en == NULL)
            throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist", edge_i, v1, v2);

          std::string edge_marker = parsed_xml_mesh->edges().edge().at(edge_i).marker();

          // Trim whitespaces.
          unsigned int begin = edge_marker.find_first_not_of(" \t\n");
          unsigned int end = edge_marker.find_last_not_of(" \t\n");
          edge_marker.erase(end + 1, edge_marker.length());
          edge_marker.erase(0, begin);

          // This functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, edge_marker);
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
            this->warn("Boundary edge node does not have a boundary marker");

        // Curves //
        // Arcs & NURBSs //
        unsigned int arc_count = parsed_xml_mesh->curves().present() ? parsed_xml_mesh->curves()->arc().size() : 0;
        unsigned int nurbs_count = parsed_xml_mesh->curves().present() ? parsed_xml_mesh->curves()->NURBS().size() : 0;

        for (unsigned int curves_i = 0; curves_i < arc_count + nurbs_count; curves_i++)
        {
          // load the control points, knot vector, etc.
          Node* en;
          int p1, p2;

          // first do arcs, then NURBSs.
          Nurbs* nurbs;
          if(curves_i < arc_count)
          {
            // read the end point indices
            p1 = vertex_is.find(parsed_xml_mesh->curves()->arc().at(curves_i).v1())->second;
            p2 = vertex_is.find(parsed_xml_mesh->curves()->arc().at(curves_i).v2())->second;

            nurbs = load_arc(mesh, parsed_xml_mesh, curves_i, &en, p1, p2);
          }
          else
          {
            // read the end point indices
            p1 = vertex_is.find(parsed_xml_mesh->curves()->NURBS().at(curves_i - arc_count).v1())->second;
            p2 = vertex_is.find(parsed_xml_mesh->curves()->NURBS().at(curves_i - arc_count).v2())->second;
            nurbs = load_nurbs(mesh, parsed_xml_mesh, curves_i - arc_count, &en, p1, p2);
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
            for (unsigned j = 0; j < e->get_num_surf(); j++)
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
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }

      return true;
    }

    bool MeshReaderH2DXML::load(std::auto_ptr<XMLSubdomains::domain> & parsed_xml_domain, Mesh *mesh, int* vertex_is, int* element_is, int* edge_is)
    {
      try
      {
        // Variables //
        unsigned int variables_count = parsed_xml_domain->variables().present() ? parsed_xml_domain->variables()->variable().size() : 0;
        std::map<std::string, double> variables;
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
          variables.insert(std::make_pair<std::string, double>(parsed_xml_domain->variables()->variable().at(variables_i).name(), parsed_xml_domain->variables()->variable().at(variables_i).value()));

        // Vertices //
        int vertices_count = parsed_xml_domain->vertices().vertex().size();

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

          // variables matching.
          std::string x = parsed_xml_domain->vertices().vertex().at(vertex_i).x();
          std::string y = parsed_xml_domain->vertices().vertex().at(vertex_i).y();

          if(parsed_xml_domain->vertices().vertex().at(vertex_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of vertex in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          // insert
          vertex_is[parsed_xml_domain->vertices().vertex().at(vertex_i).i()] = vertex_i;

          double x_value;
          double y_value;

          // variables lookup.
          bool x_found = false;
          bool y_found = false;
          if(variables.find(x) != variables.end())
          {
            x_value = variables.find(x)->second;
            x_found = true;
          }
          if(variables.find(y) != variables.end())
          {
            y_value = variables.find(y)->second;
            y_found = true;
          }

          // test of value if no variable found.
          if(!x_found)
            if(std::strtod(x.c_str(), NULL) != 0.0)
              x_value = std::strtod(x.c_str(), NULL);
            else
            {
              // This is a hard part, to find out if it is really zero.
              int dot_position = strchr(x.c_str(), '.') == NULL ? -1 : strchr(x.c_str(), '.') - x.c_str();
              for(int i = 0; i < dot_position; i++)
                if(strncmp(x.c_str() + i, "0", 1) != 0)
                  throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_i + 1);
              for(int i = dot_position + 1; i < x.length(); i++)
                if(strncmp(x.c_str() + i, "0", 1) != 0)
                  throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the x coordinate of vertex no. %i.", vertex_i + 1);
              x_value = std::strtod(x.c_str(), NULL);
            }

            if(!y_found)
              if(std::strtod(y.c_str(), NULL) != 0.0)
                y_value = std::strtod(y.c_str(), NULL);
              else
              {
                // This is a hard part, to find out if it is really zero.
                int dot_position = strchr(y.c_str(), '.') == NULL ? -1 : strchr(y.c_str(), '.') - y.c_str();
                for(int i = 0; i < dot_position; i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_i + 1);
                for(int i = dot_position + 1; i < y.length(); i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntay in the y coordinate of vertey no. %i.", vertex_i + 1);
                y_value = std::strtod(y.c_str(), NULL);
              }

              // assignment.
              node->x = x_value;
              node->y = y_value;
        }
        mesh->ntopvert = vertices_count;

        // Elements //
        unsigned int element_count = parsed_xml_domain->elements().element().size();
        mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

        Element* e;
        for (int element_i = 0; element_i < element_count; element_i++)
        {
          XMLSubdomains::domain::elements_type::element_type* element = &parsed_xml_domain->elements().element().at(element_i);

          // insert.
          if(parsed_xml_domain->elements().element().at(element_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of element in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          element_is[parsed_xml_domain->elements().element().at(element_i).i()] = element_i;

          // Trim whitespaces.
          unsigned int begin = element->marker().find_first_not_of(" \t\n");
          unsigned int end = element->marker().find_last_not_of(" \t\n");
          element->marker().erase(end + 1, element->marker().length());
          element->marker().erase(0, begin);

          mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, element->marker());

          if(dynamic_cast<XMLSubdomains::quad_type*>(element) != NULL)
            e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(element->marker()).marker,
            &mesh->nodes[dynamic_cast<XMLSubdomains::quad_type*>(element)->v1()],
            &mesh->nodes[dynamic_cast<XMLSubdomains::quad_type*>(element)->v2()],
            &mesh->nodes[dynamic_cast<XMLSubdomains::quad_type*>(element)->v3()],
            &mesh->nodes[dynamic_cast<XMLSubdomains::quad_type*>(element)->v4()],
            NULL);
          else
            e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(element->marker()).marker,
            &mesh->nodes[dynamic_cast<XMLSubdomains::triangle_type*>(element)->v1()],
            &mesh->nodes[dynamic_cast<XMLSubdomains::triangle_type*>(element)->v2()],
            &mesh->nodes[dynamic_cast<XMLSubdomains::triangle_type*>(element)->v3()],
            NULL);
        }

        // Boundaries //
        unsigned int edges_count = parsed_xml_domain->edges().edge().size();

        Node* en;
        for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
        {
          int v1 = parsed_xml_domain->edges().edge().at(edge_i).v1();
          int v2 = parsed_xml_domain->edges().edge().at(edge_i).v2();

          // insert
          if(parsed_xml_domain->edges().edge().at(edge_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of edge in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          edge_is[edge_i] = parsed_xml_domain->edges().edge().at(edge_i).i();

          en = mesh->peek_edge_node(v1, v2);
          if(en == NULL)
            throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist", edge_i, v1, v2);

          std::string edge_marker = parsed_xml_domain->edges().edge().at(edge_i).marker();

          // Trim whitespaces.
          unsigned int begin = edge_marker.find_first_not_of(" \t\n");
          unsigned int end = edge_marker.find_last_not_of(" \t\n");
          edge_marker.erase(end + 1, edge_marker.length());
          edge_marker.erase(0, begin);

          // This functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, edge_marker);
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
            this->warn("Boundary edge node does not have a boundary marker");

        // Curves //
        // Arcs & NURBSs //
        unsigned int arc_count = parsed_xml_domain->curves().present() ? parsed_xml_domain->curves()->arc().size() : 0;
        unsigned int nurbs_count = parsed_xml_domain->curves().present() ? parsed_xml_domain->curves()->NURBS().size() : 0;

        for (unsigned int curves_i = 0; curves_i < arc_count + nurbs_count; curves_i++)
        {
          // load the control points, knot vector, etc.
          Node* en;
          int p1, p2;

          // first do arcs, then NURBSs.
          Nurbs* nurbs;
          if(curves_i < arc_count)
          {
            // read the end point indices
            p1 = parsed_xml_domain->curves()->arc().at(curves_i).v1();
            p2 = parsed_xml_domain->curves()->arc().at(curves_i).v2();

            nurbs = load_arc(mesh, parsed_xml_domain, curves_i, &en, p1, p2);
          }
          else
          {
            // read the end point indices
            p1 = parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v1();
            p2 = parsed_xml_domain->curves()->NURBS().at(curves_i - arc_count).v2();
            nurbs = load_nurbs(mesh, parsed_xml_domain, curves_i - arc_count, &en, p1, p2);
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
            for (unsigned j = 0; j < e->get_num_surf(); j++)
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
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }

      return true;
    }

    template<typename T>
    Nurbs* MeshReaderH2DXML::load_arc(Mesh *mesh, std::auto_ptr<T>& parsed_xml_entity, int id, Node** en, int p1, int p2, bool skip_check)
    {
      Nurbs* nurbs = new Nurbs;
      nurbs->arc = true;

      *en = mesh->peek_edge_node(p1, p2);

      parsed_xml_entity.get();

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
      nurbs->angle = parsed_xml_entity->curves()->arc().at(id).angle();
      double a = (180.0 - nurbs->angle) / 180.0 * M_PI;

      // generate one inner control point
      double x = 1.0 / std::tan(a * 0.5);
      nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
      nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
      nurbs->pt[1][2] = Hermes::cos((M_PI - a) * 0.5);

      nurbs->ref = 0;

      return nurbs;
    }

    template<typename T>
    Nurbs* MeshReaderH2DXML::load_nurbs(Mesh *mesh, std::auto_ptr<T> & parsed_xml_entity, int id, Node** en, int p1, int p2, bool skip_check)
    {
      Nurbs* nurbs = new Nurbs;
      nurbs->arc = false;

      *en = mesh->peek_edge_node(p1, p2);

      if(*en == NULL)
      {
        if(!skip_check)
          throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: edge %d-%d does not exist.", id, p1, p2);
        else
          return NULL;
      }

      // degree of curved edge
      nurbs->degree = parsed_xml_entity->curves()->NURBS().at(id).degree();

      // get the number of control points
      int inner = parsed_xml_entity->curves()->NURBS().at(id).inner_point().size();

      nurbs->np = inner + 2;

      // edge endpoints are also control points, with weight 1.0
      nurbs->pt = new double3[nurbs->np];
      nurbs->pt[0][0] = mesh->nodes[p1].x;
      nurbs->pt[0][1] = mesh->nodes[p1].y;
      nurbs->pt[0][2] = 1.0;
      nurbs->pt[inner + 1][0] = mesh->nodes[p2].x;
      nurbs->pt[inner + 1][1] = mesh->nodes[p2].y;
      nurbs->pt[inner + 1][2] = 1.0;

      // read inner control points
      for (int i = 0; i < inner; i++)
      {
        nurbs->pt[i + 1][0] = parsed_xml_entity->curves()->NURBS().at(id).inner_point().at(i).x();
        nurbs->pt[i + 1][1] = parsed_xml_entity->curves()->NURBS().at(id).inner_point().at(i).y();
        nurbs->pt[i + 1][2] = parsed_xml_entity->curves()->NURBS().at(id).inner_point().at(i).weight();
      }

      // get the number of knot vector points
      inner = parsed_xml_entity->curves()->NURBS().at(id).knot().size();
      nurbs->nk = nurbs->degree + nurbs->np + 1;
      int outer = nurbs->nk - inner;
      if((outer & 1) == 1)
        throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: incorrect number of knot points.", id);

      // knot vector is completed by 0.0 on the left and by 1.0 on the right
      nurbs->kv = new double[nurbs->nk];

      for (int i = 0; i < outer/2; i++)
        nurbs->kv[i] = 0.0;

      if(inner > 0)
        for (int i = outer/2; i < inner + outer/2; i++)
          nurbs->kv[i] = parsed_xml_entity->curves()->NURBS().at(id).knot().at(i - (outer/2)).value();

      for (int i = outer/2 + inner; i < nurbs->nk; i++)
        nurbs->kv[i] = 1.0;

      nurbs->ref = 0;

      return nurbs;
    }

    void MeshReaderH2DXML::save_arc(Mesh *mesh, int p1, int p2, Nurbs* nurbs, XMLMesh::curves_type & curves)
    {
      curves.arc().push_back(XMLMesh::arc(p1, p2, nurbs->angle));
    }

    void MeshReaderH2DXML::save_nurbs(Mesh *mesh, int p1, int p2, Nurbs* nurbs, XMLMesh::curves_type & curves)
    {
      XMLMesh::NURBS nurbs_xml(p1, p2, nurbs->degree);

      int inner = nurbs->np - 2;
      int outer = nurbs->nk - inner;

      for (int i = 1; i < nurbs->np-1; i++)
        nurbs_xml.inner_point().push_back(XMLMesh::inner_point(nurbs->pt[i][0], nurbs->pt[i][1], nurbs->pt[i][2]));

      int max = nurbs->nk - (nurbs->degree + 1);
      for (int i = nurbs->degree + 1; i < max; i++)
        nurbs_xml.knot().push_back(XMLMesh::knot(nurbs->kv[i]));
    }
  }
}