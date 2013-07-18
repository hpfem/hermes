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
#include "api2d.h"
#include "mesh_reader_h2d_xml.h"

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

    void MeshReaderH2DXML::load(const char *filename, MeshSharedPtr mesh)
    {
      try
      {
        ::xml_schema::flags parsing_flags = 0;
        if(!this->validate)
          parsing_flags = xml_schema::flags::dont_validate;

        // init
        std::auto_ptr<XMLMesh::mesh> parsed_xml_mesh(XMLMesh::mesh_(filename, parsing_flags));
        
        // load
        load(parsed_xml_mesh, mesh);
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    void MeshReaderH2DXML::load(std::auto_ptr<XMLMesh::mesh> & parsed_xml_mesh, MeshSharedPtr mesh)
    {
      if(!mesh)
        throw Exceptions::NullException(1);

      mesh->free();

      try
      {
        std::map<unsigned int, unsigned int> vertex_is;

        // load
        load(parsed_xml_mesh, mesh, vertex_is);
        
        // refinements.
        if(parsed_xml_mesh->refinements().present() && parsed_xml_mesh->refinements()->ref().size() > 0)
        {
          // perform initial refinements
          for (unsigned int i = 0; i < parsed_xml_mesh->refinements()->ref().size(); i++)
          {
            int element_id = parsed_xml_mesh->refinements()->ref().at(i).element_id();
            int refinement_type = parsed_xml_mesh->refinements()->ref().at(i).refinement_type();
            if(refinement_type == -1)
              mesh->unrefine_element_id(element_id);
            else
              mesh->refine_element_id(element_id, refinement_type);
          }
        }
        mesh->initial_single_check();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    void MeshReaderH2DXML::save(const char *filename, MeshSharedPtr mesh)
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

        vertices.v().push_back(std::auto_ptr<XMLMesh::v>(new XMLMesh::v(x_stream.str(), y_stream.str(), i)));
      }

      // save elements
      XMLMesh::elements_type elements;
      for (int i = 0; i < mesh->get_num_base_elements(); i++)
      {
        e = mesh->get_element_fast(i);
        if(e->used)
          if(e->is_triangle())
            elements.el().push_back(XMLMesh::t_t(e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str()));
          else
            elements.el().push_back(XMLMesh::q_t(e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->vn[3]->id));
      }
      // save boundary markers
      XMLMesh::edges_type edges;
      for_all_base_elements(e, mesh)
        for (unsigned i = 0; i < e->get_nvert(); i++)
          if(mesh->get_base_edge_node(e, i)->marker)
            edges.ed().push_back(XMLMesh::ed(e->vn[i]->id, e->vn[e->next_vert(i)]->id, mesh->boundary_markers_conversion.get_user_marker(mesh->get_base_edge_node(e, i)->marker).marker.c_str()));

      // save curved edges
      XMLMesh::curves_type curves;
      for_all_base_elements(e, mesh)
        if(e->is_curved())
          for (unsigned i = 0; i < e->get_nvert(); i++)
            if(e->cm->nurbs[i] != NULL && !is_twin_nurbs(e, i))
              if(e->cm->nurbs[i]->arc)
                save_arc(mesh, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i], curves);
              else
                save_nurbs(mesh, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->nurbs[i], curves);

      // save refinements
      XMLMesh::refinements_type refinements;
      for(unsigned int refinement_i = 0; refinement_i < mesh->refinements.size(); refinement_i++)
        refinements.ref().push_back(XMLMesh::ref(mesh->refinements[refinement_i].first, mesh->refinements[refinement_i].second));

      XMLMesh::mesh xmlmesh(vertices, elements, edges);
      xmlmesh.curves().set(curves);
      xmlmesh.refinements().set(refinements);

      std::string mesh_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
      mesh_schema_location.append("/mesh_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_mesh("XMLMesh", mesh_schema_location);

      ::xml_schema::namespace_infomap namespace_info_map;
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("mesh", namespace_info_mesh));

      std::ofstream out(filename);
      ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;
      XMLMesh::mesh_(out, xmlmesh, namespace_info_map, "UTF-8", parsing_flags);
      out.close();
    }

    void MeshReaderH2DXML::load(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
    {
      try
      {
        ::xml_schema::flags parsing_flags = 0;
        if(!this->validate)
          parsing_flags = xml_schema::flags::dont_validate;

        // init
        std::auto_ptr<XMLSubdomains::domain> parsed_xml_domain (XMLSubdomains::domain_(filename, parsing_flags));

        this->load(parsed_xml_domain, meshes);
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    void MeshReaderH2DXML::load(std::auto_ptr<XMLSubdomains::domain> & parsed_xml_domain, Hermes::vector<MeshSharedPtr > meshes)
    {
      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      {
        meshes.at(meshes_i)->free();
      }

      MeshSharedPtr global_mesh(new Mesh);

      try
      {
        std::map<int, int> vertex_is;
        std::map<int, int> element_is;
        std::map<int, int> edge_is;

        // load
        load(parsed_xml_domain, global_mesh, vertex_is, element_is, edge_is);
        
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
        unsigned int subdomains_count = parsed_xml_domain->subdomains().subdomain().size();
        if(subdomains_count != meshes.size())
          throw Hermes::Exceptions::MeshLoadFailureException("Number of subdomains( = %u) does not equal the number of provided meshes in the vector( = %u).", subdomains_count, meshes.size());

        for(unsigned int subdomains_i = 0; subdomains_i < subdomains_count; subdomains_i++)
        {
          for (int element_i = 0; element_i < parsed_xml_domain->elements().el().size(); element_i++)
          {
            XMLSubdomains::domain::elements_type::el_type* element = &parsed_xml_domain->elements().el().at(element_i);

            // Trim whitespaces.
            unsigned int begin = element->m().find_first_not_of(" \t\n");
            unsigned int end = element->m().find_last_not_of(" \t\n");
            element->m().erase(end + 1, element->m().length());
            element->m().erase(0, begin);

            meshes[subdomains_i]->element_markers_conversion.insert_marker(element->m());
          }
          for(unsigned int edge_i = 0; edge_i < parsed_xml_domain->edges().ed().size(); edge_i++)
          {
            XMLSubdomains::domain::edges_type::ed_type* edge = &parsed_xml_domain->edges().ed().at(edge_i);

            // Trim whitespaces.
            unsigned int begin = edge->m().find_first_not_of(" \t\n");
            unsigned int end = edge->m().find_last_not_of(" \t\n");
            edge->m().erase(end + 1, edge->m().length());
            edge->m().erase(0, begin);

            meshes[subdomains_i]->boundary_markers_conversion.insert_marker(edge->m());
          }
        }

        for(unsigned int subdomains_i = 0; subdomains_i < subdomains_count; subdomains_i++)
        {
          unsigned int vertex_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices()->i().size() : 0;
          unsigned int element_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements()->i().size() : 0;
          unsigned int boundary_edge_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().size() : 0;
          unsigned int inner_edge_number_count = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges().present() ? parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges()->i().size() : 0;

          // copy the whole mesh if the subdomain is the whole mesh.
          if(element_number_count == 0 || element_number_count == parsed_xml_domain->elements().el().size())
          {
            meshes[subdomains_i]->copy(global_mesh);
            // refinements.
            if(parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements().present() && parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().size() > 0)
            {
              // perform initial refinements
              for (unsigned int i = 0; i < parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().size(); i++)
              {
                int element_id = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().at(i).element_id();
                int refinement_type = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().at(i).refinement_type();
                if(refinement_type == -1)
                  meshes[subdomains_i]->unrefine_element_id(element_id);
                else
                  meshes[subdomains_i]->refine_element_id(element_id, refinement_type);
              }
            }
          }
          else
          {
            // Variables //
            unsigned int variables_count = parsed_xml_domain->variables().present() ? parsed_xml_domain->variables()->var().size() : 0;

            std::map<std::string, double> variables;
            for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
#ifdef _MSC_VER
              variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_domain->variables()->var().at(variables_i).name(), (double&&)parsed_xml_domain->variables()->var().at(variables_i).value()));
#else
              variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_domain->variables()->var().at(variables_i).name(), parsed_xml_domain->variables()->var().at(variables_i).value()));
#endif
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
              vertex_number_count = parsed_xml_domain->vertices().v().size();
            for (unsigned int vertex_numbers_i = 0; vertex_numbers_i < vertex_number_count; vertex_numbers_i++)
            {
              unsigned int vertex_number;
              if(vertex_number_count == parsed_xml_domain->vertices().v().size())
                vertex_number = vertex_is[vertex_numbers_i];
              else
              {
                vertex_number =  parsed_xml_domain->subdomains().subdomain().at(subdomains_i).vertices()->i().at(vertex_numbers_i);
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

              // variables matching.
              std::string x = parsed_xml_domain->vertices().v().at(vertex_number).x();
              std::string y = parsed_xml_domain->vertices().v().at(vertex_number).y();
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
                        throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_number + 1);
                    for(int i = dot_position + 1; i < y.length(); i++)
                      if(strncmp(y.c_str() + i, "0", 1) != 0)
                        throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_number + 1);
                    y_value = std::strtod(y.c_str(), NULL);
                  }

                  // assignment.
                  node->x = x_value;
                  node->y = y_value;
            }
            meshes[subdomains_i]->ntopvert = vertex_number_count;

            // Element numbers //
            unsigned int element_count = parsed_xml_domain->elements().el().size();
            meshes[subdomains_i]->nbase = element_count;
            meshes[subdomains_i]->nactive = meshes[subdomains_i]->ninitial = element_number_count;

            Element* e;
            int* elements_existing = new int[element_count];
            for(int i = 0; i < element_count; i++)
              elements_existing[i] = -1;
            for (int element_number_i = 0; element_number_i < element_number_count; element_number_i++)
            {
              int elementI = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements()->i().at(element_number_i);
              if(elementI > max_element_i)
                  throw Exceptions::MeshLoadFailureException("Wrong element number:%i in subdomain %u.", elementI, subdomains_i);

              elements_existing[element_is[parsed_xml_domain->subdomains().subdomain().at(subdomains_i).elements()->i().at(element_number_i)]] = elementI;
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

              XMLSubdomains::domain::elements_type::el_type* element = NULL;
              for(int searched_element_i = 0; searched_element_i < element_count; searched_element_i++)
              {
                element = &parsed_xml_domain->elements().el().at(searched_element_i);
                if(element->i() == elements_existing[element_i])
                  break;
                else
                  element = NULL;
              }
              if(element == NULL)
                throw Exceptions::MeshLoadFailureException("Element number wrong in the mesh file.");

              XMLSubdomains::q_t* el_q = dynamic_cast<XMLSubdomains::q_t*>(element);
              XMLSubdomains::t_t* el_t = dynamic_cast<XMLSubdomains::t_t*>(element);
              if(el_q != NULL)
                e = meshes[subdomains_i]->create_quad(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element->m()).marker,
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_q->v1())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_q->v2())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_q->v3())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_q->v4())->second],
                NULL, element_i);
              if(el_t != NULL)
                e = meshes[subdomains_i]->create_triangle(meshes[subdomains_i]->element_markers_conversion.get_internal_marker(element->m()).marker,
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_t->v1())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_t->v2())->second],
                &meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(el_t->v3())->second],
                NULL, element_i);
            }

            // Boundary Edge numbers //
            if(boundary_edge_number_count == 0)
              boundary_edge_number_count = parsed_xml_domain->edges().ed().size();

            for (int boundary_edge_number_i = 0; boundary_edge_number_i < boundary_edge_number_count; boundary_edge_number_i++)
            {
              XMLSubdomains::domain::edges_type::ed_type* edge = NULL;
              int domain_edge_count = parsed_xml_domain->edges().ed().size();
              for(unsigned int to_find_i = 0; to_find_i < domain_edge_count; to_find_i++)
              {
                if(boundary_edge_number_count != domain_edge_count)
                {
                  if(parsed_xml_domain->edges().ed().at(to_find_i).i() == parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().at(boundary_edge_number_i))
                  {
                    edge = &parsed_xml_domain->edges().ed().at(to_find_i);
                    break;
                  }
                }
                else
                {
                  if(parsed_xml_domain->edges().ed().at(to_find_i).i() == edge_is[boundary_edge_number_i])
                  {
                    edge = &parsed_xml_domain->edges().ed().at(to_find_i);
                    break;
                  }
                }
              }

              if(edge == NULL)
                  throw Exceptions::MeshLoadFailureException("Wrong boundary-edge number:%i in subdomain %u.", parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().at(boundary_edge_number_i), subdomains_i);

              Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge->v1())->second, vertex_vertex_numbers.find(edge->v2())->second);
              if(en == NULL)
                throw Hermes::Exceptions::MeshLoadFailureException("Boundary data error (edge %i does not exist).", boundary_edge_number_i);

              en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge->m()).marker;

              meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge->v1())->second].bnd = 1;
              meshes[subdomains_i]->nodes[vertex_vertex_numbers.find(edge->v2())->second].bnd = 1;
              en->bnd = 1;
            }

            // Inner Edge numbers //
            for (int inner_edge_number_i = 0; inner_edge_number_i < inner_edge_number_count; inner_edge_number_i++)
            {
              XMLSubdomains::domain::edges_type::ed_type* edge = NULL;

              for(unsigned int to_find_i = 0; to_find_i < parsed_xml_domain->edges().ed().size(); to_find_i++)
              {
                if(parsed_xml_domain->edges().ed().at(to_find_i).i() == parsed_xml_domain->subdomains().subdomain().at(subdomains_i).inner_edges()->i().at(inner_edge_number_i))
                {
                  edge = &parsed_xml_domain->edges().ed().at(to_find_i);
                  break;
                }
              }

              if(edge == NULL)
                  throw Exceptions::MeshLoadFailureException("Wrong inner-edge number:%i in subdomain %u.", parsed_xml_domain->subdomains().subdomain().at(subdomains_i).boundary_edges()->i().at(inner_edge_number_i), subdomains_i);

              Node* en = meshes[subdomains_i]->peek_edge_node(vertex_vertex_numbers.find(edge->v1())->second, vertex_vertex_numbers.find(edge->v2())->second);
              if(en == NULL)
                throw Hermes::Exceptions::MeshLoadFailureException("Inner data error (edge %i does not exist).", inner_edge_number_i);

              en->marker = meshes[subdomains_i]->boundary_markers_conversion.get_internal_marker(edge->m()).marker;
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
            for_all_elements(e, meshes[subdomains_i])
              if(e->cm != NULL)
                e->cm->update_refmap_coeffs(e);

            // refinements.
            if(parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements().present() && parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().size() > 0)
            {
              // perform initial refinements
              for (unsigned int i = 0; i < parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().size(); i++)
              {
                int element_id = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().at(i).element_id();
                int refinement_type = parsed_xml_domain->subdomains().subdomain().at(subdomains_i).refinements()->ref().at(i).refinement_type();
                if(refinement_type == -1)
                  meshes[subdomains_i]->unrefine_element_id(element_id);
                else
                  meshes[subdomains_i]->refine_element_id(element_id, refinement_type);
              }
            }

            delete [] elements_existing;
          }
          meshes[subdomains_i]->seq = g_mesh_seq++;
          meshes[subdomains_i]->initial_single_check();
        }
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    static bool elementCompare (XMLSubdomains::el_t* el_i, XMLSubdomains::el_t* el_j) { return ( el_i->i() < el_j->i() ); }

    void MeshReaderH2DXML::save(const char *filename, Hermes::vector<MeshSharedPtr > meshes)
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
      Hermes::vector<XMLSubdomains::el_t*> elements;
      // Global boudnary edges list.
      XMLSubdomains::edges_type edges;
      // Global curves list.
      XMLMesh::curves_type curves;

      bool* baseElementsSaved = new bool[meshes[0]->get_num_base_elements()];
      memset(baseElementsSaved, 0, sizeof(bool) * meshes[0]->get_num_base_elements());

      // Subdomains.
      XMLSubdomains::subdomains subdomains;

      for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      {
        bool hasAllElements = (meshes[meshes_i]->get_num_used_base_elements() == meshes[meshes_i]->get_num_base_elements());

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
          std::map<std::pair<double, double>, unsigned int>::iterator it = points_to_vertices.find(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y));
          if(it != points_to_vertices.end())
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, it->second));
          else
          {
            int new_i = points_to_vertices.size();
            vertices_to_vertices.insert(std::pair<unsigned int, unsigned int>(i, new_i));
            points_to_vertices.insert(std::pair<std::pair<double, double>, unsigned int>(std::pair<double, double>(meshes[meshes_i]->nodes[i].x, meshes[meshes_i]->nodes[i].y), points_to_vertices.size()));
            std::ostringstream x_stream;
            x_stream << meshes[meshes_i]->nodes[i].x;

            std::ostringstream y_stream;
            y_stream << meshes[meshes_i]->nodes[i].y;

            vertices.v().push_back(std::auto_ptr<XMLMesh::v>(new XMLMesh::v(x_stream.str(), y_stream.str(), new_i)));
          }
          if(!hasAllElements)
            subdomain.vertices()->i().push_back(vertices_to_vertices.find(i)->second);
        }

        // save elements
        subdomain.elements().set(XMLSubdomains::subdomain::elements_type());
        for (int i = 0; i < meshes[meshes_i]->get_num_base_elements(); i++)
        {
          e = &(meshes[meshes_i]->elements[i]);
          if(!e->used)
            continue;
          if(!baseElementsSaved[e->id])
          {
            if(e->nvert == 3)
            {
              elements.push_back(new XMLSubdomains::t_t(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->id));
            }
            else
            {
              elements.push_back(new XMLSubdomains::q_t(vertices_to_vertices.find(e->vn[0]->id)->second, vertices_to_vertices.find(e->vn[1]->id)->second, vertices_to_vertices.find(e->vn[2]->id)->second, meshes[meshes_i]->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str(), e->id, vertices_to_vertices.find(e->vn[3]->id)->second));
            }
            baseElementsSaved[e->id] = true;
          }

          if(!hasAllElements)
            subdomain.elements()->i().push_back(e->id);
        }

        // save boundary edge markers
        subdomain.boundary_edges().set(XMLSubdomains::subdomain::boundary_edges_type());
        bool has_inner_edges = false;
        for_all_base_elements(e, meshes[meshes_i])
        {
          for (unsigned i = 0; i < e->get_nvert(); i++)
          {
            if(meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
            {
              if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
              {
                unsigned int edge_i = edges.ed().size();
                vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                edges.ed().push_back(XMLSubdomains::ed(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker.c_str(), edge_i));
              }
              if(!hasAllElements)
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
            for (unsigned i = 0; i < e->get_nvert(); i++)
            {
              if(!meshes[meshes_i]->get_base_edge_node(e, i)->bnd)
              {
                if(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second))) == vertices_to_boundaries.end())
                {
                  unsigned int edge_i = edges.ed().size();
                  vertices_to_boundaries.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), edge_i));
                  edges.ed().push_back(XMLSubdomains::ed(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, meshes[meshes_i]->boundary_markers_conversion.get_user_marker(meshes[meshes_i]->get_base_edge_node(e, i)->marker).marker.c_str(), edge_i));
                }
              if(!hasAllElements)
                subdomain.inner_edges()->i().push_back(vertices_to_boundaries.find(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)))->second);
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
                    save_arc(meshes[meshes_i], vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, e->cm->nurbs[i], curves);
                  else
                    save_nurbs(meshes[meshes_i], vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second, e->cm->nurbs[i], curves);
                  vertices_to_curves.insert(std::pair<std::pair<unsigned int, unsigned int>, bool>(std::pair<unsigned int, unsigned int>(std::min(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second), std::max(vertices_to_vertices.find(e->vn[i]->id)->second, vertices_to_vertices.find(e->vn[e->next_vert(i)]->id)->second)), true));
                }

        // save refinements
        for(unsigned int refinement_i = 0; refinement_i < meshes[meshes_i]->refinements.size(); refinement_i++)
          refinements.ref().push_back(XMLMesh::ref(meshes[meshes_i]->refinements[refinement_i].first, meshes[meshes_i]->refinements[refinement_i].second));

        subdomain.refinements().set(refinements);
        subdomains.subdomain().push_back(subdomain);
      }

      delete [] baseElementsSaved;

      std::sort (elements.begin(), elements.end(), elementCompare);

      XMLSubdomains::elements_type elementsToPass;
      for(int i = 0; i < elements.size(); i++)
        elementsToPass.el().push_back(*elements[i]);

      for(int i = 0; i < elements.size(); i++)
        delete elements[i];
      XMLSubdomains::domain xmldomain(vertices, elementsToPass, edges, subdomains);
      xmldomain.curves().set(curves);

      std::string mesh_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
      mesh_schema_location.append("/mesh_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_mesh("XMLMesh", mesh_schema_location);

      std::string domain_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
      domain_schema_location.append("/subdomains_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_domain("XMLSubdomains", domain_schema_location);

      ::xml_schema::namespace_infomap namespace_info_map;
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("mesh", namespace_info_mesh));
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("domain", namespace_info_domain));

      std::ofstream out(filename);
      ::xml_schema::flags parsing_flags = ::xml_schema::flags::base;
      XMLSubdomains::domain_(out, xmldomain, namespace_info_map, "UTF-8", parsing_flags);
      out.close();
    }

    void MeshReaderH2DXML::load(std::auto_ptr<XMLMesh::mesh> & parsed_xml_mesh, MeshSharedPtr mesh, std::map<unsigned int, unsigned int>& vertex_is)
    {
      try
      {
        // Variables //
        unsigned int variables_count = parsed_xml_mesh->variables().present() ? parsed_xml_mesh->variables()->var().size() : 0;
        std::map<std::string, double> variables;
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
#ifdef _MSC_VER
          variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_mesh->variables()->var().at(variables_i).name(), (double&&)parsed_xml_mesh->variables()->var().at(variables_i).value()));
#else
          variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_mesh->variables()->var().at(variables_i).name(), parsed_xml_mesh->variables()->var().at(variables_i).value()));
#endif

        // Vertices //
        int vertices_count = parsed_xml_mesh->vertices().v().size();

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
          std::string x = parsed_xml_mesh->vertices().v().at(vertex_i).x();
          std::string y = parsed_xml_mesh->vertices().v().at(vertex_i).y();

          // insert into the map.
          vertex_is.insert(std::pair<unsigned int, unsigned int>(parsed_xml_mesh->vertices().v().at(vertex_i).i(), vertex_i));

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
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_i + 1);
                for(int i = dot_position + 1; i < y.length(); i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_i + 1);
                y_value = std::strtod(y.c_str(), NULL);
              }

              // assignment.
              node->x = x_value;
              node->y = y_value;
        }
        mesh->ntopvert = vertices_count;

        // Elements //
        unsigned int element_count = parsed_xml_mesh->elements().el().size();
        mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

        Element* e;
        for (int element_i = 0; element_i < element_count; element_i++)
        {
          XMLMesh::mesh::elements_type::el_type* element = &parsed_xml_mesh->elements().el().at(element_i);

          // Trim whitespaces.
          unsigned int begin = element->m().find_first_not_of(" \t\n");
          unsigned int end = element->m().find_last_not_of(" \t\n");
          element->m().erase(end + 1, element->m().length());
          element->m().erase(0, begin);

          mesh->element_markers_conversion.insert_marker(element->m());

          XMLMesh::q_t* el_q = dynamic_cast<XMLMesh::q_t*>(element);
          XMLMesh::t_t* el_t = dynamic_cast<XMLMesh::t_t*>(element);
          if(el_q != NULL)
            e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(element->m()).marker,
            &mesh->nodes[vertex_is.find(el_q->v1())->second],
            &mesh->nodes[vertex_is.find(el_q->v2())->second],
            &mesh->nodes[vertex_is.find(el_q->v3())->second],
            &mesh->nodes[vertex_is.find(el_q->v4())->second],
            NULL);
          if(el_t != NULL)
            e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(element->m()).marker,
            &mesh->nodes[vertex_is.find(el_t->v1())->second],
            &mesh->nodes[vertex_is.find(el_t->v2())->second],
            &mesh->nodes[vertex_is.find(el_t->v3())->second],
            NULL);
        }

        // Boundaries //
        unsigned int edges_count = parsed_xml_mesh->edges().ed().size();

        Node* en;
        for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
        {
          int v1 = vertex_is.find(parsed_xml_mesh->edges().ed().at(edge_i).v1())->second;
          int v2 = vertex_is.find(parsed_xml_mesh->edges().ed().at(edge_i).v2())->second;

          en = mesh->peek_edge_node(v1, v2);
          if(en == NULL)
            throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist.", edge_i, v1, v2);

          std::string edge_marker = parsed_xml_mesh->edges().ed().at(edge_i).m();

          // Trim whitespaces.
          unsigned int begin = edge_marker.find_first_not_of(" \t\n");
          unsigned int end = edge_marker.find_last_not_of(" \t\n");
          edge_marker.erase(end + 1, edge_marker.length());
          edge_marker.erase(0, begin);

          // This functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->boundary_markers_conversion.insert_marker(edge_marker);

          en->marker = mesh->boundary_markers_conversion.get_internal_marker(edge_marker).marker;

          // This is extremely important, as in DG, it is assumed that negative boundary markers are reserved
          // for the inner edges.
          if(en->marker > 0)
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
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    void MeshReaderH2DXML::load(std::auto_ptr<XMLSubdomains::domain> & parsed_xml_domain, MeshSharedPtr mesh, std::map<int, int>& vertex_is, std::map<int, int>& element_is, std::map<int, int>& edge_is)
    {
      try
      {
        // Variables //
        unsigned int variables_count = parsed_xml_domain->variables().present() ? parsed_xml_domain->variables()->var().size() : 0;
        std::map<std::string, double> variables;
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
#ifdef _MSC_VER
          variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_domain->variables()->var().at(variables_i).name(), (double&&)parsed_xml_domain->variables()->var().at(variables_i).value()));
#else
          variables.insert(std::make_pair<std::string, double>((std::string)parsed_xml_domain->variables()->var().at(variables_i).name(), parsed_xml_domain->variables()->var().at(variables_i).value()));
#endif


        // Vertices //
        int vertices_count = parsed_xml_domain->vertices().v().size();

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
          std::string x = parsed_xml_domain->vertices().v().at(vertex_i).x();
          std::string y = parsed_xml_domain->vertices().v().at(vertex_i).y();

          if(parsed_xml_domain->vertices().v().at(vertex_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of vertex in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          // insert
          vertex_is.insert(std::pair<int,int>(parsed_xml_domain->vertices().v().at(vertex_i).i(), vertex_i));

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
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_i + 1);
                for(int i = dot_position + 1; i < y.length(); i++)
                  if(strncmp(y.c_str() + i, "0", 1) != 0)
                    throw Hermes::Exceptions::MeshLoadFailureException("Wrong syntax in the y coordinate of vertex no. %i.", vertex_i + 1);
                y_value = std::strtod(y.c_str(), NULL);
              }

              // assignment.
              node->x = x_value;
              node->y = y_value;
        }
        mesh->ntopvert = vertices_count;

        // Elements //
        unsigned int element_count = parsed_xml_domain->elements().el().size();
        mesh->nbase = mesh->nactive = mesh->ninitial = element_count;

        Element* e;
        for (int element_i = 0; element_i < element_count; element_i++)
        {
          XMLSubdomains::domain::elements_type::el_type* element = &parsed_xml_domain->elements().el().at(element_i);

          // insert.
          if(parsed_xml_domain->elements().el().at(element_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of element in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          element_is.insert(std::pair<int,int>(parsed_xml_domain->elements().el().at(element_i).i(), element_i));

          // Trim whitespaces.
          unsigned int begin = element->m().find_first_not_of(" \t\n");
          unsigned int end = element->m().find_last_not_of(" \t\n");
          element->m().erase(end + 1, element->m().length());
          element->m().erase(0, begin);

          mesh->element_markers_conversion.insert_marker(element->m());

          XMLSubdomains::q_t* el_q = dynamic_cast<XMLSubdomains::q_t*>(element);
          XMLSubdomains::t_t* el_t = dynamic_cast<XMLSubdomains::t_t*>(element);
          if(el_q != NULL)
            e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker(element->m()).marker,
            &mesh->nodes[el_q->v1()],
            &mesh->nodes[el_q->v2()],
            &mesh->nodes[el_q->v3()],
            &mesh->nodes[el_q->v4()],
            NULL);
          if(el_t != NULL)
            e = mesh->create_triangle(mesh->element_markers_conversion.get_internal_marker(element->m()).marker,
            &mesh->nodes[el_t->v1()],
            &mesh->nodes[el_t->v2()],
            &mesh->nodes[el_t->v3()],
            NULL);
        }

        // Boundaries //
        unsigned int edges_count = parsed_xml_domain->edges().ed().size();

        Node* en;
        for (unsigned int edge_i = 0; edge_i < edges_count; edge_i++)
        {
          int v1 = parsed_xml_domain->edges().ed().at(edge_i).v1();
          int v2 = parsed_xml_domain->edges().ed().at(edge_i).v2();

          // insert
          if(parsed_xml_domain->edges().ed().at(edge_i).i() > H2D_MAX_NODE_ID - 1)
            throw Exceptions::MeshLoadFailureException("The index 'i' of edge in the mesh file must be lower than %i.", H2D_MAX_NODE_ID);

          edge_is.insert(std::pair<int, int>(edge_i, parsed_xml_domain->edges().ed().at(edge_i).i()));

          en = mesh->peek_edge_node(v1, v2);
          if(en == NULL)
            throw Hermes::Exceptions::MeshLoadFailureException("Boundary data #%d: edge %d-%d does not exist.", edge_i, v1, v2);

          std::string edge_marker = parsed_xml_domain->edges().ed().at(edge_i).m();

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
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::MeshLoadFailureException(e.what());
      }
    }

    template<typename T>
    Nurbs* MeshReaderH2DXML::load_arc(MeshSharedPtr mesh, std::auto_ptr<T>& parsed_xml_entity, int id, Node** en, int p1, int p2, bool skip_check)
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
    Nurbs* MeshReaderH2DXML::load_nurbs(MeshSharedPtr mesh, std::auto_ptr<T> & parsed_xml_entity, int id, Node** en, int p1, int p2, bool skip_check)
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
      nurbs->degree = parsed_xml_entity->curves()->NURBS().at(id).deg();

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

    void MeshReaderH2DXML::save_arc(MeshSharedPtr mesh, int p1, int p2, Nurbs* nurbs, XMLMesh::curves_type & curves)
    {
      curves.arc().push_back(XMLMesh::arc(p1, p2, nurbs->angle));
    }

    void MeshReaderH2DXML::save_nurbs(MeshSharedPtr mesh, int p1, int p2, Nurbs* nurbs, XMLMesh::curves_type & curves)
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
