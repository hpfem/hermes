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
#include "mesh_reader_h1d_xml.h"
#include <iostream>

using namespace std;

namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderH1DXML::MeshReaderH1DXML()
    {
    }

    MeshReaderH1DXML::~MeshReaderH1DXML()
    {
    }

    bool MeshReaderH1DXML::load(const char *filename, Mesh *mesh)
    {
      mesh->free();

      try
      {
        std::auto_ptr<XMLMesh1D::mesh> parsed_xml_mesh(XMLMesh1D::mesh_(filename));
        double a = parsed_xml_mesh->vertex().at(0).x();
        double b = parsed_xml_mesh->vertex().at(0).x();
        for(unsigned int vertex_i = 0; vertex_i < parsed_xml_mesh->vertex().size(); vertex_i++)
        {
          if(parsed_xml_mesh->vertex().at(vertex_i).x() > b)
            b = parsed_xml_mesh->vertex().at(vertex_i).x();
          if(parsed_xml_mesh->vertex().at(vertex_i).x() < a)
            a = parsed_xml_mesh->vertex().at(vertex_i).x();
        }
        
        // Vertices //
        int vertices_count = parsed_xml_mesh->vertex().size();

        // Initialize mesh.
        int size = HashTable::H2D_DEFAULT_HASH_SIZE;
        while (size < 8 * vertices_count)
          size *= 2;
        mesh->init(size);

        // Create top-level vertex nodes.
        for (int vertices_i = 0; vertices_i < 2 *vertices_count; vertices_i++)
        {
          Node* node = mesh->nodes.add();
          assert(node->id == vertices_i);
          node->ref = TOP_LEVEL_REF;
          node->type = HERMES_TYPE_VERTEX;
          node->bnd = 0;
          node->p1 = node->p2 = -1;          
          node->next_hash = NULL;

          // variables matching.
          double x = parsed_xml_mesh->vertex().at(vertices_i).x();
          double y;
          if(vertices_i < vertices_count)
            y = 0;
          else
            y = (b-a) / 100;

          // assignment.
          node->x = x;
          node->y = y;
        }
        mesh->ntopvert = 2 * vertices_count;

        // Elements //
        mesh->nbase = mesh->nactive = mesh->ninitial = vertices_count - 1;

        
        mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, "H1DMarker");

        Element* e;
        for (int element_i = 0; element_i < vertices_count; element_i++)
        {
          e = mesh->create_quad(mesh->element_markers_conversion.get_internal_marker("H1DMarker").marker, 
            &mesh->nodes[element_i], 
            &mesh->nodes[element_i + 1],
            &mesh->nodes[element_i + vertices_count + 1],
            &mesh->nodes[element_i + vertices_count],
            NULL);
        }

        // Boundaries //

        Node* en;
        int v1_1 = 0;
        int v2_1 = vertices_count;
        int v1_2 = vertices_count - 1;
        int v2_2 = 2 * vertices_count - 1;

        en = mesh->peek_edge_node(v1_1, v2_1);
        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, "H1DBndMarker");
        int marker = mesh->boundary_markers_conversion.get_internal_marker("H1DBndMarker").marker;
        en->marker = marker;
        en->bnd = 1;

        en = mesh->peek_edge_node(v1_2, v2_2);
        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, "H1DBndMarker");
        marker = mesh->boundary_markers_conversion.get_internal_marker("H1DBndMarker").marker;
        en->marker = marker;
        en->bnd = 1;

        mesh->nodes[v1_1].bnd = 1;
        mesh->nodes[v2_1].bnd = 1;

        mesh->nodes[v1_2].bnd = 1;
        mesh->nodes[v2_2].bnd = 1;

      }
      catch (const xml_schema::exception& e)
      {
        std::cerr << e << std::endl;
        std::exit(1);
      }
    }

    bool MeshReaderH1DXML::save(const char *filename, Mesh *mesh)
    {
      /// \todo Is this necessary? It is a valid H2D mesh afterall.
      return true;
    }
  }
}
