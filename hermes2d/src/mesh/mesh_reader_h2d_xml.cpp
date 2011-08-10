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
      try
      {
        std::auto_ptr<mesh_h2d> parsed_h2d_xml_mesh (mesh_h2d_(filename));

        // Variables //.
        unsigned int variables_count = parsed_h2d_xml_mesh->variables().variable().size();
        double *variables = new double[variables_count];
        for (unsigned int variables_i = 0; variables_i < variables_count; variables_i++)
          variables[variables_i] = parsed_h2d_xml_mesh->variables().variable().at(variables_i).value();

        // Vertices //.
        unsigned int vertices_count = parsed_h2d_xml_mesh->vertices().vertex().size();

        // Initialize mesh.
        unsigned int size = HashTable::H2D_DEFAULT_HASH_SIZE;
        while (size < 8 * vertices_count)
          size *= 2;
        mesh->init(size);

        // Create top-level vertex nodes.
        for (unsigned int vertices_i = 0; vertices_i < vertices_count; vertices_i++)
        {
          Node* node = mesh->nodes.add();
          assert(node->id == vertices_i);
          node->ref = TOP_LEVEL_REF;
          node->type = HERMES_TYPE_VERTEX;
          node->bnd = 0;
          node->p1 = node->p2 = -1;
          node->next_hash = NULL;
          node->x = parsed_h2d_xml_mesh->vertices().vertex().at(vertices_i).x();
          node->y = parsed_h2d_xml_mesh->vertices().vertex().at(vertices_i).y();
        }
        mesh->ntopvert = vertices_count;

        // Elements.//
        unsigned int triangles_count = parsed_h2d_xml_mesh->elements().triangle().size();
        unsigned int quads_count = parsed_h2d_xml_mesh->elements().quad().size();

        // Create elements.
        mesh->nactive = 0;
        for (unsigned int triangles_i = 0; triangles_i < triangles_count; triangles_i++)
        {
          // read and check vertex indices
          int* idx = new int[3];
          std::string el_marker;
          if (parsed_h2d_xml_mesh->elements().triangle().at(triangles_i).empty()) { 
            mesh->elements.skip_slot();
            continue;
          }

          idx[0] = parsed_h2d_xml_mesh->elements().triangle().at(triangles_i).vertices().vertex().at(0);
          idx[1] = parsed_h2d_xml_mesh->elements().triangle().at(triangles_i).vertices().vertex().at(1);
          idx[2] = parsed_h2d_xml_mesh->elements().triangle().at(triangles_i).vertices().vertex().at(2);

          el_marker = parsed_h2d_xml_mesh->elements().triangle().at(triangles_i).marker();
          
          // Trim whitespaces.
          unsigned int begin = el_marker.find_first_not_of(" \t\n");
          unsigned int end = el_marker.find_last_not_of(" \t\n");
          el_marker.erase(end + 1, el_marker.length());
          el_marker.erase(0, begin);

          for (unsigned int vertex_i = 0; vertex_i < 3; vertex_i++)
            if (idx[vertex_i] < 0 || idx[vertex_i] >= mesh->ntopvert)
              error("File %s: error creating triangle #%d: vertex #%d does not exist.", filename, triangles_i, idx[vertex_i]);

          Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]];

          int marker;

          // These functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, el_marker);
          marker = mesh->element_markers_conversion.get_internal_marker(el_marker);

          Mesh::check_triangle(triangles_i, v0, v1, v2);
          mesh->create_triangle(marker, v0, v1, v2, NULL);

          mesh->nactive++;

          delete [] idx;
        }

        for (unsigned int quads_i = 0; quads_i < quads_count; quads_i++)
        {
          // read and check vertex indices
          int* idx = new int[4];
          std::string el_marker;
          if (parsed_h2d_xml_mesh->elements().quad().at(quads_i).empty()) { 
            mesh->elements.skip_slot();
            continue;
          }

          idx[0] = parsed_h2d_xml_mesh->elements().quad().at(quads_i).vertices().vertex().at(0);
          idx[1] = parsed_h2d_xml_mesh->elements().quad().at(quads_i).vertices().vertex().at(1);
          idx[2] = parsed_h2d_xml_mesh->elements().quad().at(quads_i).vertices().vertex().at(2);
          idx[3] = parsed_h2d_xml_mesh->elements().quad().at(quads_i).vertices().vertex().at(3);

          el_marker = parsed_h2d_xml_mesh->elements().quad().at(quads_i).marker();
          
          // Trim whitespaces.
          unsigned int begin = el_marker.find_first_not_of(" \t\n");
          unsigned int end = el_marker.find_last_not_of(" \t\n");
          el_marker.erase(end + 1, el_marker.length());
          el_marker.erase(0, begin);

          for (unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
            if (idx[vertex_i] < 0 || idx[vertex_i] >= mesh->ntopvert)
              error("File %s: error creating quad #%d: vertex #%d does not exist.", filename, quads_i, idx[vertex_i]);

          Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]], *v3 = &mesh->nodes[idx[3]];

          int marker;

          // This functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, el_marker);
          marker = mesh->element_markers_conversion.get_internal_marker(el_marker);

          Mesh::check_quad(quads_i, v0, v1, v2, v3);
          mesh->create_quad(marker, v0, v1, v2, v3, NULL);

          mesh->nactive++;

          delete [] idx;
        }

        mesh->nbase = triangles_count + quads_count;

        // Boundaries. //
        unsigned int boundaries_count = parsed_h2d_xml_mesh->boundaries().boundary_edge().size();

        Node* en;
        for (unsigned int boundaries_i = 0; boundaries_i < boundaries_count; boundaries_i++)
        {
          int v1 = parsed_h2d_xml_mesh->boundaries().boundary_edge().at(boundaries_i).vertex().at(0);
          int v2 = parsed_h2d_xml_mesh->boundaries().boundary_edge().at(boundaries_i).vertex().at(1);

          en = mesh->peek_edge_node(v1, v2);
          if (en == NULL)
            error("File %s: boundary data #%d: edge %d-%d does not exist", filename, boundaries_i, v1, v2);

          std::string bnd_marker = parsed_h2d_xml_mesh->boundaries().boundary_edge().at(boundaries_i).marker();
          
          // Trim whitespaces.
          unsigned int begin = bnd_marker.find_first_not_of(" \t\n");
          unsigned int end = bnd_marker.find_last_not_of(" \t\n");
          bnd_marker.erase(end + 1, bnd_marker.length());
          bnd_marker.erase(0, begin);

          // This functions check if the user-supplied marker on this element has been
          // already used, and if not, inserts it in the appropriate structure.
          mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, bnd_marker);
          int marker = mesh->boundary_markers_conversion.get_internal_marker(bnd_marker);

          en->marker = marker;

          // This is extremely important, as in DG, it is assumed that negative boundary markers are reserved
          // for the inner edges.
          if (marker > 0)
          {
            mesh->nodes[v1].bnd = 1;
            mesh->nodes[v2].bnd = 1;
            en->bnd = 1;
          }
        }

        // check that all boundary edges have a marker assigned
        for_all_edge_nodes(en, mesh)
          if (en->ref < 2 && en->marker == 0)
            warn("Boundary edge node does not have a boundary marker");

        // Curves. //
        unsigned int curves_count = parsed_h2d_xml_mesh->curves().curve().size();

        for (unsigned int curves_i = 0; curves_i < curves_count; curves_i++)
        {
          // load the control points, knot vector, etc.
          Node* en;
          int p1, p2;

          Nurbs* nurbs = load_nurbs(mesh, parsed_h2d_xml_mesh, curves_i, &en, p1, p2);

          // assign the nurbs to the elements sharing the edge node
          for (unsigned int node_i = 0; node_i < 2; node_i++)
          {
            Element* e = en->elem[node_i];
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

        // update refmap coeffs of curvilinear elements
        Element* e;
        for_all_elements(e, mesh)
          if (e->cm != NULL)
            e->cm->update_refmap_coeffs(e);

      }
      catch (const xml_schema::exception& e)
      {
        cerr << e << endl;
        return 1;
      }

      return true;
    }

    Nurbs* MeshReaderH2DXML::load_nurbs(Mesh *mesh, std::auto_ptr<mesh_h2d>& m, int id, Node** en, int &p1, int &p2)
    {
      Nurbs* nurbs = new Nurbs;

      /// \todo We can only parse arcs now.
      nurbs->arc = true;

      // read the end point indices
      p1 = m->curves().curve().at(id).vertex().at(0);
      p2 = m->curves().curve().at(id).vertex().at(1);

      *en = mesh->peek_edge_node(p1, p2);
      if (*en == NULL)
        error("Curve #%d: edge %d-%d does not exist.", id, p1, p2);

      // degree of curved edge
      nurbs->degree = 2;

      // get the number of control points
      std::vector<double> vCP;

      int inner = 1, outer;

      nurbs->np = inner + 2;

      // edge endpoints are also control points, with weight 1.0
      nurbs->pt = new double3[nurbs->np];
      nurbs->pt[0][0] = mesh->nodes[p1].x;
      nurbs->pt[0][1] = mesh->nodes[p1].y;
      nurbs->pt[0][2] = 1.0;
      nurbs->pt[inner+1][0] = mesh->nodes[p2].x;
      nurbs->pt[inner+1][1] = mesh->nodes[p2].y;
      nurbs->pt[inner+1][2] = 1.0;

      // read the arc angle
      nurbs->angle = m->curves().curve().at(id).angle();
      double a = (180.0 - nurbs->angle) / 180.0 * M_PI;

      // generate one control point
      double x = 1.0 / tan(a * 0.5);
      nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
      nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
      nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

      // get the number of knot vector points
      std::vector<double> vKnots;

      inner = 0;

      nurbs->nk = nurbs->degree + nurbs->np + 1;
      outer = nurbs->nk - inner;
      if ((outer & 1) == 1)
        error("Curve #%d: incorrect number of knot points.", id);

      // knot vector is completed by 0.0 on the left and by 1.0 on the right
      nurbs->kv = new double[nurbs->nk];

      for (int i = 0; i < outer/2; i++)
        nurbs->kv[i] = 0.0;

      if (inner) {
        for (int i = outer/2; i < inner + outer/2; i++) {
          nurbs->kv[i] = vKnots[i - (outer/2)];
        }
      }

      for (int i = outer/2 + inner; i < nurbs->nk; i++)
        nurbs->kv[i] = 1.0;

      nurbs->ref = 0;

      return nurbs;
    }
  }
}