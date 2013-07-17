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

#ifndef _MESH_READER_H2D_BSON_H_
#define _MESH_READER_H2D_BSON_H_

#include "mesh_reader.h"

#ifdef WITH_BSON
#include "bson.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /// Mesh reader from BSON format
    ///
    /// @ingroup mesh_readers
    /// Typical usage:
    /// MeshSharedPtr mesh;
    /// Hermes::Hermes2D::MeshReaderH2DBSON mloader;
    /// try
    /// {
    ///&nbsp;mloader.load("mesh.bson", &mesh);
    /// }
    /// catch(Exceptions::MeshLoadFailureException& e)
    /// {
    ///&nbsp;e.print_msg();
    ///&nbsp;return -1;
    /// }
    /// 
    class HERMES_API MeshReaderH2DBSON : public MeshReader
    {
    public:
      MeshReaderH2DBSON();
      virtual ~MeshReaderH2DBSON();

      /// This method loads a single mesh from a file.
      virtual void load(const char *filename, MeshSharedPtr mesh);

      /// This method saves a single mesh to a file.
      void save(const char *filename, MeshSharedPtr mesh);

      /// This method loads multiple meshes according to subdomains described in the meshfile.
      /// \param[in] meshes Meshes to be loaded, the number must correspond to the subdomains described in the file.
      ///&nbsp;         also the order is determined by the order in the file.
      void load(const char *filename, Hermes::vector<MeshSharedPtr > meshes);

      /// This method saves multiple meshes according to subdomains in the vector meshes.
      void save(const char *filename, Hermes::vector<MeshSharedPtr > meshes);

    private:
      /// Loads one circular arc.
      /// \param[in] skip_check Skip check that the edge exists, in case of subdomains.
      Nurbs* load_arc(MeshSharedPtr mesh, int id, Node** en, int p1, int p2, double angle, bool skip_check = false);

      struct vertex_BSON
      {
        vertex_BSON(){}
        vertex_BSON(double x, double y, int i) : x(x), y(y), i(i) {}
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "vertex");
          bson_append_double(&bw, "x", x);
          bson_append_double(&bw, "y", y);
          bson_append_int(&bw, "i", i);
          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson_iterator sub_it;
          bson_find(&sub_it, &sub, "x");
          this->x = bson_iterator_double(&sub_it);
          bson_find(&sub_it, &sub, "y");
          this->y = bson_iterator_double(&sub_it);
          bson_find(&sub_it, &sub, "i");
          this->i = bson_iterator_int(&sub_it);
        }
        double x;
        double y;
        int i;
      };

      struct edge_BSON
      {
        edge_BSON(){}
        edge_BSON(int v1, int v2, std::string marker, int i) : v1(v1), v2(v2), marker(marker), i(i) {}
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "edge");
          bson_append_int(&bw, "v1", v1);
          bson_append_int(&bw, "v2", v2);
          bson_append_string(&bw, "marker", marker.c_str());
          bson_append_int(&bw, "i", i);
          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson_iterator sub_it;
          bson_find(&sub_it, &sub, "v1");
          this->v1 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "v2");
          this->v2 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "marker");
          this->marker = bson_iterator_string(&sub_it);
          bson_find(&sub_it, &sub, "i");
          this->i = bson_iterator_int(&sub_it);
        }
        int v1;
        int v2;
        std::string marker;
        int i;
      };

      struct element_BSON
      {
        element_BSON(){}
        element_BSON(int v1, int v2, int v3, int v4, std::string marker, int i) : v1(v1), v2(v2), v3(v3), v4(v4), marker(marker), i(i) {}
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "element");
          bson_append_int(&bw, "v1", v1);
          bson_append_int(&bw, "v2", v2);
          bson_append_int(&bw, "v3", v3);
          bson_append_int(&bw, "v4", v4);
          bson_append_string(&bw, "marker", marker.c_str());
          bson_append_int(&bw, "i", i);
          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson_iterator sub_it;
          bson_find(&sub_it, &sub, "v1");
          this->v1 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "v2");
          this->v2 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "v3");
          this->v3 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "v4");
          this->v4 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "marker");
          this->marker = bson_iterator_string(&sub_it);
          bson_find(&sub_it, &sub, "i");
          this->i = bson_iterator_int(&sub_it);
        }
        int v1;
        int v2;
        int v3;
        int v4;
        std::string marker;
        int i;
      };

      struct arc_BSON
      {
        arc_BSON(){}
        arc_BSON(int p1, int p2, double angle) : p1(p1), p2(p2), angle(angle) {}
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "arc");
          bson_append_int(&bw, "p1", p1);
          bson_append_int(&bw, "p2", p2);
          bson_append_double(&bw, "angle", angle);
          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson_iterator sub_it;
          bson_find(&sub_it, &sub, "p1");
          this->p1 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "p2");
          this->p2 = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "angle");
          this->angle = bson_iterator_double(&sub_it);
        }
        int p1;
        int p2;
        double angle;
      };

      struct refinement_BSON
      {
        refinement_BSON() {}
        refinement_BSON(int id, int type) : id(id), type(type) {}
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "refinement");
          bson_append_int(&bw, "id", id);
          bson_append_int(&bw, "type", type);
          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson_iterator sub_it;
          bson_find(&sub_it, &sub, "id");
          this->id = bson_iterator_int(&sub_it);
          bson_find(&sub_it, &sub, "type");
          this->type = bson_iterator_int(&sub_it);
        }
        int id;
        int type;
      };

      struct subdomain_BSON
      {
        void save_to_BSON(bson& bw)
        {
          bson_append_start_object(&bw, "subdomain");
          // vertices
          bson_append_start_array(&bw, "vertices");
          for (int i = 0; i < vertices.size(); i++)
            bson_append_int(&bw, "i", vertices[i]);
          bson_append_finish_array(&bw);

          // elements
          bson_append_start_array(&bw, "elements");
          for (int i = 0; i < elements.size(); i++)
            bson_append_int(&bw, "i", elements[i]);
          bson_append_finish_array(&bw);

          // boundary edges
          bson_append_start_array(&bw, "boundary_edges");
          for (int i = 0; i < boundary_edges.size(); i++)
            bson_append_int(&bw, "i", boundary_edges[i]);
          bson_append_finish_array(&bw);

          // inner edges
          bson_append_start_array(&bw, "inner_edges");
          for (int i = 0; i < inner_edges.size(); i++)
            bson_append_int(&bw, "i", inner_edges[i]);
          bson_append_finish_array(&bw);

          // refinements
          bson_append_start_array(&bw, "refinements");
          for (int i = 0; i < refinements.size(); i++)
            refinements[i].save_to_BSON(bw);
          bson_append_finish_array(&bw);

          bson_append_finish_object(&bw);
        }
        void load_from_BSON(bson& sub)
        {
          bson sub_sub;
          bson_iterator sub_it, sub_sub_it;
          // vertices
          bson_find(&sub_it, &sub, "vertices");
          bson_iterator_subobject_init(&sub_it, &sub_sub, 0);
          bson_iterator_init(&sub_sub_it, &sub_sub);
          while (bson_iterator_next(&sub_sub_it))
            this->vertices.push_back(bson_iterator_int(&sub_sub_it));

          // elements
          bson_find(&sub_it, &sub, "elements");
          bson_iterator_subobject_init(&sub_it, &sub_sub, 0);
          bson_iterator_init(&sub_sub_it, &sub_sub);
          while (bson_iterator_next(&sub_sub_it))
            this->elements.push_back(bson_iterator_int(&sub_sub_it));

          // boundary_edges
          bson_find(&sub_it, &sub, "boundary_edges");
          bson_iterator_subobject_init(&sub_it, &sub_sub, 0);
          bson_iterator_init(&sub_sub_it, &sub_sub);
          while (bson_iterator_next(&sub_sub_it))
            this->boundary_edges.push_back(bson_iterator_int(&sub_sub_it));

          // inner_edges
          bson_find(&sub_it, &sub, "inner_edges");
          bson_iterator_subobject_init(&sub_it, &sub_sub, 0);
          bson_iterator_init(&sub_sub_it, &sub_sub);
          while (bson_iterator_next(&sub_sub_it))
            this->inner_edges.push_back(bson_iterator_int(&sub_sub_it));

          // refinements
          bson_find(&sub_it, &sub, "refinements");
          bson_iterator_subobject_init(&sub_it, &sub_sub, 0);
          bson_iterator_init(&sub_sub_it, &sub_sub);
          bson b_sub_sub;
          while (bson_iterator_next(&sub_sub_it))
          {
            refinement_BSON refinement;
            bson_iterator_subobject_init(&sub_sub_it, &b_sub_sub, 0);
            refinement.load_from_BSON(b_sub_sub);
            this->refinements.push_back(refinement);
          }
        }
        Hermes::vector<int> vertices;
        Hermes::vector<int> elements;
        Hermes::vector<int> boundary_edges;
        Hermes::vector<int> inner_edges;
        Hermes::vector<refinement_BSON> refinements;
      };

      static bool elementCompare (element_BSON el_i, element_BSON el_j) { return ( el_i.i < el_j.i ); }

      void load_domain(bson& br, MeshSharedPtr mesh, std::map<int, int>& vertex_is, std::map<int, int>& element_is, std::map<int, int>& edge_is,
        Hermes::vector<element_BSON>& elements, Hermes::vector<edge_BSON>& edges, Hermes::vector<vertex_BSON>& vertices, Hermes::vector<arc_BSON>& arcs, Hermes::vector<subdomain_BSON>& subdomains);
    };
  }
}
#endif
#endif

