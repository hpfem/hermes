// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ORDERIZER_H
#define __H2D_ORDERIZER_H

#include "../space/space.h"
#include "global.h"
#ifndef NOGLUT
  #include <pthread.h>
#endif
#include "../quadrature/quad_all.h"
#include "../mesh/traverse.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// Like the Linearizer, but generates a triangular mesh showing polynomial
      /// orders in a space, hence the funky name.
      ///
      class HERMES_API Orderizer
      {
      public:

        Orderizer();
        ~Orderizer();

        /// Saves the polynomial orders.
        template<typename Scalar>
        void save_orders_vtk(SpaceSharedPtr<Scalar> space, const char* file_name);

        /// Saves the mesh with markers.
        template<typename Scalar>
        void save_markers_vtk(SpaceSharedPtr<Scalar> space, const char* file_name);

        /// Saves the mesh - edges.
        template<typename Scalar>
        void save_mesh_vtk(SpaceSharedPtr<Scalar> space, const char* file_name);

        /// Returns axis aligned bounding box (AABB) of vertices. Assumes lock.
        void calc_vertices_aabb(double* min_x, double* max_x,
          double* min_y, double* max_y) const;

        /// Internal.
        template<typename Scalar>
        void process_space(SpaceSharedPtr<Scalar> space, bool show_edge_orders = false);

        int get_labels(int*& lvert, char**& ltext, double2*& lbox) const;

        int get_num_vertices();
        double3* get_vertices();

        void lock_data() const;
        void unlock_data() const;

        int3* get_triangles();
        int* get_triangle_markers();
        int get_num_triangles();

        int2* get_edges();
        int* get_edge_markers();
        int get_num_edges();

        void free();
      protected:
#ifndef NOGLUT
        mutable pthread_mutex_t data_mutex;
#endif
        int2* edges;     ///< edges: pairs of vertex indices
        int* edge_markers;     ///< edge_markers: edge markers, ordering equal to edges
        void add_edge(int iv1, int iv2, int marker);

        /// Reallocation at the beginning of process_*.
        /// Specific for Linearizer
        void reallocate(MeshSharedPtr mesh);

        char  buffer[1000];
        char* labels[11][11];

        double3* verts;  ///< vertices: (x, y, value) triplets
        int  label_size, label_count;
        int* lvert;
        char** ltext;
        double2* lbox;

        int3* tris;      ///< triangles: vertex index triplets
        int* tri_markers;///< triangle_markers: triangle markers, ordering equal to tris

        int vertex_count, triangle_count, edges_count; ///< Real numbers of vertices, triangles and edges
        int vertex_size, triangle_size, edges_size; ///< Size of arrays of vertices, triangles and edges

        void add_triangle(int iv0, int iv1, int iv2, int marker);

        static void calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y);
        
        int add_vertex();

        void make_vert(int & index, double x, double y, double val);
      };
    }
  }
}
#endif