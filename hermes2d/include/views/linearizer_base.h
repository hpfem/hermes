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

#ifndef __H2D_LINEARIZER_SUPPORT_H
#define __H2D_LINEARIZER_SUPPORT_H

#include "global.h"
#include "../quadrature/quad_all.h"
#include "../mesh/mesh.h"
#include "mixins2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      const double HERMES_EPS_LOW      = 0.007;
      const double HERMES_EPS_NORMAL   = 0.0004;
      const double HERMES_EPS_HIGH     = 0.0001;
      const double HERMES_EPS_VERYHIGH = 0.000002;

      //// linearization "quadrature" ////////////////////////////////////////////////////////////////////

      /// The tables with index zero are for obtaining solution values at the element
      /// vertices. Index one tables serve for the retrieval of interior values. Index one tables
      /// are used for adaptive approximation of the solution by transforming their points to sub-elements.
      /// Actually, the tables contain two levels of refinement -- this is an optimization to reduce
      /// the number of calls to sln->get_values().
      extern double3 lin_pts_0_tri[];

      extern double3 lin_pts_0_quad[];

      extern double3 lin_pts_1_tri[12];

      extern double3 lin_pts_1_quad[21];

      extern int quad_indices[9][5];

      extern int tri_indices[5][3];

      extern int lin_np_tri[2];
      extern int lin_np_quad[2];
      extern int* lin_np[2];

      extern double3*  lin_tables_tri[2];
      extern double3*  lin_tables_quad[2];
      extern double3** lin_tables[2];

      class Quad2DLin : public Quad2D
      {
      public:
        Quad2DLin();
      };

      extern HERMES_API Quad2DLin g_quad_lin;

      /// Base class for Linearizer, Orderizer, Vectorizer.

      class HERMES_API LinearizerBase : 
        public Hermes::Mixins::TimeMeasurable,
        public Hermes::Mixins::Loggable,
        public Hermes::Hermes2D::Mixins::Parallel
      {
      public:
        void set_max_absolute_value(double max_abs);

        double get_min_value() const;
        double get_max_value() const;

        void lock_data() const;
        void unlock_data() const;

        int3* get_triangles();
        int* get_triangle_markers();
        int get_num_triangles();
        int2* get_edges();
        int* get_edge_markers();
        int get_num_edges();

        /// The instance is empty. Either process_solution has not been called so far, or
        /// the method free() has been called.
        virtual bool is_empty();

        /// Initializing of a linearizing process.
        virtual void init_linearizer_base(MeshFunctionSharedPtr<double> sln);
        virtual void deinit_linearizer_base();

        /// Frees the instance.
        void free();
        
        /// Experimental upper-limiting of the maximum refinement level.
        int get_max_level(Element* e, int polynomial_order, MeshSharedPtr mesh);
        static double large_elements_fraction_of_mesh_size_threshold;
        int* level_map;
        
      protected:
        LinearizerBase(bool auto_max = true);
        ~LinearizerBase();

        void process_edge(int iv1, int iv2, int marker);

        bool empty;

        double max;

        void regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2, int marker);

        bool auto_max;

        int3* tris;      ///< triangles: vertex index triplets
        int* tri_markers;///< triangle_markers: triangle markers, ordering equal to tris
        int2* edges;     ///< edges: pairs of vertex indices
        int* edge_markers;     ///< edge_markers: edge markers, ordering equal to edges
        int* hash_table; ///< hash table
        int4 * info; ///< info[0] = p1, info[1] = p2, info[2] = next vertex in hash

        int vertex_count, triangle_count, edges_count; ///< Real numbers of vertices, triangles and edges
        int vertex_size, triangle_size, edges_size; ///< Size of arrays of vertices, triangles and edges

        double eps;

        double min_val, max_val;

        int del_slot;   ///< free slot index after a triangle which was deleted

        int peek_vertex(int p1, int p2);

        void add_edge(int iv1, int iv2, int marker);
        void add_triangle(int iv0, int iv1, int iv2, int marker);

        int hash(int p1, int p2);

        mutable pthread_mutex_t data_mutex;

        /// Calculates AABB from an array of X-axis and Y-axis coordinates. The distance between values in the array is stride bytes.
        static void calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y);
        friend class MeshView;
        friend class OrderView;
        friend class ScalarView;
        friend class VectorView;
        friend class StreamView;
      };

      const int LIN_MAX_LEVEL = 5;
    }
  }
}

#endif
