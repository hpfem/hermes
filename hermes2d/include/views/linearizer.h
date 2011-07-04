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

#ifndef __H2D_LINEARIZER_H
#define __H2D_LINEARIZER_H

#include "../hermes2d_common_defs.h"
#include "../function/solution.h"

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
        Quad2DLin()
        {
          mode = HERMES_MODE_TRIANGLE;
          max_order[0]  = max_order[1]  = 1;
          num_tables[0] = num_tables[1] = 2;
          tables = lin_tables;
          np = lin_np;
        };
      };

      /// Linearizer<Scalar> is a utility class which converts a higher-order FEM solution defined on
      /// a curvilinear, irregular mesh to a linear FEM solution defined on a straight-edged,
      /// regular mesh. This is done by adaptive refinement of the higher-order mesh and its
      /// subsequent regularization. The linearized mesh can then be easily displayed or
      /// exported to standard formats. The class correctly handles discontinuities in the
      /// solution (e.g., gradients or in Hcurl) by inserting double vertices where necessary.
      /// Linearizer<Scalar> also serves as a container for the resulting linearized mesh.
      ///
      template<typename Scalar>
      class HERMES_API Linearizer
      {
      public:

        Linearizer();
        ~Linearizer();

        void process_solution(MeshFunction<Scalar>* sln, int item = H2D_FN_VAL_0,
          double eps = HERMES_EPS_NORMAL, double max_abs = -1.0,
          MeshFunction<Scalar>* xdisp = NULL, MeshFunction<Scalar>* ydisp = NULL,
          double dmult = 1.0);

        void lock_data() const { pthread_mutex_lock(&data_mutex); }
        void unlock_data() const { pthread_mutex_unlock(&data_mutex); }

        double3* get_vertices() const { return verts; }
        int get_num_vertices() const { return nv; }

        int3* get_triangles() const { return tris; }
        int get_num_triangles() const { return nt; }

        int3* get_edges() const { return edges; }
        int get_num_edges() const { return ne; }

        double get_min_value() const { return min_val; }
        double get_max_value() const { return max_val; }
        virtual void calc_vertices_aabb(double* min_x, double* max_x, 
          double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        /// Saves data in a binary format.
        virtual void save_data(const char* filename);

        /// Loads data in a binary format.
        virtual void load_data(const char* filename);

        /// Saves a MeshFunction (Solution, Filter) in VTK format.
        virtual void save_solution_vtk(MeshFunction<Scalar>* meshfn, const char* file_name, const char* quantity_name,
          bool mode_3D = true, int item = H2D_FN_VAL_0, 
          double eps = HERMES_EPS_NORMAL, double max_abs = -1.0,
          MeshFunction<Scalar>* xdisp = NULL, MeshFunction<Scalar>* ydisp = NULL,
          double dmult = 1.0);

        /// This function is used by save_solution_vtk().
        virtual void save_data_vtk(const char* file_name, const char* quantity_name, bool mode_3D);

        void free();

      protected:

        Quad2DLin quad_lin;
        MeshFunction<Scalar>* sln;
        int item, ia, ib;

        double eps, max, cmax;
        bool auto_max;

        MeshFunction<Scalar> *xdisp, *ydisp;
        double dmult;

        double3* verts;  ///< vertices: (x, y, value) triplets
        int4* info;      ///< info[0]=p1, info[1]=p2, info[2]=next vertex in hash
        int3* tris;      ///< triangles: vertex index triplets
        int3* edges;     ///< edges: pairs of vertex indices
        int* hash_table; ///< hash table

        int nv, nt, ne; ///< numbers of vertices, triangles and edges
        int cv, ct, ce; ///< capacities of vertex, triangle and edge arrays
        int del_slot;   ///< free slot index after a triangle which was deleted
        int mask;       ///< hash table mask = size-1

        bool curved, disp;
        double min_val, max_val;

        int get_vertex(int p1, int p2, double x, double y, double value);
        int get_top_vertex(int id, double value);
        int peek_vertex(int p1, int p2);

        int hash(int p1, int p2) { return (984120265*p1 + 125965121*p2) & mask; }

        int add_vertex()
        {
          if (nv >= cv)
          {
            cv *= 2;
            verts = (double3*) realloc(verts, sizeof(double3) * cv);
            info = (int4*) realloc(info, sizeof(int4) * cv);
            verbose("Linearizer::add_vertex(): realloc to %d", cv);
          }
          return nv++;
        }

        void add_triangle(int iv0, int iv1, int iv2);

        void del_triangle(int index)
        {
          del_slot = index;
        }

        void add_edge(int iv1, int iv2, int marker)
        {
          if (ne >= ce) edges = (int3*) realloc(edges, sizeof(int3) * (ce = ce * 3 / 2));
          edges[ne][0] = iv1;
          edges[ne][1] = iv2;
          edges[ne++][2] = marker;
        }

        void get_gv_a_b(int mask, int& a, int& b)
        {
          a = b = 0;
          if (mask >= 0x40) { a = 1; mask >>= 6; }
          while (!(mask & 1)) { mask >>= 1; b++; }
        }

        void process_triangle(int iv0, int iv1, int iv2, int level,
          Scalar* val, double* phx, double* phy, int* indices);

        void process_quad(int iv0, int iv1, int iv2, int iv3, int level,
          Scalar* val, double* phx, double* phy, int* indices);

        void process_edge(int iv1, int iv2, int marker);
        void regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2);
        void find_min_max();
        void print_hash_stats();

        mutable pthread_mutex_t data_mutex;

        static void calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y); ///< Calculates AABB from an array of X-axis and Y-axis coordinates. The distance between values in the array is stride bytes.
      };

      /// maximum subdivision level (2^N)
      const int LIN_MAX_LEVEL = 6;

#define lin_init_array(array, type, c, e) \
  if (c < e) { \
  if (array != NULL) ::free(array); \
  array = (type*) malloc(sizeof(type) * (c = e)); }

#define lin_free_array(array, n, c) \
  if (array != NULL) { \
  ::free(array); array = NULL; \
  n = c = 0; }

    }
  }
}
#endif
