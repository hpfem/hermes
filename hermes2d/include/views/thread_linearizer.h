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

#ifndef __H2D_THREAD_LINEARIZER_H
#define __H2D_THREAD_LINEARIZER_H

#include "linearizer.h"
#include "../function/mesh_function.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      class HERMES_API LinearizerNew;

      /// ThreadLinearizerNew is a utility class for linearizing a mesh function on a single thread
      class HERMES_API ThreadLinearizerNew
      {
      public:
        /// Free the instance.
        void free();
        
      private:
        ThreadLinearizerNew(LinearizerNew* linearizer);
        ~ThreadLinearizerNew();

        void init_linearizer_data(LinearizerNew* linearizer);

        void init_processing(MeshFunctionSharedPtr<double> sln, LinearizerNew* linearizer);
        void deinit_processing();
        
        void process_state(Traverse::State* current_state);

        void set_min_value(double min);

        void set_max_value(double max);

        int hash(int p1, int p2);
        int peek_vertex(int p1, int p2);

        void process_edge(int iv1, int iv2, int marker);
        void add_edge(int iv1, int iv2, int marker);
        void add_triangle(int iv0, int iv1, int iv2, int marker);

        int add_vertex();
        int get_vertex(int p1, int p2, double x, double y, double value);

        void process_triangle(int iv0, int iv1, int iv2, int level, double* val, double* phx, double* phy, int* indices);

        void process_quad(int iv0, int iv1, int iv2, int iv3, int level, double* val, double* phx, double* phy, int* indices);

        void split_decision(int& split, int iv0, int iv1, int iv2, int iv3, ElementMode2D mode, double* val, double* phx, double* phy, int* indices) const;

        bool quad_flip(int iv0, int iv1, int iv2, int iv3) const;

        void find_min_max();
        
        /// Internal.
        void push_transforms(int transform);

        /// Internal.
        void pop_transforms();

        void reallocate(MeshSharedPtr mesh);

        /// Thread-owned clones.
        MeshFunction<double>* fns[3];

        // OpenGL part.
        triangle_t* triangles;
        edge_t* edges;
        /// - edge_markers: edge markers, ordering equal to edges
        int* edge_markers;

        // FileExport part.
        /// Vertices: (x, y, value) triplets
        /// - triangles: vertex index triplets
        triangle_indices_t* triangle_indices;

        // Common part.
        vertex_t* vertices;
        /// - triangle_markers: triangle markers, ordering equal to triangles, triangle_indices
        int* triangle_markers;
        /// - hash table
        int* hash_table;
        /// info[0] = p1, info[1] = p2, info[2] = next vertex in hash
        internal_vertex_info_t* info;

        /// Real counts of vertices, triangles and edges
        int vertex_count, triangle_count, edges_count;
        /// Size of arrays of vertices, triangles and edges
        int vertex_size, triangle_size, edges_size;

        
        /// Temporary storage - per state processing.
        double midval[3][5];
        Element* rep_element;
        bool curved;

        
        /// From LinearizerNew - for convenience & speed.
        LinearizerOutputType linearizerOutputType;

        /// The information do we want to get out of the solution.
        int item, component, value_type;
        bool user_xdisp, user_ydisp;
        double dmult;
        /// Standard and curvature epsilon.
        double epsilon, curvature_epsilon;


        /// Keep?
        bool user_specified_max, user_specified_min;
        double user_specified_max_value, user_specified_min_value;

        friend class LinearizerNew;
      };
    }
  }
}
#endif
