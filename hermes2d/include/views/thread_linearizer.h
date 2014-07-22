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
/*! \file thread_linearizer.h
\brief File containing ThreadLinearizerMultidimensional class.
*/

#ifndef __H2D_THREAD_LINEARIZER_H
#define __H2D_THREAD_LINEARIZER_H

#include "linearizer.h"
#include "linearizer_utils.h"
#include <pthread.h>
#include "../function/mesh_function.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      template<typename LinearizerDataDimensions>
      class HERMES_API LinearizerMultidimensional;

      /// ThreadLinearizerMultidimensional is a utility class for linearizing a mesh function on a single thread
      /// The main refinement (splitting) decision is contained in split_decision().
      /// Important - the instance of LinearizerCriterion (from LinearizerMultidimensional instance) - see the method split_decision() for the adaptive criterion, process_[triangle|quad] for the fixed one.
      template<typename LinearizerDataDimensions>
      class HERMES_API ThreadLinearizerMultidimensional
      {
      public:
        /// Free the instance.
        void free();

      private:
        /// Constructor
        /// \param[in] The linearizer this instance is being created for.
        ThreadLinearizerMultidimensional(LinearizerMultidimensional<LinearizerDataDimensions>* linearizer);
        ~ThreadLinearizerMultidimensional();

        /// Get the data from the "parent" LinearizerMultidimensional class.
        void init_linearizer_data(LinearizerMultidimensional<LinearizerDataDimensions>* linearizer);

        /// Initialize arrays, clone functions etc for this run of processing.
        void init_processing(MeshFunctionSharedPtr<double>* sln, LinearizerMultidimensional<LinearizerDataDimensions>* linearizer);
        /// Deinitialize the temporary data for this run of processing.
        void deinit_processing();

        /// Completely process the state current_state
        void process_state(Traverse::State* current_state);

        /// Return the hash value of the couple of vertices with indices p1, p2.
        int hash(int p1, int p2);
        /// Return the index of the vertex between vertices with indices p1, p2.
        int peek_vertex(int p1, int p2);

        /// Process the edge between vertices with indices iv1, iv2.
        void process_edge(int iv1, int iv2, int marker);
        /// Add the edge to the resulting data.
        void add_edge(int iv1, int iv2, int marker);
        /// Add the triangle to the resulting data.
        void add_triangle(int iv0, int iv1, int iv2, int marker);

        /// Add a blank new vertex.
        int add_vertex();
        /// Return the [existing|new] vertex between p1 and p2, uses add_vertex() for new vertex creation.
        int get_vertex(int p1, int p2, double x, double y, double* value);

        /// Process a triangle with vertices iv0, iv1, iv2.
        /// Recursive.
        /// \param[in] level The current level of refinement
        void process_triangle(int iv0, int iv1, int iv2, int level);
        /// Process a quad with vertices iv0, iv1, iv2, iv3.
        /// Recursive.
        /// \param[in] level The current level of refinement
        void process_quad(int iv0, int iv1, int iv2, int iv3, int level);

        /// Main method for deciding whether or not to split the currently evaluated element.
        /// \param[in/out] split Integer filled with the resulting split of the element:
        /// triangle:: split > 0 ? split to four triangles : no split.
        /// quad:: split == 0 -> no split.
        /// quad:: split == 1 -> horizontal split.
        /// quad:: split == 2 -> vertical split.
        /// quad:: split == 3 -> split to four quads.
        void split_decision(int& split, int iv0, int iv1, int iv2, int iv3, ElementMode2D mode, const double** val, double* phx, double* phy, unsigned short* indices) const;

        /// Utility - check of the orientation of tthe two triangles outputted for a quad.
        bool quad_flip(int iv0, int iv1, int iv2, int iv3) const;

        /// Internal.
        double get_max_value(Traverse::State* current_state);
        double max_value_approx;

        /// Push transforms to all necessary functions (including displacement).
        void push_transforms(int transform);

        /// Pop transforms from all necessary functions (including displacement).
        void pop_transforms();

        /// Reallocation of the data fields.
        /// Done in the initialization of a processing run.
        void reallocate(MeshSharedPtr mesh);

        /// Thread-owned clones.
        MeshFunction<double>* fns[LinearizerDataDimensions::dimension + 2];

        /// Assigned criterion.
        /// See the class LinearizerCriterion.
        /// Used by split_decision().
        LinearizerCriterion criterion;

        /// Data - OpenGL part.
        typename LinearizerDataDimensions::triangle_t* triangles;
        typename LinearizerDataDimensions::edge_t* edges;
        /// - edge_markers: edge markers, ordering equal to edges
        int* edge_markers;

        /// Data - FileExport part.
        /// Vertices: (x, y, value) triplets
        /// - triangles: vertex index triplets
        triangle_indices_t* triangle_indices;

        /// Data - Common part.
        typename LinearizerDataDimensions::vertex_t* vertices;
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
        double midval[LinearizerDataDimensions::dimension + 2][5];
        /// Temporary storage - per state processing.
        const double* val[LinearizerDataDimensions::dimension];
        /// Representing element of the currently processed state.
        Element* rep_element;
        /// The current state is curved.
        bool curved;

        /// From LinearizerMultidimensional - for convenience & speed.
        LinearizerOutputType linearizerOutputType;

        /// The information do we want to get out of the solution.
        int item[LinearizerDataDimensions::dimension], component[LinearizerDataDimensions::dimension], value_type[LinearizerDataDimensions::dimension];
        /// User displacement is present.
        bool user_xdisp, user_ydisp;
        /// Multiplication factor of the displacement function.
        double dmult;
        /// Standard and curvature epsilon.
        double curvature_epsilon;

        friend class LinearizerMultidimensional<LinearizerDataDimensions>;
      };
    }
  }
}
#endif
