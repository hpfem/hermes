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

#ifndef __H2D_VECTORIZER_H
#define __H2D_VECTORIZER_H

#include "linearizer.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// \brief "Vectorizer" is a Linearizer for vector solutions.
      /// The only difference is that linearized vertices are vector-valued. Also, regularization of the
      /// resulting mesh is not attempted. The class can handle different meshes in
      /// both X and Y components.
      ///
      class HERMES_API Vectorizer : public LinearizerBase
      {
      public:

        Vectorizer();
        ~Vectorizer();

        void process_solution(MeshFunction<double>* xsln, MeshFunction<double>* ysln, int xitem, int yitem, double eps);

        int get_num_vertices();
        double4* get_vertices();

        int2* get_dashes();
        int get_num_dashes();

        void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        /// Frees the instance.
        void free();
      protected:

        MeshFunction<double>*xsln, *ysln;
        int xitem, component_x, value_type_x;
        int yitem, component_y, value_type_y;

        double4* verts;  ///< vertices: (x, y, xvalue, yvalue) quadruples
        int2* dashes;

        int dashes_count; ///< Real numbers of vertices, triangles and edges, dashes
        int dashes_size; ///< Size of arrays of vertices, triangles and edges, dashes

        int get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue);
        int create_vertex(double x, double y, double xvalue, double yvalue);
        void process_dash(int iv1, int iv2);

        int add_vertex();

        double cmax;

        void add_dash(int iv1, int iv2);

        void push_transform(int son);

        void pop_transform();

        void process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
          double* xval, double* yval, double* phx, double* phy, int* indices, bool curved);

        void process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
          double* xval, double* yval, double* phx, double* phy, int* indices, bool curved);

        void find_min_max();
      };
    }
  }
}
#endif