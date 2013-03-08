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

        /// Main method - processes the solution and stores the data obtained by the process.
        /// \param[in] xsln the first solution (in the x-direction)
        /// \param[in] ysln the second solution (in the y-direction)
        /// \param[in] xitem what item (function value, derivative wrt. x, ..) to use in the first solution.
        /// \param[in] yitem what item (function value, derivative wrt. x, ..) to use in the second solution.
        /// \param[in] eps - tolerance parameter controlling how fine the resulting linearized approximation of the solution is.
        void process_solution(MeshFunctionSharedPtr<double> xsln, MeshFunctionSharedPtr<double> ysln, int xitem = H2D_FN_VAL_0, int yitem = H2D_FN_VAL_0, double eps = HERMES_EPS_NORMAL);

        /// Set the displacement, i.e. set two functions that will deform the domain for visualization, in the x-direction, and the y-direction.
        void set_displacement(MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult = 1.0);

        /// Get the number of vertices of this instance.
        int get_num_vertices();

        /// Get the vertices of this instance.
        double4* get_vertices();

        /// Get the dashes (the little arrows) of this instance.
        int2* get_dashes();

        /// Get the number of dashes (the little arrows) of this instance.
        int get_num_dashes();

        void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        /// Sets the threshold for how fine the output for curved elements.
        /// \param[in] curvature_epsilon The 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        /// The smaller, the finer.
        /// Default value = 1e-3.
        void set_curvature_epsilon(double curvature_epsilon);

        /// Gets the 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        double get_curvature_epsilon();

        /// Frees the instance.
        void free();

      protected:
        /// The 'curvature' epsilon.
        double curvature_epsilon;

        /// Information if user-supplied displacement functions have been provided.
        bool user_xdisp, user_ydisp;
        
        /// Displacement functions, default to ZeroFunctions, may be supplied by set_displacement();
        MeshFunctionSharedPtr<double> xdisp, ydisp;
        double dmult;

        int xitem, component_x, value_type_x;
        int yitem, component_y, value_type_y;

        double4* verts;  ///< vertices: (x, y, xvalue, yvalue) quadruples
        int2* dashes;

        int dashes_count; ///< Real numbers of vertices, triangles and edges, dashes
        int dashes_size; ///< Size of arrays of vertices, triangles and edges, dashes

        int get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue);
        void process_dash(int iv1, int iv2);

        int add_vertex();

        void add_dash(int iv1, int iv2);

        void process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
          double* xval, double* yval, double* phx, double* phy, int* indices, bool curved);

        void process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
          double* xval, double* yval, double* phx, double* phy, int* indices, bool curved);

        void find_min_max();

        /// Internal.
        void push_transforms(MeshFunction<double>** fns, int transform);

        /// Internal.
        void pop_transforms(MeshFunction<double>** fns);
      };
    }
  }
}
#endif
