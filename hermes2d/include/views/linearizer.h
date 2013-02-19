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

#include "../global.h"
#include "../function/solution.h"
#include "linearizer_base.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// Linearizer is a utility class which converts a higher-order FEM solution defined on
      /// a curvilinear, irregular mesh to a linear FEM solution defined on a straight-edged,
      /// regular mesh. This is done by adaptive refinement of the higher-order mesh and its
      /// subsequent regularization. The linearized mesh can then be easily displayed or
      /// exported to standard formats. The class correctly handles discontinuities in the
      /// solution (e.g., gradients or in Hcurl) by inserting double vertices where necessary.
      /// Linearizer also serves as a container for the resulting linearized mesh.
      ///
      class HERMES_API Linearizer : public LinearizerBase
      {
      public:

        Linearizer(bool auto_max = true);
        ~Linearizer();

        /// Main method - processes the solution and stores the data obtained by the process.
        /// \param[in] sln the solution
        /// \param[in] item what item (function value, derivative wrt. x, ..) to use in the solution.
        /// \param[in] eps - tolerance parameter controlling how fine the resulting linearized approximation of the solution is.
        void process_solution(MeshFunction<double>* sln, int item = H2D_FN_VAL_0, double eps = HERMES_EPS_NORMAL);

        /// Save a MeshFunction (Solution, Filter) in VTK format.
        void save_solution_vtk(MeshFunction<double>* sln, const char* filename, const char* quantity_name,
          bool mode_3D = true, int item = H2D_FN_VAL_0,
          double eps = HERMES_EPS_NORMAL);

        /// Set the displacement, i.e. set two functions that will deform the domain for visualization, in the x-direction, and the y-direction.
        void set_displacement(MeshFunction<double>* xdisp, MeshFunction<double>* ydisp, double dmult = 1.0);

        void calc_vertices_aabb(double* min_x, double* max_x,
          double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        /// Get the number of vertices of this instance.
        int get_num_vertices();

        /// Get the vertices of this instance.
        double3* get_vertices();

        /// Get the contours (the isolines) of this instance.
        int3* get_contour_triangles();

        /// Get the number of contours (the isolines) of this instance.
        int get_num_contour_triangles();

        /// Set the threshold for how fine the output for curved elements.
        /// \param[in] curvature_epsilon The 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        /// The smaller, the finer.
        /// Default value = 1e-3.
        void set_curvature_epsilon(double curvature_epsilon);

        /// Get the 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        double get_curvature_epsilon();

        /// Free the instance.
        void free();

      protected:
        /// The 'curvature' epsilon.
        double curvature_epsilon;

        /// Information if user-supplied displacement functions have been provided.
        bool user_xdisp, user_ydisp;

        /// Displacement functions, default to ZeroFunctions, may be supplied by set_displacement();
        MeshFunction<double> *xdisp, *ydisp;
        double dmult;

        int3* tris_contours;      ///< triangles: vertex index triplets
        int triangle_contours_count;
        double3* verts;  ///< vertices: (x, y, value) triplets

        /// What kind of information do we want to get out of the solution.
        int item, component, value_type;

        int add_vertex();
        int get_vertex(int p1, int p2, double x, double y, double value);

        void process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
          double* val, double* phx, double* phy, int* indices, bool curved);

        void process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
          double* val, double* phx, double* phy, int* indices, bool curved);

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
