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

        void process_solution(MeshFunction<double>* sln, int item = H2D_FN_VAL_0, double eps = HERMES_EPS_NORMAL);

        /// Saves a MeshFunction (Solution, Filter) in VTK format.
        void save_solution_vtk(MeshFunction<double>* sln, const char* filename, const char* quantity_name,
          bool mode_3D = true, int item = H2D_FN_VAL_0,
          double eps = HERMES_EPS_NORMAL);

        void set_displacement(MeshFunction<double>* xdisp, MeshFunction<double>* ydisp, double dmult = 1.0);

        void calc_vertices_aabb(double* min_x, double* max_x,
          double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        int get_num_vertices();
        double3* get_vertices();
      protected:
        void free();

        MeshFunction<double>* sln;

        double cmax;

        /// Information if user-supplied displacement functions have been provided.
        bool user_xdisp, user_ydisp;
        /// Displacement functions, default to ZeroFunctions, may be supplied by set_displacement();
        MeshFunction<double> *xdisp, *ydisp;
        double dmult;

        double3* verts;  ///< vertices: (x, y, value) triplets

        /// What kind of information do we want to get out of the solution.
        int item, component, value_type;

        int del_slot;   ///< free slot index after a triangle which was deleted

        int add_vertex();
        int get_vertex(int p1, int p2, double x, double y, double value);
        int get_top_vertex(int id, double value);

        void process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
          double* val, double* phx, double* phy, int* indices);

        void process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
          double* val, double* phx, double* phy, int* indices);

        void regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2);

        void find_min_max();
      };
    }
  }
}
#endif
