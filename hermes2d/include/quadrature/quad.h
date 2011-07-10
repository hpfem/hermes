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

#ifndef __H2D_QUAD_H
#define __H2D_QUAD_H

#include "../hermes2d_common_defs.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /// Indices of values in the value returned by Quad2D::get_points().
    enum GIP2DIndices {
      H2D_GIP2D_X = 0, ///< X-axis coordinate.
      H2D_GIP2D_Y = 1, ///< Y-axis coordinate.
      H2D_GIP2D_W = 2  ///< A weight.
    };

#define H2D_GIP1D_X 0
#define H2D_GIP1D_W 1

#ifdef EXTREME_QUAD
    const int g_max_quad = 99;
#else
    const int g_max_quad = 24;
    const int g_max_tri = 20;
#endif

    /// Quad1D is a base class for all 1D quadrature points.
    ///
    class HERMES_API Quad1D
    {
    public:

      double2* get_points(int order) const { return tables[order]; }
      int get_num_points(int order) const { return np[order]; };

      int get_max_order() const { return max_order; }
      double get_ref_vertex(int n) const { return ref_vert[n]; }

    protected:

      double2** tables;
      int* np;

      double ref_vert[2];
      int max_order;

      virtual void dummy_fn() = 0; // to prevent this class from being instantiated

    };


    /// Quad2D is a base class for all 2D quadrature points on triangles and quads.
    ///
    class HERMES_API Quad2D
    {
    public:
      void set_mode(int mode) { this->mode = mode; }
      int  get_mode() const { return mode; }

      int get_num_points(int order)  const { return np[mode][order]; };
      double3* get_points(int order) const { assert(order < num_tables[mode]); return tables[mode][order]; }
      int get_edge_points(int edge)  const { return max_order[mode]+1 + (3*(1-mode) + 4*mode)*max_order[mode] + edge; }
      int get_edge_points(int edge, int order) {return  max_order[mode]+1 + (3*(1-mode) + 4*mode)*order + edge;}

      int get_max_order() const { return max_order[mode]; }
      int get_safe_max_order() const { return safe_max_order[mode]; }
      int get_num_tables() const { return num_tables[mode]; }

      double2* get_ref_vertex(int n) { return &ref_vert[mode][n]; }

    protected:

      int mode;

      double3*** tables;
      int** np;

      int num_tables[2];
      int max_order[2], safe_max_order[2];
      int max_edge_order;

      double2 ref_vert[2][4];
    };
  }
}
#endif
