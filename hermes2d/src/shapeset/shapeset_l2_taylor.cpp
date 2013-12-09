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

#include "global.h"
#include "shapeset.h"
#include "shapeset_common.h"
#include "shapeset_l2_all.h"

namespace Hermes
{
  namespace Hermes2D
  {

#define minus_triangle_x_c 0.3333333333333333333333333
#define minus_triangle_y_c 0.3333333333333333333333333

#pragma region taylor_0
    static double taylor_fn_0(double x, double y)
    {
      return 1.;
    }

    static double taylor_dx_0(double x, double y)
    {
      return 0.;
    }

    static double taylor_dy_0(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxx_0(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_0(double x, double y)
    {
      return 0.;
    }

    static double taylor_dyy_0(double x, double y)
    {
      return 0.;
    }
#pragma endregion

#pragma region taylor_1
    static double taylor_fn_1_tri(double x, double y)
    {
      return (x + minus_triangle_x_c) / ELEMENT_DELTA_X;
    }

    static double taylor_fn_1_quad(double x, double y)
    {
      return x / ELEMENT_DELTA_X;
    }

    static double taylor_dx_1(double x, double y)
    {
      return 1. / ELEMENT_DELTA_X;
    }

    static double taylor_dy_1(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxx_1(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_1(double x, double y)
    {
      return 0.;
    }

    static double taylor_dyy_1(double x, double y)
    {
      return 0.;
    }
#pragma endregion

#pragma region taylor_2
    static double taylor_fn_2_tri(double x, double y)
    {
      return (y + minus_triangle_y_c) / ELEMENT_DELTA_Y;
    }

    static double taylor_fn_2_quad(double x, double y)
    {
      return y / ELEMENT_DELTA_Y;
    }

    static double taylor_dx_2(double x, double y)
    {
      return 0.;
    }

    static double taylor_dy_2(double x, double y)
    {
      return 1. / ELEMENT_DELTA_Y;
    }

    static double taylor_dxx_2(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_2(double x, double y)
    {
      return 0.;
    }

    static double taylor_dyy_2(double x, double y)
    {
      return 0.;
    }
#pragma endregion

#pragma region taylor_3
    static double taylor_fn_3_tri(double x, double y)
    {
      // Mean value: 0.5 * 0.125 * Integrate[(x+c)*(x+c), {y, -1, 1}, {x, -1, -y}] = 2/3 (1 - 2c + 3c^2)
      // = 0.25 / 9.
      return ((x + minus_triangle_x_c) * (x + minus_triangle_x_c) / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_X)) - (0.25 / 9.);
    }

    static double taylor_fn_3_quad(double x, double y)
    {
      // (1/24) is the mean value.
      return (x * x / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_X)) - (1. / 24.);
    }

    static double taylor_dx_3_tri(double x, double y)
    {
      return ((2. * x) + (2. * minus_triangle_x_c)) / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_X);
    }

    static double taylor_dx_3_quad(double x, double y)
    {
      return (2. * x) / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_X);
    }

    static double taylor_dy_3(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxx_3(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_3(double x, double y)
    {
      return 0.;
    }

    static double taylor_dyy_3(double x, double y)
    {
      return 0.;
    }
#pragma endregion

#pragma region taylor_4
    static double taylor_fn_4_tri(double x, double y)
    {
      // Mean value: 0.5 * 0.125 * Integrate[(x+c)*(x+c), {y, -1, 1}, {x, -1, -y}] = 2/3 (1 - 2c + 3c^2)
      // = 0.25 / 9.
      return ((y + minus_triangle_y_c) * (y + minus_triangle_y_c) / (2. * ELEMENT_DELTA_Y * ELEMENT_DELTA_Y)) - (0.25 / 9.);
    }

    static double taylor_fn_4_quad(double x, double y)
    {
      // (1/24) is the mean value.
      return (y * y / (2. * ELEMENT_DELTA_Y * ELEMENT_DELTA_Y)) - (1. / 24.);
    }

    static double taylor_dx_4(double x, double y)
    {
      return 0.;
    }

    static double taylor_dy_4_tri(double x, double y)
    {
      return ((2. * y) + (2. * minus_triangle_y_c)) / (2. * ELEMENT_DELTA_Y * ELEMENT_DELTA_Y);
    }

    static double taylor_dy_4_quad(double x, double y)
    {
      return (2. * y) / (2. * ELEMENT_DELTA_Y * ELEMENT_DELTA_Y);
    }

    static double taylor_dxx_4(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_4(double x, double y)
    {
      return 0.;
    }

    static double taylor_dyy_4(double x, double y)
    {
      return 0.;
    }
#pragma endregion

#pragma region taylor_5
    static double taylor_fn_5_tri(double x, double y)
    {
      // Mean value: 0.5 * 0.125 * Integrate[(x+c)*(y+c), {y, -1, 1}, {x, -1, -y}] = 2/3 (c * (-2 + 3c))
      // = - 0.125 / 9.
      return ((x + minus_triangle_x_c) * (y + minus_triangle_y_c) / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_Y)) + (0.125 / 9.);
    }

    static double taylor_fn_5_quad(double x, double y)
    {
      // (0.) is the mean value.
      return (x * y / (2. * ELEMENT_DELTA_X * ELEMENT_DELTA_Y));
    }

    static double taylor_dx_5_tri(double x, double y)
    {
      return (y + minus_triangle_y_c) / (2 * ELEMENT_DELTA_X * ELEMENT_DELTA_Y);
    }

    static double taylor_dx_5_quad(double x, double y)
    {
      return y / (2 * ELEMENT_DELTA_X * ELEMENT_DELTA_Y);
    }

    static double taylor_dy_5_tri(double x, double y)
    {
      return (x + minus_triangle_x_c) / (2 * ELEMENT_DELTA_X * ELEMENT_DELTA_Y);
    }

    static double taylor_dy_5_quad(double x, double y)
    {
      return x / (2 * ELEMENT_DELTA_X * ELEMENT_DELTA_Y);
    }

    static double taylor_dxx_5(double x, double y)
    {
      return 0.;
    }

    static double taylor_dxy_5(double x, double y)
    {
      return 1. / (2 * ELEMENT_DELTA_X * ELEMENT_DELTA_Y);
    }

    static double taylor_dyy_5(double x, double y)
    {
      return 0.;
    }
#pragma endregion

    static Shapeset::shape_fn_t fn_tri[] =
    {
      taylor_fn_0, taylor_fn_1_tri, taylor_fn_2_tri, taylor_fn_3_tri, taylor_fn_4_tri, taylor_fn_5_tri
    };

    static Shapeset::shape_fn_t fn_quad[] =
    {
      taylor_fn_0, taylor_fn_1_quad, taylor_fn_2_quad, taylor_fn_3_quad, taylor_fn_4_quad, taylor_fn_5_quad
    };

    static Shapeset::shape_fn_t dx_tri[] =
    {
      taylor_dx_0, taylor_dx_1, taylor_dx_2, taylor_dx_3_tri, taylor_dx_4, taylor_dx_5_tri
    };

    static Shapeset::shape_fn_t dx_quad[] =
    {
      taylor_dx_0, taylor_dx_1, taylor_dx_2, taylor_dx_3_quad, taylor_dx_4, taylor_dx_5_quad
    };

    static Shapeset::shape_fn_t dy_tri[] =
    {
      taylor_dy_0, taylor_dy_1, taylor_dy_2, taylor_dy_3, taylor_dy_4_tri, taylor_dy_5_tri
    };

    static Shapeset::shape_fn_t dy_quad[] =
    {
      taylor_dy_0, taylor_dy_1, taylor_dy_2, taylor_dy_3, taylor_dy_4_quad, taylor_dy_5_quad
    };

    static Shapeset::shape_fn_t fn_dxx[] =
    {
      taylor_dxx_0, taylor_dxx_1, taylor_dxx_2, taylor_dxx_3, taylor_dxx_4, taylor_dxx_5
    };

    static Shapeset::shape_fn_t fn_dxy[] =
    {
      taylor_dxy_0, taylor_dxy_1, taylor_dxy_2, taylor_dxy_3, taylor_dxy_4, taylor_dxy_5
    };

    static Shapeset::shape_fn_t fn_dyy[] =
    {
      taylor_dyy_0, taylor_dyy_1, taylor_dyy_2, taylor_dyy_3, taylor_dyy_4, taylor_dyy_5
    };

    Shapeset::shape_fn_t* shape_fn_table_tri[1]     = { fn_tri };
    Shapeset::shape_fn_t* shape_fn_table_quad[1]     = { fn_quad };
    Shapeset::shape_fn_t* shape_fn_table_dx_tri[1]  = { dx_tri };
    Shapeset::shape_fn_t* shape_fn_table_dx_quad[1]  = { dx_quad };
    Shapeset::shape_fn_t* shape_fn_table_dy_tri[1]  = { dy_tri };
    Shapeset::shape_fn_t* shape_fn_table_dy_quad[1]  = { dy_quad };
    Shapeset::shape_fn_t* mode_shape_fn_table_dxx[1]  = { fn_dxx };
    Shapeset::shape_fn_t* mode_shape_fn_table_dxy[1]  = { fn_dxy };
    Shapeset::shape_fn_t* mode_shape_fn_table_dyy[1]  = { fn_dyy };

    static int qb_0[] = { 0, };
    static int qb_1_contains_means[] = { 0, 1, 2, };
    static int qb_1[] = { 1, 2, };
    static int qb_2_contains_means[] = { 0, 1, 2, 3, 4, 5 };
    static int qb_2[] = { 1, 2, 3, 4, 5 };

    int* mode_bubble_indices_contains_means[3] = {  qb_0,   qb_1_contains_means,  qb_2_contains_means };
    int* mode_bubble_indices[3] = {  qb_0,   qb_1,  qb_2 };

    int mode_bubble_count_contains_means[3] = { 1,  3,  6 };
    int mode_bubble_count[3] = { 1,  2,  5 };

    int mode_vertex_indices[4] = { -1, -1, -1, -1 };

    static int edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };
    static int edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };

    int* mode_edge_indices[4] =
    {
      edge_indices_0,
      edge_indices_1,
      edge_indices_2,
      edge_indices_3
    };

    int mode_index_to_order[] =
    {
      0, 1, 1, 2, 2, 2
    };

    static Shapeset::shape_fn_t** shape_fn_table[2] =
    {
      shape_fn_table_tri,
      shape_fn_table_quad
    };

    static Shapeset::shape_fn_t** shape_fn_table_dx[2] =
    {
      shape_fn_table_dx_tri,
      shape_fn_table_dx_quad
    };

    static Shapeset::shape_fn_t** shape_fn_table_dy[2] =
    {
      shape_fn_table_dy_tri,
      shape_fn_table_dy_quad
    };

    static Shapeset::shape_fn_t** shape_fn_table_dxx[2] =
    {
      mode_shape_fn_table_dxx,
      mode_shape_fn_table_dxx
    };

    static Shapeset::shape_fn_t** shape_fn_table_dyy[2] =
    {
      mode_shape_fn_table_dyy,
      mode_shape_fn_table_dyy
    };

    static Shapeset::shape_fn_t** shape_fn_table_dxy[2] =
    {
      mode_shape_fn_table_dxy,
      mode_shape_fn_table_dxy
    };

    static int* s_vertex_indices[2] =
    {
      mode_vertex_indices,
      mode_vertex_indices
    };

    static int** s_edge_indices[2] =
    {
      mode_edge_indices,
      mode_edge_indices
    };

    static int** s_bubble_indices_contains_means[2] =
    {
      mode_bubble_indices_contains_means,
      mode_bubble_indices_contains_means
    };

    static int** s_bubble_indices[2] =
    {
      mode_bubble_indices,
      mode_bubble_indices
    };

    static int* s_bubble_count_contains_means[2] =
    {
      mode_bubble_count_contains_means,
      mode_bubble_count_contains_means
    };

    static int* s_bubble_count[2] =
    {
      mode_bubble_count,
      mode_bubble_count
    };

    static int* s_index_to_order[2] =
    {
      mode_index_to_order,
      mode_index_to_order
    };

    L2ShapesetTaylor::L2ShapesetTaylor(bool contains_means)
    {
      shape_table[0] = shape_fn_table;
      shape_table[1] = shape_fn_table_dx;
      shape_table[2] = shape_fn_table_dy;
      shape_table[3] = shape_fn_table_dxx;
      shape_table[4] = shape_fn_table_dyy;
      shape_table[5] = shape_fn_table_dxy;

      vertex_indices = s_vertex_indices;
      edge_indices = s_edge_indices;
      if(contains_means)
        bubble_indices = s_bubble_indices_contains_means;
      else
        bubble_indices = s_bubble_indices;

      if(contains_means)
        bubble_count = s_bubble_count_contains_means;
      else
        bubble_count = s_bubble_count;

      index_to_order = s_index_to_order;

      ref_vert[0][0][0] = -1.0;
      ref_vert[0][0][1] = -1.0;
      ref_vert[0][1][0] =  1.0;
      ref_vert[0][1][1] = -1.0;
      ref_vert[0][2][0] = -1.0;
      ref_vert[0][2][1] =  1.0;

      ref_vert[1][0][0] = -1.0;
      ref_vert[1][0][1] = -1.0;
      ref_vert[1][1][0] =  1.0;
      ref_vert[1][1][1] = -1.0;
      ref_vert[1][2][0] =  1.0;
      ref_vert[1][2][1] =  1.0;
      ref_vert[1][3][0] = -1.0;
      ref_vert[1][3][1] =  1.0;

      max_order = 2;
      min_order = 0;
      num_components = 1;

      ebias = 0;

      comb_table = NULL;
    }

    int* L2ShapesetTaylor::get_bubble_indices(int order, ElementMode2D mode) const
    {
      if(mode == HERMES_MODE_QUAD)
      {
        assert(H2D_GET_V_ORDER(order) == H2D_GET_H_ORDER(order));
        return bubble_indices[mode][H2D_GET_V_ORDER(order)];
      }
      else
        return Shapeset::get_bubble_indices(order, mode);
    }

    int L2ShapesetTaylor::get_num_bubbles(int order, ElementMode2D mode) const
    {
      if(mode == HERMES_MODE_QUAD)
      {
        assert(H2D_GET_V_ORDER(order) == H2D_GET_H_ORDER(order));
        return bubble_count[mode][H2D_GET_V_ORDER(order)];
      }
      else
        return Shapeset::get_num_bubbles(order, mode);
    }

    const int L2ShapesetTaylor::max_index[2] = { 5, 5 };
    int L2ShapesetTaylor::get_max_index(ElementMode2D mode) { return max_index[mode]; }
  }
}