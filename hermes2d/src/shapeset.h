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

#ifndef __H2D_SHAPESET_H
#define __H2D_SHAPESET_H

#include "common.h"


#define H2D_CHECK_MODE      assert(mode == H2D_MODE_TRIANGLE || mode == H2D_MODE_QUAD)
#define H2D_CHECK_VERTEX    assert(vertex >= 0 && vertex < nvert)
#define H2D_CHECK_EDGE      assert(edge >= 0 && edge < nvert)
#define H2D_CHECK_ORDER(o)  assert((o) >= 0 && (o) <= max_order)
#define H2D_CHECK_PART      assert(part >= 0)
#define H2D_CHECK_INDEX     assert(index >= 0 && index <= max_index[mode])
#define H2D_CHECK_COMPONENT assert(component >= 0 && component < num_components)

/// Index of a function expansion. Used to selected a value in Shapeset::get_value().
enum FunctionExpansionIndex {
  H2D_FEI_VALUE = 0, ///< Index of a function value f.
  H2D_FEI_DX = 1, ///< Index of df/dx.
  H2D_FEI_DY = 2, ///< Index of df/dy.
  H2D_FEI_DXX = 3, ///< Index of df/dxdx.
  H2D_FEI_DYY = 4, ///< Index of df/dydy.
  H2D_FEI_DXY = 5 ///< Index of df/dxdy.
};

/// \brief Defines a set of shape functions.
///
/// This class stores mainly the definitions of the polynomials for all shape functions,
/// but also their polynomial degrees, their types (vertex, edge, bubble), and contains
/// mechanisms for the calculation and storage of constrained shape functions.
///
/// The class returns shape function values for both triangles and quads, depending on
/// what mode it is in. Use the function set_mode() to switch between H2D_MODE_TRIANGLE and
/// H2D_MODE_QUAD.
///
/// Each shape function is assigned a unique number - 'index'. For standard shape functions,
/// index is positive and not greater than the value returned by get_max_index(). Negative
/// index values are reserved for constrained shape functions. These are special edge
/// functions designed to fit standard edge functions along a portion of an edge, used in
/// the construction of FE spaces on meshes with hanging nodes.
///
/// The usage of Shapeset is simple: you first obtain the index of the shape function you
/// are interested in, be it a vertex function for a given vertex, or an edge function of
/// a given order, etc. Then you call one of the functions get_fn_value(), get_dx_value(), etc.
/// All shape functions are defined on the reference domain. For triangles, this is the
/// standard triangle (-1,-1), (1,-1), (-1,1), and for quads this is the square (-1,1)^2.
///
/// The polynomial degree (or 'order') is an integer typically in the range [1-10] for H1
/// shapesets and [0-10] for H(curl) shapesets. Quadrilaterals are allowed to have different
/// orders in the x and y directions (of the reference domain). The 'order' for quads thus
/// has to be formed with the macro H2D_MAKE_QUAD_ORDER(), see common.h.
///
/// Vertex shape functions in H1 shapesets are also regarded as edge functions of orders 0
/// and 1. This simplifies constraint calculations and BC projections.
///
/// Shape functions are always real-valued.
///
class H2D_API Shapeset
{
public:

  ~Shapeset() { free_constrained_edge_combinations(); }

  /// Selects H2D_MODE_TRIANGLE or H2D_MODE_QUAD.
  void set_mode(int mode)
  {
    H2D_CHECK_MODE;
    this->mode = mode;
    nvert = (mode == H2D_MODE_TRIANGLE) ? 3 : 4;
  }

  /// Returns the current mode.
  int get_mode() const { return mode; }

  /// Returns the maximum poly degree for all shape functions.
  int get_max_order() const { return max_order; }

  /// Returns the highest shape function index.
  int get_max_index() const { return max_index[mode]; }

  /// Returns 2 if this is a vector shapeset, 1 otherwise.
  int get_num_components() const { return num_components; }


  /// Returns the index of a vertex shape function associated with the specified vertex.
  int get_vertex_index(int vertex) const
  {
    H2D_CHECK_VERTEX;
    return vertex_indices[mode][vertex];
  }

  /// Returns the index of an edge function associated with the specified edge and of the
  /// requested order. 'ori' can be 0 or 1 and determines edge orientation (this is for
  /// shapesets with non-symmetric edge functions).
  int get_edge_index(int edge, int ori, int order) const
  {
    H2D_CHECK_EDGE; H2D_CHECK_ORDER(order); assert(ori == 0 || ori == 1);
    return edge_indices[mode][edge][2*order + ori];
  }

  /// Returns a complete set of indices of bubble functions for an element of the given order.
  int* get_bubble_indices(int order) const
  {
    H2D_CHECK_ORDER(H2D_GET_H_ORDER(order));
    H2D_CHECK_ORDER(H2D_GET_V_ORDER(order));
    int index = order;
    if (mode == H2D_MODE_QUAD) //tables of bubble indices are transposed
      index = H2D_MAKE_QUAD_ORDER(H2D_GET_V_ORDER(order), H2D_GET_H_ORDER(order));
    return bubble_indices[mode][index];
  }

  /// Returns the number of bubble functions for an element of the given order.
  int get_num_bubbles(int order) const
  {
    H2D_CHECK_ORDER(H2D_GET_H_ORDER(order));
    H2D_CHECK_ORDER(H2D_GET_V_ORDER(order));
    return bubble_count[mode][order];
  }

  /// Returns the index of a constrained edge function. 'part' is 0 or 1 for edge
  /// halves, 2, 3, 4, 5 for edge quarters, etc. See shapeset.cpp.
  int get_constrained_edge_index(int edge, int order, int ori, int part) const
  {
    H2D_CHECK_EDGE; H2D_CHECK_ORDER(order); H2D_CHECK_PART;
    assert(order <= H2D_ORDER_MASK);
    return -1 - ((part << 7) + (order << 3) + (edge << 1) + ori);
  }

  /// Returns the polynomial degree of the specified shape function.
  /// If on quads, it returns encoded orders. The orders has to be decoded through macros
  /// H2D_GET_H_ORDER and H2D_GET_V_ORDER.
  int get_order(int index) const
  {
    if (index >= 0) {
      H2D_CHECK_INDEX;
      return index_to_order[mode][index];
    }
    else return ((-1 - index) >> 3) & 15;
  }


  /// Obtains the value of the given shape function. (x,y) is a coordinate in the reference
  /// domain, component is 0 for scalar shapesets and 0 or 1 for vector shapesets.
  inline double get_value(int n, int index, double x, double y, int component)
  {
    if (index >= 0)
    {
      H2D_CHECK_INDEX; H2D_CHECK_COMPONENT;
      Shapeset::shape_fn_t** shape_expansion = shape_table[n][mode];
      if (shape_expansion == NULL) { // requested exansion (f, df/dx, df/dy, ddf/dxdx, ...) is not defined
        static int warned_mode = -1, warned_index = -1, warned_n = 1; //just to keep the number of warnings low: warn just once about a given combinations of n, mode, and index.
        warn_if(warned_mode != mode || warned_index != index || warned_n != n, "Requested undefined expansion %d (mode: %d) of a shape %d, returning 0", n, mode, index);
        warned_mode = mode; warned_index = index; warned_n = n;
        return 0;
      }
      else
        return shape_expansion[component][index](x, y);
    }
    else
      return get_constrained_value(n, index, x, y, component);
  }

  inline double get_fn_value (int index, double x, double y, int component) { return get_value(0, index, x, y, component); }
  inline double get_dx_value (int index, double x, double y, int component) { return get_value(1, index, x, y, component); }
  inline double get_dy_value (int index, double x, double y, int component) { return get_value(2, index, x, y, component); }
  inline double get_dxx_value(int index, double x, double y, int component) { return get_value(3, index, x, y, component); }
  inline double get_dyy_value(int index, double x, double y, int component) { return get_value(4, index, x, y, component); }
  inline double get_dxy_value(int index, double x, double y, int component) { return get_value(5, index, x, y, component); }


  /// Returns the coordinates of the reference domain vertices.
  double2* get_ref_vertex(int vertex)
  {
    return &ref_vert[mode][vertex];
  }

  /// Shape-function function type. Internal.
  typedef double (*shape_fn_t)(double, double);

  /// Returns shapeset identifier. Internal.
  virtual int get_id() const = 0;


protected:

  int mode;
  int nvert;

  shape_fn_t*** shape_table[6];

  int**  vertex_indices;
  int*** edge_indices;
  int*** bubble_indices;
  int**  bubble_count;
  int**  index_to_order;

  double2 ref_vert[2][4];
  int max_order;
  int max_index[2];
  int num_components;

  int ebias; ///< 2 for H1 shapesets, 0 for H(curl) shapesets. It is the order of the
             ///< first edge function.

  double** comb_table;
  int table_size;

  double* calculate_constrained_edge_combination(int order, int part, int ori);
  double* get_constrained_edge_combination(int order, int part, int ori, int& nitems);

  void    free_constrained_edge_combinations();

  double get_constrained_value(int n, int index, double x, double y, int component);

};

// TODO : promyslet moznost ulozeni shapesetu jako tabulky monomialnich koeficientu
// - mozna efektivnejsi nez stavajici kod
// - monolitictejsi, elegantnejsi
// - pri ulozeni koeficientu jako long double mozna i presnejsi
// - moznost ulozeni jen zakladnich polynomu, derivace lze dopocitat automaticky


#undef H2D_CHECK_MODE
#undef H2D_CHECK_VERTEX
#undef H2D_CHECK_EDGE
#undef H2D_CHECK_ORDER
#undef H2D_CHECK_PART
#undef H2D_CHECK_INDEX
#undef H2D_CHECK_COMPONENT

#endif
