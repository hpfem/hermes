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

#include "../global.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /// Index of a function expansion. Used to selected a value in Shapeset::get_value().
    enum FunctionExpansionIndex {
      H2D_FEI_VALUE = 0, ///< Index of a function value f.
      H2D_FEI_DX = 1, ///< Index of df/dx.
      H2D_FEI_DY = 2, ///< Index of df/dy.
      H2D_FEI_DXX = 3, ///< Index of df/dxdx.
      H2D_FEI_DYY = 4, ///< Index of df/dydy.
      H2D_FEI_DXY = 5 ///< Index of df/dxdy.
    };

    namespace RefinementSelectors
    {
      template<typename Scalar> class Selector;
      template<typename Scalar> class POnlySelector;
      template<typename Scalar> class HOnlySelector;
      template<typename Scalar> class OptimumSelector;
      template<typename Scalar> class ProjBasedSelector;
      template<typename Scalar> class H1ProjBasedSelector;
      template<typename Scalar> class L2ProjBasedSelector;
      template<typename Scalar> class HcurlProjBasedSelector;
    };

    /// @ingroup spaces
    /// \brief Defines a set of shape functions.
    ///
    /// This class stores mainly the definitions of the polynomials for all shape functions,
    /// but also their polynomial degrees, their types (vertex, edge, bubble), and contains
    /// mechanisms for the calculation and storage of constrained shape functions.
    ///
    /// The class returns shape function values for both triangles and quads, depending on
    /// what mode is requested.
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
    /// The polynomial degree (or 'order') is an integer typically in the range[1-10] for H1
    /// shapesets and[0-10] for H(curl) shapesets. Quadrilaterals are allowed to have different
    /// orders in the x and y directions (of the reference domain). The 'order' for quads thus
    /// has to be formed with the macro H2D_MAKE_QUAD_ORDER(), see h2d_common.h.
    ///
    /// Vertex shape functions in H1 shapesets are also regarded as edge functions of orders 0
    /// and 1. This simplifies constraint calculations and BC projections.
    ///
    /// Shape functions are always Real-valued.
    ///
    class HERMES_API Shapeset : public Hermes::Mixins::Loggable
    {
    public:
      ~Shapeset();

      /// Shape-function function type. Internal.
      typedef double (*shape_fn_t)(double, double);

      /// Returns the polynomial degree of the specified shape function.
      /// If on quads, it returns encoded orders. The orders has to be decoded through macros
      /// H2D_GET_H_ORDER and H2D_GET_V_ORDER.
      int get_order(int index, ElementMode2D mode) const;

      virtual Shapeset* clone() = 0;

      /// Returns 2 if this is a vector shapeset, 1 otherwise.
      int get_num_components() const;

      /// Returns the maximum poly degree for all shape functions.
      int get_max_order() const;
      int get_min_order() const;

      /// Returns the highest shape function index.
      virtual int get_max_index(ElementMode2D mode) = 0;

      /// Returns the index of a vertex shape function associated with the specified vertex.
      int get_vertex_index(int vertex, ElementMode2D mode) const;

      /// Returns the index of an edge function associated with the specified edge and of the
      /// requested order. 'ori' can be 0 or 1 and determines edge orientation (this is for
      /// shapesets with non-symmetric edge functions).
      int get_edge_index(int edge, int ori, int order, ElementMode2D mode) const;

      /// Returns space type.
      virtual SpaceType get_space_type() const = 0;

    protected:
      /// Returns a complete set of indices of bubble functions for an element of the given order.
      int* get_bubble_indices(int order, ElementMode2D mode) const;

      /// Returns the number of bubble functions for an element of the given order.
      int get_num_bubbles(int order, ElementMode2D mode) const;

      /// Returns the index of a constrained edge function. 'part' is 0 or 1 for edge
      /// halves, 2, 3, 4, 5 for edge quarters, etc. See shapeset.cpp.
      int get_constrained_edge_index(int edge, int order, int ori, int part, ElementMode2D mode) const;

      /// Obtains the value of the given shape function. (x,y) is a coordinate in the reference
      /// domain, component is 0 for Scalar shapesets and 0 or 1 for vector shapesets.
      double get_value(int n, int index, double x, double y, int component, ElementMode2D mode);

      double get_fn_value (int index, double x, double y, int component, ElementMode2D mode);
      double get_dx_value (int index, double x, double y, int component, ElementMode2D mode);
      double get_dy_value (int index, double x, double y, int component, ElementMode2D mode);
      double get_dxx_value(int index, double x, double y, int component, ElementMode2D mode);
      double get_dyy_value(int index, double x, double y, int component, ElementMode2D mode);
      double get_dxy_value(int index, double x, double y, int component, ElementMode2D mode);

      /// Returns the coordinates of the reference domain vertices.
      double2* get_ref_vertex(int vertex, ElementMode2D mode);

      /// Returns shapeset identifier. Internal.
      virtual int get_id() const = 0;

      shape_fn_t*** shape_table[6];

      int**  vertex_indices;
      int*** edge_indices;
      int*** bubble_indices;
      int**  bubble_count;
      int**  index_to_order;

      double2 ref_vert[H2D_MAX_SOLUTION_COMPONENTS][H2D_MAX_NUMBER_VERTICES];
      int max_order, min_order;
      int num_components;

      int ebias; ///< 2 for H1 shapesets, 0 for H(curl) shapesets. It is the order of the
      ///< first edge function.

      double** comb_table;
      int table_size;
      /**    numbering of edge intervals: (the variable 'part')
      -+-        -+-         -+-
      |          |        13 |
      |        5 |          -+-
      |          |        12 |
      1 |         -+-         -+-                  finer interval:
      |          |        11 |
      |        4 |          -+-                  p = (p + 1) * 2   (+1)
      |          |        10 |
      -+-        -+-         -+-   ... etc.
      |          |         9 |
      |        3 |          -+-
      |          |         8 |
      0 |         -+-         -+-
      |          |         7 |
      |        2 |          -+-
      |          |         6 |
      -+-        -+-         -+-            **/

      /// Constrained edge functions are constructed by subtracting the linear part (ie., two
      /// vertex functions) from the constraining edge function and expressing the rest as a
      /// linear combination of standard edge functions. This function determines the coefficients
      /// of such linear combination by forming and solving a simple linear system.
      ///
      double* calculate_constrained_edge_combination(int order, int part, int ori, ElementMode2D mode);

      /// Returns the coefficients for the linear combination forming a constrained edge function.
      /// This function performs the storage (caching) of these coefficients, so that they can be
      /// calculated only once.
      ///
      double* get_constrained_edge_combination(int order, int part, int ori, int& nitems, ElementMode2D mode);

      /// Releases all cached coefficients.
      void free_constrained_edge_combinations();

      /// Constructs the linear combination of edge functions, forming a constrained edge function.
      ///
      double get_constrained_value(int n, int index, double x, double y, int component, ElementMode2D mode);

      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class Solution;
      friend class CurvMap; friend class RefMap;
      template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::OptimumSelector;
      friend class PrecalcShapeset;
      friend void check_leg_tri(Shapeset* shapeset);
      friend void check_gradleg_tri(Shapeset* shapeset);
      template<typename Scalar> friend class Form;
      template<typename Scalar> friend class MatrixForm;
      template<typename Scalar> friend class VectorForm;
      template<typename Scalar> friend class Space;
      template<typename Scalar> friend class H1Space;
      template<typename Scalar> friend class L2Space;
      template<typename Scalar> friend class HcurlSpace;
      template<typename Scalar> friend class HdivSpace;
      template<typename Scalar> friend class L2MaterialWiseConstSpace;
    };
  }
}
#endif