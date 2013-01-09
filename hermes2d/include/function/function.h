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

#ifndef __H2D_FUNCTION_H
#define __H2D_FUNCTION_H

#include "transformable.h"
#include "../quadrature/quad.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    struct SurfPos;
    class PrecalcShapeset;
    namespace RefinementSelectors{
      template<typename Scalar> class Selector;
      template<typename Scalar> class HOnlySelector;
      template<typename Scalar> class POnlySelector;
      template<typename Scalar> class OptimumSelector;
      template<typename Scalar> class ProjBasedSelector;
      template<typename Scalar> class L2ProjBasedSelector;
      template<typename Scalar> class H1ProjBasedSelector;
      template<typename Scalar> class HcurlProjBasedSelector;
    };
    template<typename Scalar> class Geom;

    namespace Views{
      class Orderizer;
      class Linearizer;
      class Vectorizer;
    };

    /// Precalculation masks
    enum
    {
      H2D_FN_VAL_0 = 0x0001, H2D_FN_VAL_1 = 0x0040, // Function values
      H2D_FN_DX_0  = 0x0002, H2D_FN_DX_1  = 0x0080, // First derivative
      H2D_FN_DY_0  = 0x0004, H2D_FN_DY_1  = 0x0100, // First derivative
      H2D_FN_DXX_0 = 0x0008, H2D_FN_DXX_1 = 0x0200, // Second derivative
      H2D_FN_DYY_0 = 0x0010, H2D_FN_DYY_1 = 0x0400, // Second derivative
      H2D_FN_DXY_0 = 0x0020, H2D_FN_DXY_1 = 0x0800  // Second mixed derivative
    };

    /// Both components are usually requested together...
    const int H2D_FN_VAL = H2D_FN_VAL_0 | H2D_FN_VAL_1;
    const int H2D_FN_DX  = H2D_FN_DX_0  | H2D_FN_DX_1;
    const int H2D_FN_DY  = H2D_FN_DY_0  | H2D_FN_DY_1;
    const int H2D_FN_DXX = H2D_FN_DXX_0 | H2D_FN_DXX_1;
    const int H2D_FN_DYY = H2D_FN_DYY_0 | H2D_FN_DYY_1;
    const int H2D_FN_DXY = H2D_FN_DXY_0 | H2D_FN_DXY_1;

    const int H2D_FN_DEFAULT = H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY;            ///< default precalculation mask
    const int H2D_FN_ALL = H2D_FN_DEFAULT | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY; ///< precalculate everything

    const int H2D_FN_COMPONENT_0 = H2D_FN_VAL_0 | H2D_FN_DX_0 | H2D_FN_DY_0 | H2D_FN_DXX_0 | H2D_FN_DYY_0 | H2D_FN_DXY_0;
    const int H2D_FN_COMPONENT_1 = H2D_FN_VAL_1 | H2D_FN_DX_1 | H2D_FN_DY_1 | H2D_FN_DXX_1 | H2D_FN_DYY_1 | H2D_FN_DXY_1;

    /// @ingroup meshFunctions
    /// \brief Represents an arbitrary function defined on an element.
    ///
    /// The Function class is an abstraction of a function defined in integration points on an
    /// element. You first specify what quadrature tables you want to use (set_quad_2d()) and select
    /// an element (Transformable::set_active_element()). Then you select concrete integration points
    /// (set_quad_order()) and obtain the function values by calling one of the functions get_fn_values(),
    /// get_dx_values(), etc.
    ///
    /// This class is a template for RealFunction and ScalarFunction, depending of which type the
    /// function values are. For example, shape functions are always Real (see PrecalcShapeset), while
    /// the solution can be complex (see Solution).
    ///
    /// The design goal for this class is to define a single common interface for functions used as
    /// integrands in the weak formulation. It should not matter whether you are integrating a shape
    /// function or, for example, a previous solution of the PDE in time-dependent problems.
    /// Ideally, you should also be able to apply the bilinear form not only to shape functions
    /// during assembling, but also to the solution when calculating energy norms etc. The last
    /// feature is unfortunately limited to Real code, because a PDE solution can be complex (hence
    /// Solution inherits from ScalarFunction), but shape functions are Real and for efficiency
    /// the bilinear form only takes RealFunction arguments.
    ///
    /// Since this class inherits from Transformable, you can obtain function values in integration
    /// points transformed to sub-areas of the current element (see push_transform(), pop_transform()).
    ///
    template<typename Scalar>
    class HERMES_API Function : public Transformable
    {
    public:

      /// Default constructor.
      Function();

      /// Default destructor.
      /// Added to delete Function::sub_tables.
      virtual ~Function();

      /// \brief Returns the number of components of the function being represented by the class.
      int get_num_components() const;

      /// \brief Returns function values.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The values of the function at all points of the current integration rule.
      Scalar* get_fn_values(int component = 0);

      /// \brief Returns the x partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The x partial derivative of the function at all points of the current integration rule.
      Scalar* get_dx_values(int component = 0);

      /// \brief Returns the y partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The y partial derivative of the function at all points of the current integration rule.
      Scalar* get_dy_values(int component = 0);

      /// \brief Returns both x and y partial derivatives.
      /// This function provides the both often-used dx and dy values in one call.
      /// \param dx[out] Variable which receives the pointer to the first partial derivatives by x
      /// \param dy[out] Variable which receives the pointer to the first partial derivatives by y
      /// \param component[in] The component of the function (0 or 1).
      void get_dx_dy_values(Scalar*& dx, Scalar*& dy, int component = 0);

      /// \brief Returns the second x partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The x second partial derivative of the function at all points of the current integration rule.
      Scalar* get_dxx_values(int component = 0);

      /// \brief Returns the second y partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The y second partial derivative of the function at all points of the current integration rule.
      Scalar* get_dyy_values(int component = 0);

      /// \brief Returns the second mixed derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The second mixed derivative of the function at all points of the current integration rule.
      Scalar* get_dxy_values(int component = 0);

      /// \brief Returns the current quadrature points.
      Quad2D* get_quad_2d() const;

      /// Activates an integration rule of the specified order. Subsequent calls to
      /// get_values(), get_dx_values() etc. will be returning function values at these points.
      /// \param order[in] Integration rule order.
      /// \param mask[in] A combination of one or more of the constants H2D_FN_VAL, H2D_FN_DX, H2D_FN_DY,
      ///   H2D_FN_DXX, H2D_FN_DYY, H2D_FN_DXY specifying the values which should be precalculated. The default is
      ///   H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY. You can also use H2D_FN_ALL to precalculate everything.
      void set_quad_order(unsigned int order, int mask = H2D_FN_DEFAULT);

      Scalar* get_values(int a, int b);

      /// \brief Returns the polynomial degree of the function being represented by the class.
      int get_fn_order() const;

    protected:
      /// \brief Selects the quadrature points in which the function will be evaluated.
      /// \details It is possible to switch back and forth between different quadrature
      /// points: no precalculated values are freed. The standard quadrature is
      /// always selected by default already.
      /// \param quad_2d[in] The quadrature points.
      virtual void set_quad_2d(Quad2D* quad_2d);

      /// \brief Frees all precalculated tables.
      virtual void free() = 0;

      struct Node
      {
        int mask;           ///< a combination of H2D_FN_XXX: specifies which tables are present

        int size;           ///< size in bytes of this struct (for maintaining total_mem)

        Scalar* values[2][6]; ///< pointers to 'data'

        Scalar data[1];       ///< value tables. The length may vary.

      private:
        Node(const Node& org) {}; ///< Copy constructor is disabled.

        Node& operator=(const Node& other) { return *this; }; ///< Assignment is not allowed.
      };

      /// \brief Returns the polynomial degree of the function at given edge. To be overridden in derived classes.
      /// \param edge[in] Edge at which the order should be evaluated. (0-3)
      virtual int get_edge_fn_order(int edge) const;

      /// precalculates the current function at the current integration points.
      virtual void precalculate(int order, int mask) = 0;

      int order;          ///< current function polynomial order

      int num_components; ///< number of vector components

      /// Table of Node tables, for each possible transformation there can be a different Node table.
      std::map<uint64_t, LightArray<Node*>*>* sub_tables;

      /// Table of nodes.
      LightArray<Node*>* nodes;

      /// Current Node.
      Node* cur_node;

      /// Nodes for the overflow sub-element transformation.
      LightArray<Node*>* overflow_nodes;

      /// With changed sub-element mapping, there comes the need for a change of the current
      /// Node table nodes.
      void update_nodes_ptr();

      /// For internal use only.
      void force_transform(uint64_t sub_idx, Trf* ctm);

      Quad2D* quads[8]; ///< list of available quadratures

      int cur_quad;     ///< active quadrature (index into 'quads')

      int total_mem;    ///< total memory in bytes used by the tables

      int max_mem;      ///< peak memory usage

      Node* new_node(int mask, int num_points); ///< allocates a new Node structure

      virtual void  handle_overflow_idx() = 0;

      void replace_cur_node(Node* node);

      static void check_params(int component, Node* cur_node, int num_components);

      static void check_table(int component, Node* cur_node, int n, const char* msg);

      static int idx2mask[6][2];  ///< index to mask table
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::HcurlProjBasedSelector;
      template<typename T> friend class Adapt;
      friend class Views::Orderizer;
      friend class Views::Vectorizer;
      friend class Views::Linearizer;

      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class Global;
      friend class CurvMap;

      template<typename T> friend class Func;
      template<typename T> friend class Geom;
      template<typename T> friend class Filter;
      template<typename T> friend class SimpleFilter;
      template<typename T> friend class DXDYFilter;
      friend class ComplexFilter;
      friend class VonMisesFilter;
      friend HERMES_API Geom<double>* init_geom_vol(RefMap *rm, const int order);
      friend HERMES_API Geom<double>* init_geom_surf(RefMap *rm, SurfPos* surf_pos, const int order);
      friend HERMES_API Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
      template<typename T> friend HERMES_API Func<T>* init_fn(MeshFunction<T>*fu, const int order);
    };
  }
}
#endif
