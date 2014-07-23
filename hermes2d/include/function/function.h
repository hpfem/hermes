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

    /// Precalculation masks
    enum
    {
      H2D_FN_VAL_0 = 0x0001, H2D_FN_VAL_1 = 0x0040, // Function values
      H2D_FN_DX_0 = 0x0002, H2D_FN_DX_1 = 0x0080, // First derivative
      H2D_FN_DY_0 = 0x0004, H2D_FN_DY_1 = 0x0100, // First derivative
#ifdef H2D_USE_SECOND_DERIVATIVES
      H2D_FN_DXX_0 = 0x0008, H2D_FN_DXX_1 = 0x0200, // Second derivative
      H2D_FN_DYY_0 = 0x0010, H2D_FN_DYY_1 = 0x0400, // Second derivative
      H2D_FN_DXY_0 = 0x0020, H2D_FN_DXY_1 = 0x0800  // Second mixed derivative
#endif
    };

    /// Both components are usually requested together...
    const int H2D_FN_VAL = H2D_FN_VAL_0 | H2D_FN_VAL_1;
    const int H2D_FN_DX = H2D_FN_DX_0 | H2D_FN_DX_1;
    const int H2D_FN_DY = H2D_FN_DY_0 | H2D_FN_DY_1;
    /// default precalculation mask
    const int H2D_FN_DEFAULT = H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY;

#ifdef H2D_USE_SECOND_DERIVATIVES
    const int H2D_FN_COMPONENT_0 = H2D_FN_VAL_0 | H2D_FN_DX_0 | H2D_FN_DY_0 | H2D_FN_DXX_0 | H2D_FN_DYY_0 | H2D_FN_DXY_0;
    const int H2D_FN_COMPONENT_1 = H2D_FN_VAL_1 | H2D_FN_DX_1 | H2D_FN_DY_1 | H2D_FN_DXX_1 | H2D_FN_DYY_1 | H2D_FN_DXY_1;
#else
    const int H2D_FN_COMPONENT_0 = H2D_FN_VAL_0 | H2D_FN_DX_0 | H2D_FN_DY_0;
    const int H2D_FN_COMPONENT_1 = H2D_FN_VAL_1 | H2D_FN_DX_1 | H2D_FN_DY_1;
#endif

    const int H2D_GRAD = H2D_FN_DX_0 | H2D_FN_DY_0;
    const int H2D_CURL = H2D_FN_DX | H2D_FN_DY;

#ifdef H2D_USE_SECOND_DERIVATIVES
    const int H2D_FN_DXX = H2D_FN_DXX_0 | H2D_FN_DXX_1;
    const int H2D_FN_DYY = H2D_FN_DYY_0 | H2D_FN_DYY_1;
    const int H2D_FN_DXY = H2D_FN_DXY_0 | H2D_FN_DXY_1;

    const int H2D_SECOND = H2D_FN_DXX_0 | H2D_FN_DXY_0 | H2D_FN_DYY_0;

    /// precalculate everything
    const int H2D_FN_ALL = H2D_FN_DEFAULT | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY;
#endif

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
      unsigned char get_num_components() const;

      /// \brief Returns function values.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The values of the function at all points of the current integration rule.
      virtual const Scalar* get_fn_values(int component = 0) const;

      /// \brief Returns the x partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The x partial derivative of the function at all points of the current integration rule.
      virtual const Scalar* get_dx_values(int component = 0) const;

      /// \brief Returns the y partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The y partial derivative of the function at all points of the current integration rule.
      virtual const Scalar* get_dy_values(int component = 0) const;

#ifdef H2D_USE_SECOND_DERIVATIVES
      /// \brief Returns the second x partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The x second partial derivative of the function at all points of the current integration rule.
      virtual const Scalar* get_dxx_values(int component = 0) const;

      /// \brief Returns the second y partial derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The y second partial derivative of the function at all points of the current integration rule.
      virtual const Scalar* get_dyy_values(int component = 0) const;

      /// \brief Returns the second mixed derivative.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The second mixed derivative of the function at all points of the current integration rule.
      virtual const Scalar* get_dxy_values(int component = 0) const;
#endif

      /// \brief Returns function values.
      /// \param component[in] The component of the function (0 or 1).
      /// \return The values of the function at all points of the current integration rule.
      Scalar* deep_copy_array(int component = 0, int item = 0) const;

      /// \brief Returns the current quadrature points.
      Quad2D* get_quad_2d() const;

      /// Activates an integration rule of the specified order. Subsequent calls to
      /// get_values(), get_dx_values() etc. will be returning function values at these points.
      /// \param order[in] Integration rule order.
      /// \param mask[in] A combination of one or more of the constants H2D_FN_VAL, H2D_FN_DX, H2D_FN_DY, ...
      void set_quad_order(unsigned short order, unsigned short mask = H2D_FN_DEFAULT);

      virtual const Scalar* get_values(int component, int item) const;

      /// \brief Returns the polynomial degree of the function being represented by the class.
      virtual int get_fn_order() const;

      /// See Transformable::push_transform.
      /// Internal.
      virtual void push_transform(int son);

      /// See Transformable::pop_transform.
      /// Internal.
      virtual void pop_transform();

      /// Sets the active element
      virtual void set_active_element(Element* e);

      /// Sets the current transform at once as if it was created by multiple calls to push_transform().
      /// \param idx[in] The number of the sub-element, as returned by get_transform().
      virtual void set_transform(uint64_t idx);

    protected:
      /// \brief Selects the quadrature points in which the function will be evaluated.
      /// \details It is possible to switch back and forth between different quadrature
      /// points: no precalculated values are freed. The standard quadrature is
      /// always selected by default already.
      /// \param quad_2d[in] The quadrature points.
      virtual void set_quad_2d(Quad2D* quad_2d);

      /// Empties the stack, loads identity transform.
      virtual void reset_transform();
      /// For internal use only.
      virtual void force_transform(uint64_t sub_idx, Trf* ctm);

      /// \brief Frees all precalculated tables.
      virtual void free() = 0;

      /// The data.
      Scalar values[H2D_MAX_SOLUTION_COMPONENTS][H2D_NUM_FUNCTION_VALUES][H2D_MAX_INTEGRATION_POINTS_COUNT];
      /// Flag that the data are not 'dirty'
      bool values_valid;

      /// \brief Returns the polynomial degree of the function at given edge. To be overridden in derived classes.
      /// \param edge[in] Edge at which the order should be evaluated. (0-3)
      virtual int get_edge_fn_order(unsigned char edge) const;

      /// precalculates the current function at the current integration points.
      virtual void precalculate(unsigned short order, unsigned short mask);

      /// Current function polynomial order
      int order;

      /// Number of vector components
      unsigned char num_components;

      /// With changed sub-element mapping, or an element, or anything else there comes the need for a change of the current values.
      /// This invalidates the current values.
      /// See values_valid
      void invalidate_values();

      /// List of available quadratures
      Quad2D* quads[H2D_MAX_QUADRATURES];
      /// Active quadrature (index into 'quads')
      int cur_quad;

      /// Index to mask table
      static int idx2mask[H2D_NUM_FUNCTION_VALUES][2];

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename T> friend class RefinementSelectors::HcurlProjBasedSelector;
      template<typename T> friend class Adapt;

      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      friend class CurvMap;

      template<typename T> friend class Func;
      template<typename T> friend class Filter;
      template<typename T> friend class SimpleFilter;
      template<typename T> friend class DXDYFilter;
      friend class ComplexFilter;
      friend class VonMisesFilter;
      friend HERMES_API Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
      template<typename T> friend HERMES_API Func<T>* init_fn(MeshFunction<T>*fu, const int order);
    };
  }
}
#endif
