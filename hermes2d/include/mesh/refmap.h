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

#ifndef __H2D_REFMAP_H
#define __H2D_REFMAP_H

#include "../global.h"
#include "../shapeset/precalc.h"
#include "../mesh/mesh.h"
#include "../quadrature/quad_all.h"
#include "shapeset/shapeset_h1_all.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class Element;
    class Mesh;

    /// @ingroup meshFunctions
    /// \brief Represents the reference mapping.
    ///
    /// RefMap represents the mapping from the reference to the physical element.
    /// Its main task is to provide both the (inverse) reference mapping matrix and
    /// its jacobian, precalculated in quadrature points. Other functions include
    /// the calculation of integration points positions in the physical domain and
    /// the calculation of edge tangents in 1D integration points.
    ///
    class HERMES_API RefMap : public Transformable
    {
    public:
      RefMap();

      ~RefMap();

      /// Sets the quadrature points in which the reference map will be evaluated.
      /// \param quad_2d[in] The quadrature points.
      void set_quad_2d(Quad2D* quad_2d);

      /// Returns the current quadrature points.
      Quad2D* get_quad_2d() const;

      /// Initializes the reference map for the specified element.
      /// Must be called prior to using all other functions in the class.
      virtual void set_active_element(Element* e);

      /// Returns the triples[x, y, norm] of the tangent to the specified (possibly
      /// curved) edge at the 1D integration points along the edge. The maximum
      /// 1D quadrature rule is used by default, but the user may specify his own
      /// order. In this case, the edge pseudo-order is expected (as returned by
      /// Quad2D::get_edge_points).
      double3* get_tangent(int edge, int order = -1);

      /// Transforms physical coordinates x, y from the element e back to the reference domain.
      /// If the point (x, y) does not lie in e, then (xi1, xi2) will not lie in the reference domain.
      static void untransform(Element* e, double x, double y, double& xi1, double& xi2);

      /// Returns the element pointer located at physical coordinates x, y.
      /// \param[in] x Physical x-coordinate.
      /// \param[in] y Physical y-coordinate.
      /// \param[in] x_reference Optional parameter, in which the x-coordinate of x in the reference domain will be returned.
      /// \param[in] y_reference Optional parameter, in which the y-coordinate of y in the reference domain will be returned.
      static Element* element_on_physical_coordinates(bool use_MeshHashGrid, MeshSharedPtr mesh, double x, double y, double* x_reference = nullptr, double* y_reference = nullptr);

      /// Find out if the coordinatex [x,y] lie in the element e.
      /// \param[in] x Physical x-coordinate.
      /// \param[in] y Physical y-coordinate.
      /// \param[in] x_reference Optional parameter, in which the x-coordinate of x in the reference domain will be returned.
      /// \param[in] y_reference Optional parameter, in which the y-coordinate of y in the reference domain will be returned.
      static bool is_element_on_physical_coordinates(Element* e, double x, double y, double* x_reference = nullptr, double* y_reference = nullptr);

      /// Returns the x-coordinates of the integration points transformed to the
      /// physical domain of the element. Intended for integrals containing spatial
      /// variables.
      double* get_phys_x(int order);

      /// Returns he y-coordinates of the integration points transformed to the
      /// physical domain of the element. Intended for integrals containing spatial
      /// variables.
      double* get_phys_y(int order);

      /// Returns true if the jacobian of the reference map is constant (which
      /// is the case for non-curvilinear triangular elements), false otherwise.
      inline bool is_jacobian_const() const
      {
        return is_const;
      }

      /// Returns the increase in the integration order due to the reference map.
      inline int get_inv_ref_order() const
      {
        return inv_ref_order;
      }

      /// If the jacobian of the reference map is constant, this is the fast
      /// way to obtain it.
      inline double get_const_jacobian() const
      {
        return const_jacobian;
      }

      /// If the reference map is constant, this is the fast way to obtain
      /// its inverse matrix.
      inline double2x2* get_const_inv_ref_map()
      {
        return &const_inv_ref_map;
      }

      /// Returns the jacobian of the reference map precalculated at the integration
      /// points of the specified order. Intended for non-constant jacobian elements.
      inline double* get_jacobian(int order)
      {
        if (this->is_const)
          throw Hermes::Exceptions::Exception("RefMap::get_jacobian() called with a const jacobian.");
        if (order != this->jacobian_calculated)
          this->calc_inv_ref_map(order);
        return this->jacobian;
      }

      /// Returns the inverse matrices of the reference map precalculated at the
      /// integration points of the specified order. Intended for non-constant
      /// jacobian elements.
      inline double2x2* get_inv_ref_map(int order)
      {
        if (this->is_const)
          throw Hermes::Exceptions::Exception("RefMap::get_inv_ref_map() called with a const jacobian.");
        if (order != this->inv_ref_map_calculated)
          this->calc_inv_ref_map(order);
        return this->inv_ref_map;
      }

      /// Calculates the inverse Jacobi matrix of reference map at a particular point (xi1, xi2).
      void inv_ref_map_at_point(double xi1, double xi2, double& x, double& y, double2x2& m);

#ifdef H2D_USE_SECOND_DERIVATIVES
      /// Returns coefficients for weak forms with second derivatives.
      double3x2* get_second_ref_map(int order);

      /// Calculates the second reference map at a particular point (xi1, xi2).
      void second_ref_map_at_point(double xi1, double xi2, double& x, double& y, double3x2& mm);
#endif

      /// See Transformable::push_transform()
      virtual void push_transform(int son);

      /// See Transformable::pop_transform()
      virtual void pop_transform();

      /// For internal use only.
      void force_transform(uint64_t sub_idx, Trf* ctm);

      static bool is_parallelogram(Element* e);

      void set_element_iro_cache(Element* element);

    private:
      /// re-init the storage
      void reinit_storage();

      H1ShapesetJacobi ref_map_shapeset;
      PrecalcShapesetAssembling ref_map_pss;

      /// Constant reference mapping.
      bool is_const;

      /// For constant ref. map.
      double const_jacobian;
      double2x2 const_inv_ref_map;
      double3x2 const_second_ref_map;
      double const_direct_ref_map[2][2];

      /// For non-constant ref. map.
      double jacobian[H2D_MAX_INTEGRATION_POINTS_COUNT];
      int jacobian_calculated;
      double2x2 inv_ref_map[H2D_MAX_INTEGRATION_POINTS_COUNT];
      int inv_ref_map_calculated;
      double3x2 second_ref_map[H2D_MAX_INTEGRATION_POINTS_COUNT];
      int second_ref_map_calculated;
      double direct_ref_map[2][2][H2D_MAX_INTEGRATION_POINTS_COUNT];
      int direct_ref_map_calculated;
      int inv_ref_order;

      /// Data
      double phys_x[H2D_MAX_INTEGRATION_POINTS_COUNT];
      int phys_x_calculated;
      double phys_y[H2D_MAX_INTEGRATION_POINTS_COUNT];
      int phys_y_calculated;
      double3 tan[H2D_MAX_NUMBER_EDGES][H2D_MAX_INTEGRATION_POINTS_COUNT];
      int tan_calculated[H2D_MAX_NUMBER_EDGES];

      Quad2D* quad_2d;

      void calc_inv_ref_map(int order);

      /// Quickly calculates the (hard-coded) reference mapping for elements with constant jacobians
      /// (ie., linear triangles and linear parallelogram quads).
      void calc_const_inv_ref_map();

#ifdef H2D_USE_SECOND_DERIVATIVES
      void calc_second_ref_map(int order);
#endif
      void calc_phys_x(int order);

      void calc_phys_y(int order);

      void calc_tangent(int edge, int eo);

      /// Finds the necessary quadrature degree needed to integrate the inverse reference mapping
      /// matrix alone. This is added to the total integration order in weak form itegrals.
      int calc_inv_ref_order();

      unsigned short indices[70];

      int nc;

      double2* coeffs;

      double2  lin_coeffs[H2D_MAX_NUMBER_EDGES];
    };
  }
}
#endif