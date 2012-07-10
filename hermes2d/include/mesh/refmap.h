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
#include "../quadrature/quad_all.h"
#include "shapeset/shapeset_h1_all.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class Element;
    namespace Views{
      class Orderizer;
      class Linearizer;
      class Vectorizer;
    };

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

      /// Returns the 1D quadrature for use in surface integrals.
      const Quad1D* get_quad_1d() const;

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
      void untransform(Element* e, double x, double y, double& xi1, double& xi2);

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
      bool is_jacobian_const() const;

      /// If the jacobian of the reference map is constant, this is the fast
      /// way to obtain it.
      double get_const_jacobian() const;

      /// Returns the jacobian of the reference map precalculated at the integration
      /// points of the specified order. Intended for non-constant jacobian elements.
      double* get_jacobian(int order);

      /// Returns the increase in the integration order due to the reference map.
      int get_inv_ref_order() const;

      H1ShapesetJacobi ref_map_shapeset;
      PrecalcShapeset ref_map_pss;
    private:
      /// If the reference map is constant, this is the fast way to obtain
      /// its inverse matrix.
      double2x2* get_const_inv_ref_map();

      /// Returns the inverse matrices of the reference map precalculated at the
      /// integration points of the specified order. Intended for non-constant
      /// jacobian elements.
      double2x2* get_inv_ref_map(int order);

      /// Returns coefficients for weak forms with second derivatives.
      double3x2* get_second_ref_map(int order);

      /// Calculates the inverse Jacobi matrix of reference map at a particular point (xi1, xi2).
      void inv_ref_map_at_point(double xi1, double xi2, double& x, double& y, double2x2& m);

      /// Calculates the second reference map at a particular point (xi1, xi2).
      void second_ref_map_at_point(double xi1, double xi2, double& x, double& y, double3x2& mm);

      /// See Transformable::push_transform()
      virtual void push_transform(int son);

      /// See Transformable::pop_transform()
      virtual void pop_transform();

      /// Frees all data associated with the instance.
      void free();

      /// For internal use only.
      void force_transform(uint64_t sub_idx, Trf* ctm)
      {
        this->sub_idx = sub_idx;
        stack[top] = *ctm;
        this->ctm = stack + top;
        update_cur_node();
        if(is_const) calc_const_inv_ref_map();
      }

      Quad2D* quad_2d;

      int num_tables;

      bool is_const;

      int inv_ref_order;

      double const_jacobian;

      double2x2 const_inv_ref_map;

      static const int H2D_MAX_TABLES = g_max_quad + 1 + 4 * g_max_quad + 4;

      /// This structure represents one complete piece of information about the reference mapping
      /// taking into account the sub-element mapping.
      struct Node
      {
        double* jacobian[H2D_MAX_TABLES];
        double2x2* inv_ref_map[H2D_MAX_TABLES];
        double3x2* second_ref_map[H2D_MAX_TABLES];
        double* phys_x[H2D_MAX_TABLES];
        double* phys_y[H2D_MAX_TABLES];
        double3* tan[4];
      };

      /// Table of RefMap::Nodes, indexed by a sub-element mapping.
      std::map<uint64_t, Node*> nodes;

      Node* cur_node;

      Node* overflow;

      void update_cur_node()
      {
        Node* updated_node = new Node;

        if(sub_idx > H2D_MAX_IDX) {
          delete updated_node;
          cur_node = handle_overflow();
        }
        else {
          if(nodes.insert(std::make_pair(sub_idx, updated_node)).second == false)
            /// The value had already existed.
            delete updated_node;
          else
            /// The value had not existed.
            init_node(updated_node);
          cur_node = nodes[sub_idx];
        }
      }

      void calc_inv_ref_map(int order);

      /// Quickly calculates the (hard-coded) reference mapping for elements with constant jacobians
      /// (ie., linear triangles and linear parallelogram quads).
      void calc_const_inv_ref_map();

      void calc_second_ref_map(int order);

      bool is_parallelogram();

      void calc_phys_x(int order);

      void calc_phys_y(int order);

      void calc_tangent(int edge, int eo);

      /// Finds the necessary quadrature degree needed to integrate the inverse reference mapping
      /// matrix alone. This is added to the total integration order in weak form itegrals.
      int calc_inv_ref_order();

      void init_node(Node* pp);

      void free_node(Node* node);

      Node* handle_overflow();

      Quad1DStd quad_1d;

      int indices[70];

      int nc;

      double2* coeffs;

      double2  lin_coeffs[4];
      template<typename T> friend class MeshFunction;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class Solution;
      template<typename T> friend class ExactSolution;
      template<typename T> friend class ExactSolutionScalar;
      template<typename T> friend class ExactSolutionVector;
      template<typename T> friend class Adapt;
      template<typename T> friend class KellyTypeAdapt;
      friend class Views::Orderizer;
      friend class Views::Vectorizer;
      friend class Views::Linearizer;
      template<typename T> friend class Global;
      friend class VonMisesFilter;
      template<typename T> friend class Func;
      template<typename T> friend class Geom;
      friend Geom<double>* init_geom_vol(RefMap *rm, const int order);
      friend Geom<double>* init_geom_surf(RefMap *rm, SurfPos* surf_pos, const int order);
      friend Func<double>* init_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
    };
  }
}
#endif