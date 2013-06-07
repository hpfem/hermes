// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file global.h
\brief Common definitions for Hermes2D.
\todo Put more common H2D stuff here.
*/
#ifndef __H2D_COMMON_H_
#define __H2D_COMMON_H_

#include "hermes_common.h"
#include "config.h"

/// Macros.
#define H2D_MAX_ELEMENT_SONS 4 ///< A maximum number of sons of an element.
#define H2D_MAX_NUMBER_EDGES 4 ///< A maximum number of edges of an element.
#define H2D_MAX_NUMBER_VERTICES 4 ///< A maximum number of vertices of an element.

/// Internal.
#define H2D_NUM_MODES 2 ///< A number of modes, see enum ElementMode2D.
#define H2D_SOLUTION_ELEMENT_CACHE_SIZE 4 ///< An internal parameter.
#define H2D_MAX_NODE_ID 10000000
#define H2D_MAX_SOLUTION_COMPONENTS 2
static const std::string H2D_DG_INNER_EDGE = "-1234567";

#define HERMES_ONE NULL
#define HERMES_DEFAULT_FUNCTION NULL
#define HERMES_DEFAULT_SPLINE NULL

/// A total number of valid transformation of a triangle to a sub-domain.
static const int H2D_TRF_TRI_NUM = 4;
/// A total number of valid transformation of a quad to a sub-domain.
static const int H2D_TRF_QUAD_NUM = 8;
/// A total number of transformations.
static const int H2D_TRF_NUM = (H2D_TRF_QUAD_NUM + 1);
/// An index of identity transformation.
static const int H2D_TRF_IDENTITY = H2D_TRF_QUAD_NUM;

#define H2DRS_ASSUMED_MAX_CANDS 512 ///< An estimated maximum number of candidates. Used for purpose of reserving space. \internal \ingroup g_selectors

//TODO: find out why 20 used used, should'n be there 2*(H2DRS_MAX_ORDER+1)
#define H2DRS_INTR_GIP_ORDER 20 ///< An integration order used to integrate while evaluating a candidate. \internal \ingroup g_selectors
#define H2DRS_MAX_ORDER_INC 1 ///< Maximum increase of an order in candidates. \ingroup g_selectors

#define H2DRS_ORDER_ANY -1 ///< Any order. Used as a wildcard to indicate that a given order can by any valid order. \internal \ingroup g_selectors

# define H2DRS_DEFAULT_ERR_WEIGHT_H 2.0 ///< A default multiplicative coefficient of an error of a H-candidate. \ingroup g_selectors
# define H2DRS_DEFAULT_ERR_WEIGHT_P 1.0 ///< A default multiplicative coefficient of an error of a P-candidate. \ingroup g_selectors
# define H2DRS_DEFAULT_ERR_WEIGHT_ANISO 1.414214 ///< A default multiplicative coefficient of an error of a ANISO-candidate. \ingroup g_selectors

template<typename Scalar> class MeshFunctionSharedPtr;
template<typename Scalar> class MeshFunctionSharedPtr;

namespace Hermes
{
  /// Namespace containing definitions specific for Hermes2D.
  namespace Hermes2D
  {
    class RefMap;
    template<typename Scalar> class DiscreteProblem;
    template<typename Scalar> class Space;
    template<typename Scalar> class WeakForm;
    template<typename Scalar> class MeshFunction;
    template<typename Scalar> class Solution;
    class Quad2D;
    class Quad1DStd;
    class Quad2DStd;

    /// How many bits the encoded_order number takes.
    const int H2D_ORDER_BITS = 5;
    const int H2D_ORDER_MASK = (1 << H2D_ORDER_BITS) - 1;

    /// Macros for combining quad horizontal and vertical encoded_orders.
#define H2D_GET_H_ORDER(encoded_order) ((encoded_order) & H2D_ORDER_MASK)
#define H2D_GET_V_ORDER(encoded_order) ((encoded_order) >> H2D_ORDER_BITS)
#define H2D_MAKE_QUAD_ORDER(h_encoded_order, v_encoded_order) (((v_encoded_order) << H2D_ORDER_BITS) + (h_encoded_order))
#define H2D_MAKE_EDGE_ORDER(mode, edge, order) ((mode == HERMES_MODE_TRIANGLE || edge == 0 || edge == 2) ? H2D_GET_H_ORDER(order) : H2D_GET_V_ORDER(order))

    template<typename Scalar>
    double get_l2_norm(Vector<Scalar>* vec)
    {
      Scalar val = 0;
      for (unsigned int i = 0; i < vec->length(); i++)
      {
        Scalar inc = vec->get(i);
        val = val + inc*conj(inc);
      }
      return sqrt(std::abs(val));
    };

    template<typename Scalar>
    double get_l2_norm(Scalar* vec, int count)
    {
      Scalar val = 0;
      for (unsigned int i = 0; i < count; i++)
      {
        Scalar inc = vec[i];
        val = val + inc*conj(inc);
      }
      return sqrt(std::abs(val));
    }

#define H2DRS_DEFAULT_ORDER -1 ///< A default order. Used to indicate an unkonwn order or a maximum support order.  \ingroup g_selectors
#define H2DRS_MAX_ORDER 10 ///< A maximum order suported by refinement selectors. \ingroup g_selectors
#define H2D_NUM_SHAPES_SIZE 12 ///< A maximum order suported by refinement selectors. \ingroup g_selectors

    /// Projection norms.
    /// Used in projections and adaptivity.
    enum NormType
    {
      HERMES_L2_NORM,
      HERMES_H1_NORM,
      HERMES_H1_SEMINORM,
      HERMES_HCURL_NORM,
      HERMES_HDIV_NORM,
      HERMES_UNSET_NORM
    };

    enum ElementMode2D {
      HERMES_MODE_TRIANGLE = 0,
      HERMES_MODE_QUAD = 1
    };

    enum SpaceType {
      HERMES_H1_SPACE = 0,
      HERMES_HCURL_SPACE = 1,
      HERMES_HDIV_SPACE = 2,
      HERMES_L2_SPACE = 3,
      HERMES_UTILITY_L2_SPACES = 4,
      HERMES_INVALID_SPACE = -9999
    };

    /// Geometrical type of weak forms.
    enum GeomType
    {
      HERMES_PLANAR = 0,         // Planar problem.
      HERMES_AXISYM_X = 1,       // Axisymmetric problem where x-axis is the axis of symmetry.
      HERMES_AXISYM_Y = 2        // Axisymmetric problem where y-axis is the axis of symmetry.
    };

    /// Bilinear form symmetry flag, see WeakForm::add_matrix_form
    enum SymFlag
    {
      HERMES_ANTISYM = -1,
      HERMES_NONSYM = 0,
      HERMES_SYM = 1
    };
  }
}
#endif