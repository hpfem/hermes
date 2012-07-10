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

#define H2D_NUM_MODES 2 ///< A number of modes, see enum ElementMode2D.

#define H2D_MAX_NODE_ID 10000000

#define HERMES_ONE NULL
#define HERMES_DEFAULT_FUNCTION NULL
#define HERMES_DEFAULT_SPLINE NULL

#define H2D_MAX_COMPONENTS 10 ///< The maximum number of components in Hermes2D.

    /// Constant used by Adapt::calc_eror().
#define HERMES_TOTAL_ERROR_REL  0x00  ///< A flag which defines interpretation of the total error. \ingroup g_adapt
    ///  The total error is divided by the norm and therefore it should be in a range[0, 1].
    ///  \note Used by Adapt::calc_errors_internal().. This flag is mutually exclusive with ::H2D_TOTAL_ERROR_ABS.
#define HERMES_TOTAL_ERROR_ABS  0x01  ///< A flag which defines interpretation of the total error. \ingroup g_adapt
    ///  The total error is absolute, i.e., it is an integral over squares of differencies.
    ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::HERMES_TOTAL_ERROR_REL.
#define HERMES_ELEMENT_ERROR_REL 0x00 ///< A flag which defines interpretation of an error of an element. \ingroup g_adapt
    ///  An error of an element is a square of an error divided by a square of a norm of a corresponding component.
    ///  When norms of 2 components are very different (e.g. microwave heating), it can help.
    ///  Navier-stokes on different meshes work only when absolute error (see ::H2D_ELEMENT_ERROR_ABS) is used.
    ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::H2D_ELEMENT_ERROR_ABS.
#define HERMES_ELEMENT_ERROR_ABS 0x10 ///< A flag which defines interpretation of of an error of an element. \ingroup g_adapt
    ///  An error of an element is a square of an asolute error, i.e., it is an integral over squares of differencies.
    ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::HERMES_ELEMENT_ERROR_REL.

    /// A total number of valid transformation of a triangle to a sub-domain.
    static const int H2D_TRF_TRI_NUM = 4;
    /// A total number of valid transformation of a quad to a sub-domain.
    static const int H2D_TRF_QUAD_NUM = 8;
    /// A total number of transformations.
    static const int H2D_TRF_NUM = (H2D_TRF_QUAD_NUM + 1);
    /// An index of identity transformation.
    static const int H2D_TRF_IDENTITY = H2D_TRF_QUAD_NUM;
    
    /// Enabling second derivatives in weak forms. Turned on by default. Second
    /// derivatives are employed, among others, by stabilization methods for
    /// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

#define H2DRS_ASSUMED_MAX_CANDS 512 ///< An estimated maximum number of candidates. Used for purpose of reserving space. \internal \ingroup g_selectors

//TODO: find out why 20 used used, should'n be there 2*(H2DRS_MAX_ORDER+1)
#define H2DRS_INTR_GIP_ORDER 20 ///< An integration order used to integrate while evaluating a candidate. \internal \ingroup g_selectors
#define H2DRS_MAX_ORDER_INC 2 ///< Maximum increase of an order in candidates. \ingroup g_selectors

#define H2DRS_SCORE_DIFF_ZERO 1E-13 ///< A threshold of difference between scores. Anything below this values is considered zero. \internal \ingroup g_selectors

#define H2DRS_ORDER_ANY -1 ///< Any order. Used as a wildcard to indicate that a given order can by any valid order. \internal \ingroup g_selectors

# define H2DRS_DEFAULT_ERR_WEIGHT_H 2.0 ///< A default multiplicative coefficient of an error of a H-candidate. \ingroup g_selectors
# define H2DRS_DEFAULT_ERR_WEIGHT_P 1.0 ///< A default multiplicative coefficient of an error of a P-candidate. \ingroup g_selectors
# define H2DRS_DEFAULT_ERR_WEIGHT_ANISO 1.414214 ///< A default multiplicative coefficient of an error of a ANISO-candidate. \ingroup g_selectors

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

    /// Class for global functions.
    template<typename Scalar>
    class HERMES_API Global : public Hermes::Mixins::Loggable
    {
    public:
      Global() : Hermes::Mixins::Loggable() {};
      friend void warn_order();
    public:
      /// Error calculation in Hermes, useful for non-adaptive computations.
      // Note: coarse mesh sln has to be first, then
      // ref_sln (because the abs. error is divided
      // by the norm of the latter).
      static double calc_rel_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type);

      static double calc_abs_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type);

      static double calc_norm(MeshFunction<Scalar>* sln, int norm_type);

      /// Calculate norm of a (possibly vector-valued) solution.
      /// Take norm from spaces where these solutions belong.
      static double calc_norms(Hermes::vector<Solution<Scalar>*> slns);
      static double calc_abs_errors(Hermes::vector<Solution<Scalar>*> slns1, Hermes::vector<Solution<Scalar>*> slns2);
      static double calc_rel_errors(Hermes::vector<Solution<Scalar>*> slns1, Hermes::vector<Solution<Scalar>*> slns2);

      static double error_fn_l2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv);
      static double norm_fn_l2(MeshFunction<Scalar>* sln, RefMap* ru);

      static double error_fn_h1(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv);
      static double norm_fn_h1(MeshFunction<Scalar>* sln, RefMap* ru);

      static double error_fn_hc(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv);
      static double norm_fn_hc(MeshFunction<Scalar>* sln, RefMap* ru);

      static double error_fn_hcl2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv);
      static double norm_fn_hcl2(MeshFunction<Scalar>* sln, RefMap* ru);

      static double error_fn_hdiv(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv);
      static double norm_fn_hdiv(MeshFunction<Scalar>* sln, RefMap* ru);

      static double get_l2_norm(Vector<Scalar>* vec);
    };

    /// Projection norms.
    /// Used in projections and adaptivity.
    enum ProjNormType
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
  }
}
#endif