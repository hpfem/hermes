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

#ifndef __H2D_REFINEMENT_SELECTOR_H
#define __H2D_REFINEMENT_SELECTOR_H

#include "element_to_refine.h"
#include "../mesh/mesh.h"
#include "adapt/error_calculator.h"

/** \defgroup g_selectors Refinement Selectors
*  \brief Refinement selectors allows to select a refinement
*  according to an error of a candidate.
*
*  The error is calculate by comparing a candidate projected
*  to a reference solution with the reference solution. All selectors
*  has to be derived from the class Selector. An instance of
*  a selector should be created outside the adaptivity loop
*  in order to save initialization time.
*
*  Currently available selectors:
*  - static selectors: Selectors that does not support
*    a mechanism for selecting an optimal refinement from
*    a list of refinements.
*    -# HOnlySelector
*    -# POnlySelector
*  - optimum selectors: Selectors that are able to select
*    a refinement from a provided list of possibilities
*    which minimizes the error. All optimum selectors have
*    to be derived from the class OptimumSelector.
*    -# H1ProjBasedSelector
*    -# L2ProjBasedSelector
*    \if H2D_COMPLEX # HcurlProjBasedSelector \endif
*/
namespace Hermes
{
  namespace Hermes2D
  {
    /// Namespace which encapsulates all refinement selectors. \ingroup g_selectors
    namespace RefinementSelectors {
      /// A parent of all refinement selectors. Abstract class. \ingroup g_selectors
      /** All refinement selectors have to derive from this class or its children.
      *  The interface of the class provides methods for:
      *  - selecting a refinement based on a reference soution,
      *  - updating orders of a mesh shared among components. */
      template<typename Scalar>
      class HERMES_API Selector : public Hermes::Mixins::Loggable, public Hermes::Mixins::TimeMeasurable
      {
      public:
        virtual ~Selector() {};
          /// Selects a refinement.
          /** This methods has to be implemented.
          *  \param[in] element An element which is being refined.
          *  \param[in] quad_order An encoded order of the element.
          *  \param[in] rsln A reference solution which is used to select a refinement.
          *  \param[out] refinement A selected refinement. It contains a valid contents if and only if the method returns true.
          *  \return True if a refinement was proposed. False if the selector is unable to select a refinement or it suggest that the element should not be refined. */
          virtual bool select_refinement(Element* element, int quad_order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement) = 0;
          
      protected:
        const int min_order; ///< A minimum allowed order.
        const int max_order; ///< A maximum allowed order.

        /// Constructor
        /** \param[in] max_order A maximum order used by this selector. If it is ::H2DRS_DEFAULT_ORDER, a maximum supported order is used. */
        Selector(int min_order = 1, int max_order = H2DRS_DEFAULT_ORDER) : min_order(min_order), max_order(max_order) {};

        template<typename T> friend class Adapt;
        template<typename T> friend class KellyTypeAdapt;
      };

      /// A selector that selects H-refinements only. \ingroup g_selectors
      template<typename Scalar>
      class HERMES_API HOnlySelector : public Selector<Scalar> {
      public:
        /// Constructor.
        HOnlySelector() : Selector<Scalar>() {};
      protected:
        /// Selects a refinement.
        /** Selects a H-refienements. For details, see Selector::select_refinement. */
        virtual bool select_refinement(Element* element, int quad_order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement);
        
        template<typename T> friend class Adapt;
        template<typename T> friend class KellyTypeAdapt;
      };

      /// A selector that increases order (i.e., it selects P-refinements only). \ingroup g_selectors
      template<typename Scalar>
      class HERMES_API POnlySelector : public Selector<Scalar> {
        const int order_h_inc; ///< Increase along the horizontal direction in a quadrilateral or increase of an order in a triangle.
        const int order_v_inc; ///< Increase along the vertical direction in a quadrilateral.
      public:
        /// Constructor.
        /** \param[in] max_order A maximum order used by this selector. If it is ::H2DRS_DEFAULT_ORDER, a maximum supported order is used.
        *  \param[in] order_h_inc An increase of the horizontal order in a quadrilateral and an order in a triangle. The increase has to be greater or equal to 0.
        *  \param[in] order_v_inc An increase of the vertical order in a quadrilateral. The increase has to be greater or equal to 0. */
        POnlySelector(int max_order, int order_h_inc, int order_v_inc);

      protected:
        /// Selects a refinement.
        /** Increases an order ising POnlySelector::order_h_inc and POnlySelector::order_v_inc. Fails if the order cannot be increased due to the maximum order. For details, see Selector::select_refinement. */
        virtual bool select_refinement(Element* element, int quad_order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement);

        template<typename T> friend class Adapt;
        template<typename T> friend class KellyTypeAdapt;
      };
    }
  }
}
#endif
