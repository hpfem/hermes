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

#ifndef _MSC_VER
#include "../mesh/refinement_type.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class ElementToRefine;
    struct Element;
    template<typename Scalar> class Solution;
  }
}
#else
#include "../mesh/element_to_refine.h"
#endif
#include "../mesh/mesh.h"

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
#define H2DRS_DEFAULT_ORDER -1 ///< A default order. Used to indicate an unkonwn order or a maximum support order.  \ingroup g_selectors
#define H2DRS_MAX_ORDER 9 ///< A maximum order suported by refinement selectors. \ingroup g_selectors

    /// Namespace which encapsulates all refinement selectors. \ingroup g_selectors
    namespace RefinementSelectors {
      /// A parent of all refinement selectors. Abstract class. \ingroup g_selectors
      /** All refinement selectors have to derive from this class or its children.
      *  The interface of the class provides methods for:
      *  - selecting a refinement based on a reference soution,
      *  - updating orders of a mesh shared among components. */
      template<typename Scalar>
      class HERMES_API Selector : public Hermes::Mixins::Loggable
      {
      protected:
        const int max_order; ///< A maximum allowed order.
        /// Constructor
        /** \param[in] max_order A maximum order used by this selector. If it is ::H2DRS_DEFAULT_ORDER, a maximum supported order is used. */
        Selector(int max_order = H2DRS_DEFAULT_ORDER) : max_order(max_order) {};

        /// Cloning for paralelism.
        virtual Selector<Scalar>* clone() = 0;

        /// Selects a refinement.
        /** This methods has to be implemented.
        *  \param[in] element An element which is being refined.
        *  \param[in] quad_order An encoded order of the element.
        *  \param[in] rsln A reference solution which is used to select a refinement.
        *  \param[out] refinement A selected refinement. It contains a valid contents if and only if the method returns true.
        *  \return True if a refinement was proposed. False if the selector is unable to select a refinement or it suggest that the element should not be refined. */
        virtual bool select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement) = 0;

        /// Generates orders of elements which will be created due to a proposed refinement in another component that shares the same a mesh.
        /** \param[in] element An element which is about the be refined.
        *  \param[in] orig_quad_order An encoded order of the element.
        *  \param[in] refinement A refinement of the element in the mesh. Possible values are defined by the enum RefinementType.
        *  \param[out] tgt_quad_orders Generated encoded orders.
        *  \param[in] suggested_quad_orders Suggested encoded orders. If not NULL, the method should copy them to the output. If NULL, the method have to calculate orders. */
        virtual void generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders) = 0;

        template<typename T> friend class Adapt;
        template<typename T> friend class KellyTypeAdapt;
      };

      /// A selector that selects H-refinements only. \ingroup g_selectors
      template<typename Scalar>
      class HERMES_API HOnlySelector : public Selector<Scalar> {
      public:
        /// Constructor.
        HOnlySelector() : Selector<Scalar>() {};

        /// Cloning for paralelism.
        virtual Selector<Scalar>* clone();

      protected:
        /// Selects a refinement.
        /** Selects a H-refienements. For details, see Selector::select_refinement. */
        virtual bool select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement);

        /// Generates orders of elements which will be created due to a proposed refinement in another component that shares the same a mesh.
        /** If a parameter suggested_quad_orders is NULL, the method uses an encoded order in orig_quad_order.
        *  For details, see Selector::generate_shared_mesh_orders. */
        virtual void generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders);
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

        /// Cloning for paralelism.
        virtual Selector<Scalar>* clone();

      protected:
        /// Selects a refinement.
        /** Increases an order ising POnlySelector::order_h_inc and POnlySelector::order_v_inc. Fails if the order cannot be increased due to the maximum order. For details, see Selector::select_refinement. */
        virtual bool select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement);

        /// Generates orders of elements which will be created due to a proposed refinement in another component that shares the same a mesh.
        /** If a parameter suggested_quad_orders is NULL, the method uses an encoded order in orig_quad_order.
        *  For details, see Selector::generate_shared_mesh_orders. */
        virtual void generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders);
        template<typename T> friend class Adapt;
        template<typename T> friend class KellyTypeAdapt;
      };
    }
  }
}
#endif