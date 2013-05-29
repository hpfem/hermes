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

#ifndef __H2D_ELEMENT_TO_REFINE_H
#define __H2D_ELEMENT_TO_REFINE_H

#include "global.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors{
      template<typename Scalar> class Selector;
      template<typename Scalar> class HOnlySelector;
      template<typename Scalar> class POnlySelector;
      template<typename Scalar> class OptimumSelector;
    };

    /// A refinement record. \ingroup g_adapt
    /** Except the attribute ElementToRefine::q, the class is able to dump its content to a stringstream
    *  in a human readable form, see operator<<(std::ostream& stream, const ElementToRefine& elem_ref). */
    class HERMES_API ElementToRefine {
    public:
      /// Constructor. Creates an invalid refinement.
      ElementToRefine();

      /// Constructor.
      /** \param[in] id An ID of the element.
      *  \param[in] comp An index of a component. */
      ElementToRefine(int id, int comp);

      /// Copy-contructor.
      ElementToRefine(const ElementToRefine &orig);

      /// Assignment operator.
      ElementToRefine& operator=(const ElementToRefine& orig);
    
      /// An ID of the element.
      int id;
      /// An index of the component.
      int comp;
      /// Proposed refinement. Possible values are defined in the enum ::RefinementType.
      int split;
      /// Encoded orders of sons.
      int refinement_polynomial_order[H2D_MAX_ELEMENT_SONS];
      /// Encoded orders of the best refinement of a certaint type.
      /// Indexed by enum RefinementType.
      int best_refinement_polynomial_order_type[4][H2D_MAX_ELEMENT_SONS]; 
      /// Error of the selected candidate.
      double errors[H2D_MAX_ELEMENT_SONS];

      /// Returns a number of sons.
      /** \return A number of sons of a given refinement. */
      int get_num_sons() const;

      /// Copies array of orders.
      /** The length of the array is defubed by ::H2D_MAX_ELEMENT_SONS.
      *  \param[in] dest A destination array.
      *  \param[in] src A source arrapy. */
      static void copy_orders(int* dest, const int* src);
      static void copy_errors(double* dest, const double* src);
      template<typename T> friend class Adapt;
      template<typename T> friend class ErrorCalculator;
      template<typename T> friend class RefinementSelectors::Selector;
      template<typename T> friend class RefinementSelectors::HOnlySelector;
      template<typename T> friend class RefinementSelectors::POnlySelector;
      template<typename T> friend class RefinementSelectors::OptimumSelector;
    };
  }
}
#endif