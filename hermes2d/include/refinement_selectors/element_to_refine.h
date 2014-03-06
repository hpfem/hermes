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
    /// Possible refinements of an element.
    enum RefinementType {
      H2D_REFINEMENT_P = 0, ///< P-refinement.
      H2D_REFINEMENT_H = 1, ///< H-refinement.
      H2D_REFINEMENT_H_ANISO_H = 2, ///< ANISO-refienement. The element is split along the horizontal axis. Quadrilaterals only.
      H2D_REFINEMENT_H_ANISO_V = 3 ///< ANISO-refienement. The element is split along the vertical axis. Quadrilaterals only.
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
      ElementToRefine(int id, unsigned short comp);

      /// Copy-contructor.
      ElementToRefine(const ElementToRefine &orig);

      /// Assignment operator.
      ElementToRefine& operator=(const ElementToRefine& orig);
    
      /// Validity info.
      bool valid;

      /// An ID of the element.
      int id;
      /// An index of the component.
      unsigned short comp;
      /// Proposed refinement. Possible values are defined in the enum ::RefinementType.
      RefinementType split;
      /// Encoded orders of sons.
      unsigned short refinement_polynomial_order[H2D_MAX_ELEMENT_SONS];
      /// Encoded orders of the best refinement of a certaint type.
      /// Indexed by enum RefinementType.
      unsigned short best_refinement_polynomial_order_type[4][H2D_MAX_ELEMENT_SONS]; 
      /// Error of the selected candidate.
      double errors[H2D_MAX_ELEMENT_SONS];

      /// Returns a number of sons.
      /** \return A number of sons of a given refinement. */
      unsigned short get_num_sons() const;

      /// Copies array of orders.
      /** The length of the array is defubed by ::H2D_MAX_ELEMENT_SONS.
      *  \param[in] dest A destination array.
      *  \param[in] src A source arrapy. */
      static void copy_orders(unsigned short* dest, const unsigned short* src);
      static void copy_errors(double* dest, const double* src);
    private:
      /// This array is internal.
      bool refinement_polynomial_order_changed[H2D_MAX_ELEMENT_SONS];
      template<typename T> friend class Adapt;
      template<typename T, typename S> friend class AdaptSolver;
    };
  }
}
#endif