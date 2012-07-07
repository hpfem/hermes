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

#ifndef __H2D_REFINEMENT_ORDER_PERMUTATOR_H
#define __H2D_REFINEMENT_ORDER_PERMUTATOR_H

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors {
      /// Hermes::Order permutator. Generates all permutations of orders from a set defined by a range of orders. \ingroup g_selectors
      /** Used to generate all possible generation of horizontal and vertical orders bettween
      *  the starting order and the ending order. Both the starting order and the ending order are included.
      *  \par Suggested Use
      *  \code
      OrderPermutator order_perm(start_order, end_order, false);
      do {
      int order_h = order_perm.get_order_h();
      } while (order_perm.next());
      *  \endcode */
      class HERMES_API OrderPermutator {
      protected:
        int order_h; ///< The current horizontal order.
        int order_v; ///< The current verical order.
        int start_order_h; ///< The starting horizontal order.
        int start_order_v; ///< The starting vertical order.
        int end_order_h; ///< The ending horizontal order.
        int end_order_v; ///< The ending vertical order.
        bool iso_p; ///< True if orders is incresed in both the horizontal order and the vertical order are increased simultaneously.
        int* tgt_quad_order; ///< A pointer to which a current order is stored in encoded form. Ignored if NULL.

      public:
        /// Constructor.
        /** \param[in] start_quad_order The starting order in an encoded form.
        *  \param[in] end_quad_order The ending order in an encoded form.
        *  \param[in] iso_p True if both orders (i.e. horizontal and vertical) should be increased simultaneously.
        *  \param[in] tgt_quad_order A point to a location to which a current orders (i.e. horizontal and vertical) are stored in an encoded form. Ignored if NULL. */
        OrderPermutator (int start_quad_order = 0, int end_quad_order = 0, bool iso_p = false, int* tgt_quad_order = NULL);

        /// Moves to the next permutation of orders.
        /** \return True if there is a next permutation of orders. */
        bool next();

        /// Resets permutator to the starting order.
        void reset();

        /// Returns the current horizontal order.
        /** \return The current horizontal order. */
        int get_order_h() const;

        /// Returns the current vertical order.
        /** \return The current vertical order. */
        int get_order_v() const;

        /// Returns the current order in an encoded form.
        /** \return The current order in an encoded form. */
        int get_quad_order() const;

        /// Returns the starting order in an encoded form.
        /** \return The starting order in an encoded form. */
        int get_start_quad_order() const;

        /// Returns the ending order in an encoded form.
        /** \return The ending order in an encoded form. */
        int get_end_quad_order() const;
      };
    }
  }
}
#endif