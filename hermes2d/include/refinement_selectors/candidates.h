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

#ifndef __H2D_CANDIDATES_H
#define __H2D_CANDIDATES_H

#include <ostream>
#include "global.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Possible refinements of an element.
    enum RefinementType {
      H2D_REFINEMENT_P = 0, ///< P-refinement.
      H2D_REFINEMENT_H = 1, ///< H-refinement.
      H2D_REFINEMENT_ANISO_H = 2, ///< ANISO-refienement. The element is split along the horizontal axis. Quadrilaterals only.
      H2D_REFINEMENT_ANISO_V = 3 ///< ANISO-refienement. The element is split along the vertical axis. Quadrilaterals only.
    };

    /// Retuns true if a given refinement is an ANISO-refinement.
    /** \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
     *  \return True of a given refinement is an ANISO-refinement. */
    extern HERMES_API bool is_refin_aniso(const int refin_type);

    /// Returns a maximum number of sons that will be generated if a given refinement is applied.
    /** \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
     *  \return A number of possible sons. In a case of P-refinement, the function returns 1 even thought this refinement yields just a change in orders. */
    extern HERMES_API int get_refin_sons(const int refin_type);

    /// Returns a string representation of the refinement.
    /** Used for debugging and event logging purposes.
     *  \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
     *  \return A string representation of a given refinement. */
    extern HERMES_API const std::string get_refin_str(const int refin_type);

    namespace RefinementSelectors
    {
      /// Predefined list of candidates. \ingroup g_selectors
      enum CandList {
        H2D_NONE,  ///< No adaptivity. (Used only in modules.)
        H2D_P_ISO, ///< P-candidates only. Hermes::Orders are modified uniformly.
        H2D_P_ANISO, ///< P-candidates only. Hermes::Orders are modified non-uniformly.
        H2D_H_ISO, ///< H-candidates only. Hermes::Orders are not modified.
        H2D_H_ANISO, ///< H- and ANISO-candidates only. Hermes::Orders are not modified.
        H2D_HP_ISO, ///< H- and P-candidates only. Hermes::Orders are modified uniformly.
        H2D_HP_ANISO_H, ///< H-, ANISO- and P-candidates. Hermes::Orders are modified uniformly.
        H2D_HP_ANISO_P, ///< H- and P-candidates only. Hermes::Orders are modified non-uniformly.
        H2D_HP_ANISO ///< H-, ANISO- and P-candidates. Hermes::Orders are modified non-uniformly.
      };

      /// Returns a string representation of a predefined candidate list. \ingroup g_selectors
      /** Used for debugging and output purposes.
      *  \param cand_list A predefined list of candidates.
      *  \return A string representation of the enum value. */
      extern HERMES_API const char* get_cand_list_str(const CandList cand_list);

      /// Returns true if a predefined candidate list may contain candidates that are HP. \ingroup g_selectors
      /** \param cand_list A predefined list of candidates.
      *  \return True if a predefined candidate list may contain candidates that are HP. */
      extern HERMES_API bool is_hp(const CandList cand_list);

      /// Returns true if a predefined candidate list may contain candidates that increase P. \ingroup g_selectors
      /** \param cand_list A predefined list of candidates.
      *  \return True if a predefined candidate list may contain candidates that increase P. */
      extern HERMES_API bool is_p(const CandList cand_list);

      /// Returns true if a predefined candidate list may contain candidates with an anisotropic change of orders. \ingroup g_selectors
      /** \param cand_list A predefined list of candidates.
      *  \return True if a predefined candidate list may contain candidates with an anisotropic change of orders. */
      extern HERMES_API bool is_p_aniso(const CandList cand_list);

      /// A candidate.
      class HERMES_API Cand 
      {
      public:
        double error; ///< Error of this candidate's sons.
        double errors[H2D_MAX_ELEMENT_SONS]; ///< Error of this candidate's sons.
        int dofs;  ///< An estimated number of DOFs.
        int split; ///< A refinement, see the enum RefinementType.
        int p[H2D_MAX_ELEMENT_SONS]; ///< Encoded orders of sons, see ::H2D_MAKE_QUAD_ORDER. In a case of a triangle, the vertical order is equal to the horizontal one.
        double score; ///< A score of a candidate: the higher the better. If zero, the score is not valid and a candidate should be ignored. Evaluated in OptimumSelector::select_best_candidate.

        /// Constructor.
        /** \param[in] split A refinement, see the enum RefinementTypes.
        *  \param[in] order_elems Encoded orders for all element of candidate. If triangle, a vertical order has to be equal to the horizontal one. Unused elements of the array can be ignored. */
        Cand(const int split, const int order_elems[H2D_MAX_ELEMENT_SONS]);

        /// Constructor.
        /** \param[in] split A refinement, see the enum RefinementTypes.
        *  \param[in] order_elem0 Encoded order of the first element of the candidate. If triangle, a vertical order has to be equal to the horizontal one.
        *  \param[in] order_elem1 Encoded order of the second element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one.
        *  \param[in] order_elem2 Encoded order of the third element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one.
        *  \param[in] order_elem3 Encoded order of the fourth element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one. */
        Cand(const int split, const int order_elem0, const int order_elem1 = 0, const int order_elem2 = 0, const int order_elem3 = 0);

        /// Returns a number of elements of a candidate.
        /** \return A number of elements of a candidate. */
        int get_num_elems() const;
      };
    }
  }
}

#endif